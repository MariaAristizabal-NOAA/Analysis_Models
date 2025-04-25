#%% User input
# forecasting cycle to be used

# User input
# Helene
#cycle = '2024092500'
#storm_num = '09'
#basin = 'al'
#storm_id = '09l'
#storm_name = 'Helene'

# Milton
cycle = '2024100800'
storm_num = '14'
basin = 'al'
storm_id = '14l'
storm_name = 'Milton'

exp_names = ['HFSB_oper','hafs_20241220_v2p1b_baseline','hafs_20250306_v2p1b_hb43']
exp_labels = ['HFSB_oper','HAFSv2.1B baseline','HAFSv2.1B final']
exp_colors = ['lime','olive','darkorange']
ocean = ['hycom','hycom','mom6']
hafs = ['hfsb','hfsb','hfsb']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']
markers = ['o','o','o']

'''
exp_names = ['HFSA_oper','hafs_20241220_v2p1a_baseline','hafs_20250210_v2p1a_ha30']
exp_labels = ['HFSA_oper','HAFSv2.1A baseline','HAFSv2.1A final']
exp_colors = ['purple','dodgerblue','#00c8c8']
ocean = ['mom6','mom6','mom6']
hafs = ['hfsa','hfsa','hfsa']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']
markers = ['o','o','o']
'''

file_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2024.nc'

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

############################################################
def get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output):
    subvar_obs = []
    subvar_model = []

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]
    var_obs = var_obss[oklo]

    okt_sort = np.argsort(timestamp_obs)
    timestamp_obs_sort = timestamp_obs[okt_sort]
    lon_obs_sort = lon_obs[okt_sort]
    lat_obs_sort = lat_obs[okt_sort]
    var_obs_sort = var_obs[okt_sort]

    target_time = []

    for x,file in enumerate(files_model):
        print(x,' ',file)
        model = xr.open_dataset(file)
        if file.split('.')[-1] == 'nc':
            t =  mdates.num2date(mdates.date2num(np.asarray(model[time_name][:])))[0]
            ti = t - timedelta(hours=int(freq_output/2))
            tii = mdates.date2num(ti)
            te = t + timedelta(hours=int(freq_output/2))
            tee = mdates.date2num(te)
            timestamp = mdates.date2num(t)
            target_time.append(mdates.num2date(timestamp))
            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])
        if file.split('.')[-1] == 'grb2':
            t0 = model[var_name].attrs['initial_time']
            dt = model[var_name].attrs['forecast_time'][0]
            t = datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt))
            timestamp = mdates.date2num(t)
            target_time.append(t)
            lon_model = np.asarray(model.lon_0)
            lat_model = np.asarray(model.lat_0)

        okt = np.logical_and(timestamp_obs_sort > tii,timestamp_obs_sort < tee) 
        sublon_obs = lon_obs_sort[okt]
        sublat_obs = lat_obs_sort[okt]

        if len(sublon_obs) != 0 or len(sublat_obs) != 0:
            okposs = []
            if lon_model.ndim == 1 and lat_model.ndim == 1:
                oklonn = np.round(np.interp(sublon_obs,lon_model,np.arange(len(lon_model)))).astype(int)
                oklatt = np.round(np.interp(sublat_obs,lat_model,np.arange(len(lat_model)))).astype(int)
                for i in np.arange(len(oklonn)):
                    okposs.append([oklatt[i],oklonn[i]])
            if lon_model.ndim == 2 and lat_model.ndim == 2:
                for i in np.arange(len(sublon_obs)):
                    oklatt, oklonn = find_grid_position_hycom(lat_model,lon_model,sublat_obs[i],sublon_obs[i])
                    okposs.append([oklatt,oklonn])

            okpos,okpos_ind = np.unique(okposs,return_index=True,axis=0)

            subvar_obs.append(var_obs_sort[okt][okpos_ind].tolist())

            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name])[0,depth_level,okpos[:,0],okpos[:,1]]
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name])[0,okpos[:,0],okpos[:,1]]
            var[var==0] = np.nan
            subvar_model.append(var.tolist())

    return subvar_obs,subvar_model

#####################################################################
def find_grid_position_hycom(lat,lon,target_lat,target_lon):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklat,oklon

#####################################################################
def taylor_template(angle_lim,std_lim):

    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure(figsize=(9,7))
    tr = PolarAxes.PolarTransform()

    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))

    STDgrid=np.arange(0,std_lim,.5)
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))

    ra0, ra1 =0, angle_lim
    cz0, cz1 = 0, std_lim
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)

    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["left"].label.set_text("Normalized Standard Deviation")
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")

    ax1.axis["bottom"].set_visible(False)
    ax1 = ax1.get_aux_axes(tr)

    plt.grid(linestyle=':',alpha=0.5)

    return fig,ax1

#####################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            get_var_from_model_following_trajectory


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    if exp_names[i][0:5] == 'RTOFS':
        folder_exps.append(scratch_folder[i] + exp_names[i] + '/rtofs.' + cycle[:-2] + '/')
    else:
        folder_exps.append(scratch_folder[i] + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

################################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

#################################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
# Loop the experiments to obtain forecasted track
lon_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lat_forec_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get storm track from trak atcf files
    if hafs[i] == 'hfsa' or hafs[i] == 'hfsb':
        file_track = folder + storm_id + '.' + cycle + '.' + hafs[i] + '.trak.atcfunix'
        print(file_track)

        okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
        lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], _, _, _ = get_storm_track_and_int(file_track,storm_num)

    else:
        lon_forec_track[i,:] = np.nan
        lat_forec_track[i,:] = np.nan

#################################################################################
#%% Read drifter data

url = file_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
platform_code = np.asarray(gdata.platform_code)
wind_speed = np.asarray(gdata.windspd)
sea_level_press = np.asarray(gdata.slp)
sea_surface_temp = np.asarray(gdata.sst)
sea_surface_salt = np.asarray(gdata.salinity)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(times)
times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(times) >= mdates.date2num(tini),\
                         mdates.date2num(times) <= mdates.date2num(tend))

# Fields within time window
timeDr = times[oktimeg]
timestampDr = timestamps[oktimeg]
latDr = latitude[oktimeg]
lonDr = longitude[oktimeg]
platform_codeDr = platform_code[oktimeg]
windspdDr = wind_speed[oktimeg]
slpDr = sea_level_press[oktimeg]
sstDr = sea_surface_temp[oktimeg]
sssDr = sea_surface_salt[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
platform_codeD = platform_codeDr[oklat][oklon]
windspdD = windspdDr[oklat][oklon]
slpD = slpDr[oklat][oklon]
sstD = sstDr[oklat][oklon]
sssD = sssDr[oklat][oklon]

codes = np.unique(platform_codeD)

########################################################################
number_obs = np.empty(len(folder_exps))
number_obs[:] = np.nan
number_obs_used = np.empty(len(folder_exps))
number_obs_used[:] = np.nan
std_subsst_obs = np.empty(len(folder_exps))
std_subsst_obs[:] = np.nan
std_subsst_model = np.empty(len(folder_exps))
std_subsst_model[:] = np.nan
bias = np.empty(len(folder_exps))
bias[:] = np.nan
corr = np.empty(len(folder_exps))
corr[:] = np.nan

# Loop the experiments
for i,folder in enumerate(folder_exps):
    print(folder)

    if ocean[i] == 'hycom':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon = np.asarray(hafs_ocean['Longitude'][:])
        #lat = np.asarray(hafs_ocean['Latitude'][:])
        time_name = 'MT'
        lat_name = 'Latitude'
        lon_name = 'Longitude'
        var_name = 'temperature'
        freq_output_hours = 6

    if ocean[i] == 'mom6':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon = np.asarray(hafs_ocean['xh'][:])
        #lat = np.asarray(hafs_ocean['yh'][:])
        time_name = 'time'
        lat_name = 'yh'
        lon_name = 'xh'
        var_name = 'temp'
        freq_output_hours =3

    if ocean[i] == 'rtofs':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3dz*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon = np.asarray(hafs_ocean['Longitude'][:])
        #lat = np.asarray(hafs_ocean['Latitude'][:])
        time_name = 'MT'
        lon_name = 'Longitude'
        lat_name = 'Latitude'
        var_name = 'temperature'
        freq_output_hours = 6
   
    lonDD,latDD = geo_coord_to_HYCOM_coord(lonD,latD)
    lonDD = np.asarray(lonDD)
    latDD = np.asarray(latDD)

    files_model = files_hafs_ocean
    depth_level = 0
    timestamp_obss = timestampD

    if np.min(lon) < 0:
        lon_obss = lonD
    else: 
        lon_obss = lonDD
    lat_obss = latDD

    var_obss = sstD        

    subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output_hours)

    ########################################################################

    subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
    subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])
    ok1 = np.isfinite(subsst_obs_arr)
    subsst_obs_arra = subsst_obs_arr[ok1]
    subsst_model_arra = subsst_model_arr[ok1]
    ok2 = np.isfinite(subsst_model_arra)
    subsst_obs_array = subsst_obs_arra[ok2]
    subsst_model_array = subsst_model_arra[ok2]
    number_obs[i] = lonD.shape[0]
    number_obs_used[i] = len(subsst_obs_array)
    std_subsst_obs[i] = np.nanstd(subsst_obs_array)
    std_subsst_model[i] = np.nanstd(subsst_model_array)
    bias[i] = np.nanmean(subsst_obs_array) - np.nanmean(subsst_model_array)
    corr[i] = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]
    
    plt.figure()
    plt.plot(subsst_obs_array,subsst_model_array,'.',color=exp_colors[i],markersize=7,markeredgecolor='k')
    plt.plot(np.arange(33),np.arange(33),'-',color='silver',linewidth=2)
    plt.ylim([17,33])
    plt.xlim([17,33])
    plt.xlabel('Drifters Observations',fontsize=14)
    plt.ylabel(exp_labels[i],fontsize=14)
    plt.title('SST',fontsize=18)
    plt.text(26,23,'Total Observations = ' + str(number_obs[i]))
    plt.text(26,22,'Observations used = ' + str(number_obs_used[i]))
    plt.text(26,21,'Bias = ' + str(np.round(bias[i],2)))
    plt.text(26,20,'STD obs = ' + str(np.round(std_subsst_obs[i],2)))
    plt.text(26,19,'STD model = ' + str(np.round(std_subsst_model[i],2)))
    plt.text(26,18,'Corr = ' + str(np.round(corr[i],2)))

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax1 = taylor_template(angle_lim,std_lim)
ax1.plot(0,1,'o',label='Drifters'+' =  '+str(number_obs[0]),color='blue',markersize=10,markeredgecolor='k')

for m in np.arange(len(folder_exps)):  # loop the models
    print('m=',m)
    theta = np.arccos(corr[m])
    rr = std_subsst_model[m]/std_subsst_obs[m]
    ax1.plot(theta,rr,markers[m],color = exp_colors[m],markersize=8,markeredgecolor='k',label=exp_labels[m]+' =  '+str(np.round(number_obs_used[m],2)))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('SST ',fontsize=18)

y0 = 1.1
plt.text(1.4,y0,'Bias',fontsize=14)
for m in np.arange(len(folder_exps)):  # loop the models
    plt.text(1.4,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias[m],3)),fontsize=14)
#plt.savefig('fig7',bbox_inches = 'tight',pad_inches = 0.1)

##################################################################
# Figure track
lev = np.arange(-9000,9100,100)
fig,ax = plt.subplots(figsize=(9,4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonD, latD,'.',color='orangered',label='Drifters')
plt.legend(loc='upper right',bbox_to_anchor=[1.45,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])
