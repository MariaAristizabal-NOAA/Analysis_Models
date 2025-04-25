#%% User input
# forecasting cycle to be used

# User input
# Helene
cycle = '2024092312'
storm_num = '09'
basin = 'al'
storm_id = '09l'
storm_name = 'Helene'
fhour = 0

# Milton
#cycle = '2024100800'
#storm_num = '14'
#basin = 'al'
#storm_id = '14l'
#storm_name = 'Milton'

exp_names = ['RTOFS','RTOFS_v2.5.test01','HFSA_oper','HFSB_oper','hafsv2p0p1a_2024rt_cplinit_AOBS']
exp_labels = ['RTOFS_oper','RTOFSv2.5','HFSA_oper','HFSB_oper','HAFS_cplinit_AOBS']
exp_colors = ['magenta','salmon','purple','lime','cyan']
hafs = [' ',' ','hfsa','hfsb','hfsa']
ocean = ['rtofs','rtofs','mom6','hycom','mom6']
markers = ['p','p','o','o','o']

folder_CoralTemp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Coral_Reef_Watch/2024/'

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

############################################################
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
from scipy.interpolate import LinearNDInterpolator
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

#####################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    if exp_names[i][0:5] == 'RTOFS':
        folder_exps.append(scratch_folder + exp_names[i] + '/rtofs.' + cycle[:-2] + '/')
    else:
        folder_exps.append(scratch_folder + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

#####################################################################
# target time
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
target_time = tini + timedelta(hours=fhour)
target_timestamp = mdates.date2num(target_time)

#####################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

####################################################################
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

#####################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#####################################################################
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

#####################################################################
#%% Read Coral Reef Watch SST
files_satell = np.asarray(sorted(glob.glob(os.path.join(folder_CoralTemp,'coraltemp_v3.1_*.nc'))))

var_name = 'analysed_sst'
time_name = 'time'
lon_name = 'lon'
lat_name = 'lat'

t_satell = []
timestamp_satell = []
for f,file in enumerate(files_satell):
    print(f,' ',file)
    satell = xr.open_dataset(file)
    t_satell.append(mdates.num2date(mdates.date2num(np.asarray(satell[time_name][:])))[0])
    timestamp_satell.append(mdates.date2num(np.asarray(satell[time_name][:]))[0])

timestamp_satell = np.asarray(timestamp_satell)
okt = np.argsort(timestamp_satell)
files_satell_sorted = files_satell[okt]
timestamp_satell_sorted =  timestamp_satell[okt]

okf = [t <= target_timestamp for t in timestamp_satell_sorted]
file_satell = files_satell_sorted[okf][-1]

satell = xr.open_dataset(file_satell,decode_times=False)
lon_satell = np.asarray(satell[lon_name][:])
lat_satell = np.asarray(satell[lat_name][:])

oklatsat = np.logical_and(lat_satell >= lat_lim[0],lat_satell <= lat_lim[-1])
oklonsat = np.logical_and(lon_satell >= lon_lim[0],lon_satell <= lon_lim[-1])

lat_sat = lat_satell[oklatsat]
lon_sat = lon_satell[oklonsat]
sst_sat = np.asarray(satell[var_name][0,:,:])[oklatsat,:][:,oklonsat]

lon_satH,lat_satH = geo_coord_to_HYCOM_coord(lon_sat,lat_sat)
lon_satH = np.asarray(lon_satH)
lat_satH = np.asarray(lat_satH)

#fig = plt.figure(figsize=(7,4))
plt.figure()
levels = np.arange(18,34)
plt.contourf(lon_sat,lat_sat,sst_sat,levels=levels,cmap='turbo',extend='both')
plt.colorbar(extendrect=True)
plt.title('CoralTemp SST Cycle='+cycle+' fhour='+str(fhour),fontsize=18)
plt.axis('scaled')

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
        time_name = 'MT'
        lat_name = 'Latitude'
        lon_name = 'Longitude'
        var_name = 'temperature'

    if ocean[i] == 'mom6':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        time_name = 'time'
        lat_name = 'yh'
        lon_name = 'xh'
        var_name = 'temp'

    if ocean[i] == 'rtofs':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3dz*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        time_name = 'MT'
        lon_name = 'Longitude'
        lat_name = 'Latitude'
        var_name = 'temperature'
   
    files_model = np.asarray(files_hafs_ocean)
    
    #read time
    t_model = []
    timestamp_model = []
    for f,file in enumerate(files_model):
        print(f,' ',file)
        model = xr.open_dataset(file)
        t_model.append(mdates.num2date(mdates.date2num(np.asarray(model[time_name][:])))[0])
        timestamp_model.append(mdates.date2num(np.asarray(model[time_name][:]))[0])

    timestamp_model = np.asarray(timestamp_model)
    okt = np.argsort(timestamp_model)
    files_model_sorted = files_model[okt]
    timestamp_model_sorted =  timestamp_model[okt]

    #okf = [t <= target_timestamp for t in timestamp_model_sorted]
    okf = int(np.round(np.interp(target_timestamp,timestamp_model_sorted,np.arange(len(timestamp_model_sorted)))))

    file_model = files_model_sorted[okf]

    model = xr.open_dataset(file_model,decode_times=False)
    lon_model = np.asarray(model[lon_name][:])
    lat_model = np.asarray(model[lat_name][:])
    sst_model = np.asarray(model[var_name][0,0,:,:])
    #oklatmod = np.logical_and(lat_mod >= lat_lim[0],lat_mod <= lat_lim[-1])
    #oklonmod = np.logical_and(lon_mod >= lon_lim[0],lon_mod <= lon_lim[-1])
    #lat_model = lat_mod[oklatmod,:][:,oklonmod]
    #lon_model = lon_mod[oklonmod]
    #sst_model = np.asarray(model[var_name][0,0,:,:])[oklatmod,:][:,oklonmod]
    
    if np.min(lon_model) < 0:
        lon_obs = lon_sat
    else: 
        lon_obs = lon_satH
    lat_obs = lat_sat

    #Spatial interpolation of model to CoralTemp spatial resolution
    lono,lato = np.meshgrid(lon_obs,lat_obs)
    if lon_model.ndim == 2 and lat_model.ndim ==2:
        lonm = lon_model
        latm = lat_model
    else:
        lonm,latm = np.meshgrid(lon_model,lat_model)

    interpolator = LinearNDInterpolator(list(zip(np.ravel(lonm),np.ravel(latm))),np.ravel(sst_model))
    sst_from_model_to_sat = interpolator((lono,lato))

    fig = plt.figure(figsize=(6,4))
    levels = np.arange(18,34)
    plt.contourf(lon_sat,lat_sat,sst_from_model_to_sat,levels=levels,cmap='turbo',extend='both')
    plt.colorbar(extendrect=True)
    plt.title(exp_labels[i]+' SST Cycle='+cycle+' fhour='+str(fhour),fontsize=16)
    plt.axis('scaled')

    ##################################################################
    sst_sat_vector = np.ravel(sst_sat)
    sst_from_model_to_sat_vector = np.ravel(sst_from_model_to_sat)
    ok1 = np.isfinite(sst_sat_vector)
    sst_sat_vect = sst_sat_vector[ok1]
    sst_from_model_to_sat_vect = sst_from_model_to_sat_vector[ok1]
    ok2 = np.isfinite(sst_from_model_to_sat_vect)
    sst_sat_vec = sst_sat_vect[ok2]
    sst_from_model_to_sat_vec = sst_from_model_to_sat_vect[ok2]

    number_obs[i] = len(sst_sat_vect)
    number_obs_used[i] = len(sst_from_model_to_sat_vec)
    std_subsst_obs[i] = np.nanstd(sst_sat_vec)
    std_subsst_model[i] = np.nanstd(sst_from_model_to_sat_vec)
    bias[i] = np.nanmean(sst_sat_vec) - np.nanmean(sst_from_model_to_sat_vec)
    corr[i] = np.corrcoef(sst_sat_vec,sst_from_model_to_sat_vec)[0,1]
    
    plt.figure()
    plt.plot(sst_sat_vec,sst_from_model_to_sat_vec,'.',color=exp_colors[i],markersize=7,markeredgecolor='k')
    plt.plot(np.arange(35),np.arange(35),'-',color='silver',linewidth=2)
    plt.ylim([17,33])
    plt.xlim([17,33])
    plt.xlabel('CoralTemp',fontsize=14)
    plt.ylabel(exp_labels[i],fontsize=14)
    plt.title('SST Cycle='+cycle+' fhour='+str(fhour),fontsize=18)
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
ax1.plot(0,1,'o',label='CoralTemp'+' =  '+str(number_obs[0]),color='blue',markersize=10,markeredgecolor='k')

for m in np.arange(len(folder_exps)):  # loop the models
    print('m=',m)
    theta = np.arccos(corr[m])
    rr = std_subsst_model[m]/std_subsst_obs[m]
    ax1.plot(theta,rr,markers[m],color = exp_colors[m],markersize=8,markeredgecolor='k',label=exp_labels[m]+' =  '+str(np.round(number_obs_used[m],2)))

plt.legend(loc='upper left',bbox_to_anchor=[0.9,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('SST Cycle='+cycle+' fhour='+str(fhour),fontsize=18)

y0 = 1.1
plt.text(1.5,y0,'Bias',fontsize=14)
for m in np.arange(len(folder_exps)):  # loop the models
    plt.text(1.5,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias[m],3)),fontsize=14)
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
plt.legend(loc='upper right',bbox_to_anchor=[1.45,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim(lon_lim)
plt.ylim(lat_lim)
