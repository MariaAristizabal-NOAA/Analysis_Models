#%% User input
# forecasting cycle to be used

# Lee
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'
url_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2023.nc'

exp_names = ['HFSAv2a_baseline_MOM6_epbl/baseline','HFSAv2a_baseline_MOM6_epbl/OM4','HFSAv2a_baseline_MOM6_epbl/OM4_LT','HFSAv2a_baseline_MOM6_epbl/RH18','HFSAv2a_baseline_MOM6_epbl/RH18_LT']
exp_labels = ['epbl_baseline','epbl_OM4','epbl_OM4_LT','epbl_RH18','epbl_RH18_LT']
exp_colors = ['olivedrab','mediumspringgreen','mediumpurple','chartreuse','darkorchid']

#exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_HYCOM','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_labels = ['HFSA_oper_HYCOM_kpp','HWRF_POM_MY2.5','HFSAv1p1_HYCOM_kpp','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_colors = ['darkviolet','pink','forestgreen','cyan','royalblue']

lon_lim = [-83,-30]
lat_lim = [10,45]

home_folder = '/home/Maria.Aristizabal/'
home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

# Otis
'''
cycle = '2023102218'
storm_num = '18'
basin = 'ep'
storm_id = '18e'
storm_name = 'otis'
url_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2023.nc'

exp_names = ['HFSA_oper']
exp_labels = ['HFSA_oper_HYCOM_kpp']
exp_colors = ['darkviolet']

lon_lim = [-120,-80]
lat_lim = [0,25]

home_folder = '/home/Maria.Aristizabal/'
home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder1 + exp_names[0] + '/' + cycle + '/' + storm_num + basin[0] + '/']
'''

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
def get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level):

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
        model = xr.open_dataset(file,engine="pynio")
        if file.split('.')[-1] == 'nc':
            t = model[time_name][:]
            timestamp = mdates.date2num(t)[0]
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

        okt = int(np.interp(timestamp,timestamp_obs_sort,np.arange(len(timestamp_obs_sort))))
        oktt = timestamp_obs_sort[0:okt+1] >= timestamp
        #timestamp_obs_sort[0:okt+1][oktt]
        ind = np.arange(okt+1)[oktt]
        sublon_obs = lon_obs_sort[ind]
        sublat_obs = lat_obs_sort[ind]
        subvar_obs.append(var_obs_sort[ind].tolist())

        if lon_model.ndim == 1:
            oklon = np.round(np.interp(sublon_obs,lon_model,np.arange(len(lon_model)))).astype(int)
        if lon_model.ndim == 2:
            oklon = np.round(np.interp(sublon_obs,lon_model[0,:],np.arange(len(lon_model[0,:])))).astype(int)
        if lat_model.ndim == 1:
            oklat = np.round(np.interp(sublat_obs,lat_model,np.arange(len(lat_model)))).astype(int)
        if lat_model.ndim == 2:
            oklat = np.round(np.interp(sublat_obs,lat_model[:,0],np.arange(len(lat_model[:,0])))).astype(int)
        var = np.asarray(model[var_name])[0,depth_level,oklat,oklon]
        var[var==0] = np.nan        
        subvar_model.append(var.tolist())        

    return subvar_obs,subvar_model

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
    folder_exps.append(scratch_folder1 + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

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
'''
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    #%% Get storm track from trak atcf files
    if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSAv1p1_HYCOM' or exp_names[i] == 'HFSAv1p1_MOM6_epbl' or exp_names[i] == 'HFSAv1p1_MOM6_kpp':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if exp_names[i] == 'HFSB':
        file_track = folder + + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
    if exp_names[i] == 'HWRF_2023':
        file_track = folder + storm_name + storm_id + '.' + cycle + '.trak.hwrf.atcfunix'

    print(file_track)
    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
'''

#################################################################################
#%% Read drifter data

url = url_drifter

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

    if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSB' or exp_names[i] == 'HFSAv1p1_HYCOM':
        #%% Get list files
        files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))

        # Reading HAFS/HYCOM grid
        hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
        lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
        lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
        depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

        lonD_hyc, latD_hyc = geo_coord_to_HYCOM_coord(lonD,latD)
        lonD_hyc = np.asarray(lonD_hyc)
        latD_hyc = np.asarray(latD_hyc)

        files_model = files_hafs_hycom
        time_name = 'MT'
        lat_name = 'Latitude'
        lon_name = 'Longitude'
        var_name = 'temperature'
        depth_level = 0
        timestamp_obss = timestampD

        if np.min(lon_hafs_hycom) < 0:
            lon_obss = lonD
        else: 
            lon_obss = lonD_hyc
        lat_obss = latD_hyc

        var_obss = sstD        

        subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level)

    ########################################################################
    if exp_names[i] == 'HWRF_2023':
        #%% Get list files
        files_pom = sorted(glob.glob(os.path.join(folder,'*pom*00*.nc')))
    
        #%% Reading POM grid
        print('Retrieving coordinates from POM')
        grid_file = glob.glob(os.path.join(folder,'*pom.grid.nc'))[0]
        pom_grid = xr.open_dataset(grid_file)
        lon_pom = np.asarray(pom_grid['east_e'][:])
        lat_pom = np.asarray( pom_grid['north_e'][:])
        zlevc = np.asarray(pom_grid['zz'][:])
        topoz = np.asarray(pom_grid['h'][:])
    
        #%% Read POM time
        '''
        time_pom = []
        for n,file in enumerate(files_pom):
            print(file)
            pom = xr.open_dataset(file)
            t = pom.variables['time'][:]
            timestamp = mdates.date2num(t)[0]
            time_pom.append(mdates.num2date(timestamp))
    
        time_pom = np.asarray(time_pom)
        '''    

        ############################################################################
        #%% Retrieve POM temp. following drifters trajectory
    
        files_model = files_pom
        time_name = 'time'
        lat_name = 'north_e'
        lon_name = 'east_e'
        var_name = 't'
        timestamp_obss = timestampD
        lon_obss = lonD
        lat_obss = latD
        var_obss = sstD

        subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level)

    ############################################################################
    if exp_labels[i].split('_')[0] == 'epbl' or exp_labels[i].split('_')[0] == 'kpp' or exp_labels[i] == 'HFSAv1p1' or exp_labels[i] == 'HFSAv2a_baseline' or exp_labels[i] == 'HFSAv2a_baseline_latest':
        
        #%% Get list files
        files_hafs_mom6 = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))

        #%% Reading MOM6 grid
        hafs_mom6_grid = xr.open_dataset(files_hafs_mom6[0],decode_times=False)
        lon_hafs_mom6 = np.asarray(hafs_mom6_grid['xh'][:])
        lat_hafs_mom6 = np.asarray(hafs_mom6_grid['yh'][:])
        depth_hafs_mom6 = np.asarray(hafs_mom6_grid['z_l'][:])

        #%% Read HAFS/HYCOM time
        '''
        time_mom6 = []
        timestamp_mom6 = []
        for n,file in enumerate(files_hafs_mom6):
            print(file)
            MOM6 = xr.open_dataset(file)
            t = MOM6['time'][:]
            timestamp = mdates.date2num(t)[0]
            time_mom6.append(mdates.num2date(timestamp))
            timestamp_mom6.append(timestamp)

        time_mom6 = np.asarray(time_mom6)
        timestamp_mom6 = np.asarray(timestamp_mom6)
        '''

        #################################################################################
        #%% Retrieve HAFS_HYCOM temp. following saildrone trajectory

        files_model = files_hafs_mom6
        time_name = 'time'
        lat_name = 'yh'
        lon_name = 'xh'
        var_name = 'temp'
        depth_level = 0
        timestamp_obss = timestampD
        lon_obss = lonD
        lat_obss = latD
        var_obss = sstD

        subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level)

    #######################################################################        
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
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Drifters',color='orangered',markersize=10,markeredgecolor='k')

for m in np.arange(len(folder_exps)):  # loop the models
    print('m=',m)
    theta = np.arccos(corr[m])
    rr = std_subsst_model[m]/std_subsst_obs[m]
    ax1.plot(theta,rr,'o',color = exp_colors[m],markersize=8,markeredgecolor='k',label=exp_labels[m])

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('SST ',fontsize=18)
#plt.savefig('fig7',bbox_inches = 'tight',pad_inches = 0.1)

