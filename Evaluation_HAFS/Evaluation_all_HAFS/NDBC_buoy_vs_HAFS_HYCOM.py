#%% User input
# forecasting cycle to be used
cycle = '2021082800'
storm_num = '09'
basin = 'al'

exp_names = ['hafsv0p2a_2021rt_natl','HAFSv0p3_baseline','HAFSv0p3_HF3A']
exp_labels = ['HAFSv0.2a','H3BL','HF3A']
exp_colors = ['limegreen','indianred','c']

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle,
               scratch_folder + exp_names[1] + '/' + cycle]

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = scratch_folder + 'adeck/a' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

url_sd1045 = scratch_folder + 'Data/Saildrones/sd1045_hurricane_2021_19cb_1f86_cf24.nc'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob
import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

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
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

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
#%% Read GFS track
lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)
 
#################################################################################
#%% Read Saildrone data
# Sam
if cycle == '2021093006':
    url = url_sd1045

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
time = np.asarray(gdata.time)
dataset_id = gdata.drone_id
temp_air_mean = np.asarray(gdata.TEMP_AIR_MEAN)
rh_mean = np.asarray(gdata.RH_MEAN)
baro_pres_mean = np.asarray(gdata.BARO_PRES_MEAN)
temp_sb37_mean = np.asarray(gdata.TEMP_SBE37_MEAN)
wind_from_mean = np.asarray(gdata.WIND_FROM_MEAN)
wind_speed_mean = np.asarray(gdata.WIND_SPEED_MEAN)
sal_sb37_mean = np.asarray(gdata.SAL_SBE37_MEAN)
water_current_speed_mean = np.asarray(gdata.WATER_CURRENT_SPEED_MEAN)
water_current_direccion_mean = np.asarray(gdata.WATER_CURRENT_DIRECTION_MEAN)
wave_dominant_period = np.asarray(gdata.WAVE_DOMINANT_PERIOD)
wave_significant_height = np.asarray(gdata.WAVE_SIGNIFICANT_HEIGHT)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(time)
times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
latS = latitude[oktimeg]
lonS = longitude[oktimeg]
temp_air = temp_air_mean[oktimeg]
rh = rh_mean[oktimeg]
baro_pres =  baro_pres_mean[oktimeg]
temp_sb37 = temp_sb37_mean[oktimeg]
wind_from = wind_from_mean[oktimeg]
wind_speed = wind_speed_mean[oktimeg]
sal_sb37 = sal_sb37_mean[oktimeg]
water_speed = water_current_speed_mean[oktimeg]
water_dir = water_current_direccion_mean[oktimeg]
wave_dom_period = wave_dominant_period[oktimeg]
wave_dom_height = wave_significant_height[oktimeg]
timeS = times[oktimeg]
timestampS = timestamps[oktimeg]

#################################################################################
#%% Loop the experiments

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan
rmw_track = np.empty((len(folder_exps),43))
rmw_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.nc')))

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], rmw_track[i,0:okn] = get_storm_track_and_int(file_track,storm_num)

    #%% Read HAFS/Fv3 time
    time_fv3 = []
    for n,file in enumerate(files_hafs_fv3):
        print(file)
        FV3 = xr.open_dataset(file)
        t = FV3.variables['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_fv3.append(mdates.num2date(timestamp))

    time_fv3 = np.asarray(time_fv3)

    #%% Read HAFS/HYCOM time
    time_hycom = []
    timestamp_hycom = []
    for n,file in enumerate(files_hafs_hycom):
        print(file)
        HYCOM = xr.open_dataset(file)
        t = HYCOM['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time_hycom.append(mdates.num2date(timestamp))
        timestamp_hycom.append(timestamp)

    time_hycom = np.asarray(time_hycom)
    timestamp_hycom = np.asarray(timestamp_hycom)

    #################################################################################
    #%% Retrieve HAFS_HYCOM temp. following saildrone trajectory

    # Conversion from glider longitude and latitude to HYCOM convention
    lonS_hyc, latS_hyc = geo_coord_to_HYCOM_coord(lonS,latS)

    files_model = files_hafs_hycom
    time_name = 'MT'
    lat_name = 'Latitude'
    lon_name = 'Longitude'
    depth_level = 0
    lon_obs = lonS_hyc
    lat_obs = latS_hyc
    timestamp_obs = timestampS
    kwargs = dict(depth_level = 0)

    target_timeShyc, target_tempS = get_var_from_model_following_trajectory(files_model,'temperature',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

    target_timeShyc, target_saltS = get_var_from_model_following_trajectory(files_model,'salinity',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

    #################################################################################
    #%% Retrieve HAFS_atm press following saildrone trajectory

    files_model = files_hafs_fv3
    time_name = 'time'
    lat_name = 'latitude'
    lon_name = 'longitude'
    depth_level = 0
    lon_obs = lonS
    lat_obs = latS
    timestamp_obs = timestampS

    target_timeSfv3, target_press = get_var_from_model_following_trajectory(files_model,'PRES_surface',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv3, target_tmp_surf = get_var_from_model_following_trajectory(files_model,'TMP_surface',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv3, target_tmp_2mabove = get_var_from_model_following_trajectory(files_model,'TMP_2maboveground',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv3, target_rh_2mabove = get_var_from_model_following_trajectory(files_model,'RH_2maboveground',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_forec_track_hafs[::2], lat_forec_track_hafs[::2],'o-',color='c',markeredgecolor='k',label='HAFS-HYCOM',markersize=7)
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
plt.plot(lonS, latS,'.-',color='blue',label='Saildrone Track')
plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
# Elsa
if cycle == '2021070100':
    plt.xlim([-85,-40])
    plt.ylim([5,30])
# Ida
if cycle == '2021092800':
    plt.xlim([-95,-75])
    plt.ylim([20,40])
# Sam
if cycle == '2021093006':
    plt.xlim([-63,-40])
    plt.ylim([18,50])

#################################################################################

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,temp_sb37,'.-',color='blue',label='Saildrone')
plt.plot(target_timeShyc,target_tempS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('Water Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,sal_sb37,'.-',color='blue',label='Saildrone')
plt.plot(target_timeShyc,target_saltS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('Salinity Cycle '+ cycle,fontsize=18)
plt.ylabel(' ',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,baro_pres,'.-',color='blue',label='Saildrone')
plt.plot(target_timeSfv3,target_press/100,'o-',color='c',label='HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,temp_air,'.-',color='blue',label='Saildrone')
plt.plot(target_timeSfv3,target_tmp_2mabove-273.15,'o-',color='c',label='tmp_2maboveground HAFS-HYCOM',markeredgecolor='k',markersize=10)
plt.legend()
plt.title('Air Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,rh,'.-',color='b',label='Saildrone')
plt.plot(target_timeSfv3,target_rh_2mabove,'o-',color='c',label='2m above HAFS-HYCOM',markeredgecolor='k',markersize=10)
plt.legend()
plt.title('Relative humidity Cycle '+ cycle,fontsize=18)
plt.ylabel('%',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

