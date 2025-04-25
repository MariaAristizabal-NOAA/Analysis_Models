#%% User input
# forecasting cycle to be used

# Beryl
#cycle = '2024081312'

# Ernesto
#cycle = '2024081612'
#storm_num = '05'
#basin = 'al'
#storm_id = '05l'
#storm_name = 'Ernesto'

# Helena
#cycle = '2024092500'
#storm_num = '09'
#basin = 'al'
#storm_id = '09l'
#storm_name = 'Helena'

# Milton
cycle = '2024100800'
storm_num = '14'
basin = 'al'
storm_id = '14l'
storm_name = 'Milton'

#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1083_hurricane_2024.nc'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1057_hurricane_2024.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1031_hurricane_2024.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1068_hurricane_2024.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1069_hurricane_2024.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1091_hurricane_2024.nc'

exp_names = ['HFSB_oper','hafs_20241220_v2p1b_baseline','hafs_20250306_v2p1b_hb43']
exp_labels = ['HFSB_oper','HAFSv2.1B baseline','HAFSv2.1B final']
exp_colors = ['lime','olive','darkorange']
hafs_ab = ['hfsb','hfsb','hfsb']
ocean = ['hycom','hycom','mom6']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']

'''
exp_names = ['HFSA_oper','hafs_20241220_v2p1a_baseline','hafs_20250210_v2p1a_ha30']
exp_labels = ['HFSA_oper','HAFSv2.1A baseline','HAFSv2.1A final']
exp_colors = ['purple','dodgerblue','#00c8c8']
hafs_ab = ['hfsa','hfsa','hfsa']
ocean = ['mom6','mom6','mom6']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']
'''

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

lon_lim = [-75,-55.0]
lat_lim = [20.0,50.0]

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
#import cfgrib
import pygrib
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
#%% Read Saildrone data
# Sam

url = url_saildrone

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
#water_current_speed_mean = np.asarray(gdata.WATER_CURRENT_SPEED_MEAN)
#water_current_direccion_mean = np.asarray(gdata.WATER_CURRENT_DIRECTION_MEAN)
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
#water_speed = water_current_speed_mean[oktimeg]
#water_dir = water_current_direccion_mean[oktimeg]
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

target_timeS = []
target_tempS = np.empty((len(folder_exps),43))
target_tempS[:] = np.nan
target_saltS = np.empty((len(folder_exps),43))
target_saltS[:] = np.nan
target_timeSfv3 = []
target_surf_pres = np.empty((len(folder_exps),43))
target_surf_pres[:] = np.nan
target_tmp_surf = np.empty((len(folder_exps),43))
target_tmp_surf[:] = np.nan
target_tmp_2mabove = np.empty((len(folder_exps),43))
target_tmp_2mabove[:] = np.nan
target_rh_2mabove = np.empty((len(folder_exps),43))
target_rh_2mabove[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    if hafs_ab[i] == 'hfsa':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if hafs_ab[i] == 'hfsb':
        file_track = folder + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
    print(file_track)

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

    #%% Get list files
    if ocean[i] == 'hycom':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon_hafs_ocean = np.asarray(hafs_ocean['Longitude'][:])
        lat_hafs_ocean = np.asarray(hafs_ocean['Latitude'][:])
        depth_hafs_ocean = np.asarray(hafs_ocean['Z'][:])
        time_name = 'MT'
        temp_name = 'temperature'
        salt_name = 'salinity'
        lon_name = 'Longitude'
        lat_name = 'Latitude'

    if ocean[i] == 'mom6':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
        lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
        depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])
        time_name = 'time'
        temp_name = 'temp'
        salt_name = 'so'
        lon_name = 'xh'
        lat_name = 'yh'

    #%% Read time
    time_ocean = []
    timestamp_ocean = []
    for n,file in enumerate(files_hafs_ocean):
        print(file)
        hafs_ocean = xr.open_dataset(file)
        t = hafs_ocean[time_name][:]
    timestamp = mdates.date2num(t)[0]
    time_ocean.append(mdates.num2date(timestamp))
    timestamp_ocean.append(timestamp)
    time_ocean = np.asarray(time_ocean)
    timestamp_ocean = np.asarray(timestamp_ocean)

    #%% Get list files
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.storm.atm.*.grb2')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))

    #%% Reading FV3 grid
    #fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
    #lon_hafs_fv3 = np.asarray(fv3.lon_0)
    #lat_hafs_fv3 = np.asarray(fv3.lat_0)
    
    #fv3 = cfgrib.open_dataset(files_hafs_fv3[0],filter_by_keys={'typeOfLevel': 'meanSea'})
    #lon_hafs_fv3 = np.asarray(fv3.longitude)
    #lat_hafs_fv3 = np.asarray(fv3.latitude)

    grbindx = pygrib.index(files_hafs_fv3[0],'shortName','typeOfLevel','level')
    selected_grbs = grbindx.select(shortName='t',typeOfLevel='surface',level=0)
    lat_hafs_fv3, lon_hafs_fv3 = selected_grbs[0].latlons()

    #################################################################################
    #%% Retrieve temp. following saildrone trajectory
    ncfiles = files_hafs_ocean
    lon = lon_hafs_ocean
    lat = lat_hafs_ocean
    depth_level = 0
    timestamp_obss = timestampS
    kwargs = dict(depth_level = 0)

    # Conversion from glider longitude and latitude to HYCOM convention
    lonS_hyc, latS_hyc = geo_coord_to_HYCOM_coord(lonS,latS)
    lonS_hyc = np.asarray(lonS_hyc)
    latS_hyc = np.asarray(latS_hyc)

    if np.min(lon_hafs_ocean) < 0:
        lon_obss = lonS
    else:
        lon_obss = lonS_hyc
    lat_obss = latS_hyc

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    #target_timeSocean, target_tempS[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(ncfiles,temp_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)
    target_timeSocean, target_tempS[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,ncfiles,temp_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name,depth_level=0)

    #_, target_saltS[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(ncfiles,salt_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)
    _, target_saltS[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,ncfiles,salt_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name,depth_level=0)

    target_timeS.append(target_timeSocean)

    #################################################################################
    #%% Retrieve HAFS_atm press following saildrone trajectory
    grib2files = files_hafs_fv3
    time_name = 'time'
    lat_name = 'latitude'
    lon_name = 'longitude'
    #var_name = 'prmsl'
    if np.min(lon_hafs_fv3) < 0:
        lon_obss = lonS
    else:
        lon_obss = lonS + 360
    lat_obss = latS
    timestamp_obss = timestampS

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    #target_timeSfv, target_press_surf[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(ncfiles,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv, target_surf_pres[i,0:len(grib2files)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2files,'sp',typeoflevel='surface',level='0')

    target_timeSfv3.append(target_timeSfv)

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonS, latS,'.-',color='blue',label='Saildrone sd'+dataset_id)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend(loc='upper right',bbox_to_anchor=[1.4,1.0])
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
#plt.savefig('track1.png')    

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonS, latS,'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.legend(loc='upper right',bbox_to_anchor=[1.6,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonS)-3,np.nanmax(lonS)+3])
plt.ylim([np.nanmin(latS)-3,np.nanmax(latS)+3])
#plt.savefig('track2.png')    

#################################################################################

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,temp_sb37,'.-',color='blue',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeS[i],target_tempS[i,0:len(target_timeS[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeShyc,target_tempS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
#plt.legend(loc='upper right',bbox_to_anchor=[1.35,1.0])
plt.title('Water Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.savefig('SST.png')    

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,sal_sb37,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeS[i],target_saltS[i,0:len(target_timeS[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeShyc,target_saltS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
#plt.legend(loc='upper right',bbox_to_anchor=[1.35,1.0])
plt.title('Salinity Cycle '+ cycle,fontsize=18)
plt.ylabel(' ',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
plt.ylim([35.0,36.7])
#plt.savefig('SSS.png')    

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,baro_pres,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_surf_pres[i,0:len(target_timeSfv3[i])]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
#plt.legend(loc='lower right')
plt.legend(loc='lower left')
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.ylim([970,1020])
plt.ylim([950,1020])
#plt.savefig('sea_surf_press.png')


