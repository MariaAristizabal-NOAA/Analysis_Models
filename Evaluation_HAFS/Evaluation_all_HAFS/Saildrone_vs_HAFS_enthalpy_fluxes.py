#%% User input
# forecasting cycle to be used

# Lee
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'

#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/sd1064_hurricane_2023.nc'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_update26feb2024.nc'
url_saildrone_adcp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_hurricane_2023_adcp_cbe4_c409_047b_U1696433734196.nc'

exp_names = ['HFSAv2a_baseline_latest']
exp_labels = ['HFSAv2a_baseline']
exp_colors = ['orange']

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

lon_lim = [-75,-55.0]
lat_lim = [20.0,50.0]

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
#folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

# folder CORE-algorithm
folder_COARE = '/home/Maria.Aristizabal/Analysis/COARE-algorithm/Python/COARE3.6'

#############################################################################
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

#sys.path.append(folder_utils4hycom)
#from utils4HYCOM import readgrids,readdepth,readVar,readBinz

#sys.path.append(folder_uom)
#from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            get_var_from_model_following_trajectory

sys.path.append(folder_COARE)
from coare36vn_zrf_et import coare36vn_zrf_et

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#############################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder1 + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

#############################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#############################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

#############################################################################
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

##############################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

##############################################################################
#%% Read Saildrone data

url = url_saildrone
gdata = xr.open_dataset(url)#,decode_times=False)

# Variables in the same order as required by CORE

# u = water-relative wind speed magnitude (m/s) at height zu (m)
wind_from_mean = np.asarray(gdata.WIND_FROM_MEAN)
wind_speed_mean = np.asarray(gdata.WIND_SPEED_MEAN)
water_current_speed_mean = np.asarray(gdata.WATER_CURRENT_SPEED_MEAN)
water_current_direccion_mean = np.asarray(gdata.WATER_CURRENT_DIRECTION_MEAN)

# zu height of wind speed observations
zu = gdata.WIND_SPEED_MEAN.installed_height

# t = air temperature (degC) at height zt (m)
temp_air_mean = np.asarray(gdata.TEMP_AIR_MEAN)

# zt height of air temp. observations
zt = gdata.TEMP_AIR_MEAN.installed_height

# rh = relative humidity (#) at height zq (m)
rh_mean = np.asarray(gdata.RH_MEAN)

# zq height of relative humidity observations
zq = gdata.RH_MEAN.installed_height

# P = sea level air pressure (mb)
baro_pres_mean = np.asarray(gdata.BARO_PRES_MEAN)

# ts = seawater temperature (degC)
temp_sb37_mean = np.asarray(gdata.TEMP_SBE37_MEAN)

# sw_dn = downward (positive) shortwave radiation (W/m^2)
sw_dn = np.tile(200,len(temp_sb37_mean)) # default value inside COARE

# lat = latitude defined positive to north
latitude = np.asarray(gdata.latitude)

# lw_dn = downward (positive) longwave radiation (W/m^2)
lw_dn = 400-1.6*np.abs(latitude) # default value inside COARE

# lon = longitude defined positive to east
longitude = np.asarray(gdata.longitude)

# jd = year day or julian day, where day Jan 1 00:00 UTC = 0
time = np.asarray(gdata.time)
timestamps = mdates.date2num(time)
times = np.asarray(mdates.num2date(timestamps))
JD = np.asarray([t.timetuple().tm_yday for t in times])

# zi = PBL height (m) (default or typical value = 600m)
zi = 600

# rain = rain rate (mm/hr)
rain = 0

#  Ss = sea surface salinity (PSU)
sal_sb37_mean = np.asarray(gdata.SAL_SBE37_MEAN)

# cp = phase speed of dominant waves (m/s) computed from peak period
wave_dominant_period = np.asarray(gdata.WAVE_DOMINANT_PERIOD)

# sigH = significant wave height (m)
wave_significant_height = np.asarray(gdata.WAVE_SIGNIFICANT_HEIGHT)

dataset_id = gdata.drone_id

##############################################################################
# Obtain enthalpy fluxes using COARE algorithm

oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

u = wind_speed_mean[oktimeg]
zu = zu
t = temp_air_mean[oktimeg] 
zt = zt
rh = rh_mean[oktimeg] 
zq = zq
P = baro_pres_mean[oktimeg] 
ts = temp_sb37_mean[oktimeg] 
sw_dn = sw_dn[oktimeg]
lw_dn = lw_dn[oktimeg]
lat = latitude[oktimeg]
lon = longitude[oktimeg]
jd = JD[oktimeg]
zi = zi
rain = rain
Ss = sal_sb37_mean[oktimeg] 
cp = wave_dominant_period[oktimeg] 
sigH = wave_significant_height[oktimeg] 

A = coare36vn_zrf_et(u, zu, t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon, jd, zi, rain, Ss, cp, sigH, zrf_u=10.0, zrf_t=10.0, zrf_q=10.0)

timeS = times[oktimeg]
timestampS = timestamps[oktimeg]

tauS = A[:,1]
shtfluxS = A[:,2]
lhtfluxS = A[:,3]
CdS = A[:,12]
ChS = A[:,13]
CeS = A[:,14]

fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,shtfluxS,'.-',markersize=1,color='lightcoral',label='Sensible heat flux')
ax1.plot(timeS,lhtfluxS,'.-',markersize=1,color='lightseagreen',label='latent heat flux')
ax1.legend(loc='center right')
plt.ylabel('($W/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Enthalpy Fluxes using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,P,'.-',markersize=1,label='Sea Level Pressure')
plt.ylabel('(hPa)',fontsize=14)
ax2.legend()


##############################################################################
# Read adcp saildrone data
'''
url_adcp = url_saildrone_adcp
gdata_adcp = xr.open_dataset(url_adcp)#,decode_times=False)

latitude_adcp = np.asarray(gdata_adcp.latitude)
longitude_adcp = np.asarray(gdata_adcp.longitude)
time_adcp = np.asarray(gdata_adcp.time)
depth_adcp = np.asarray(gdata_adcp.depth)
dataset_id_adcp = gdata_adcp.drone_id
vel_east = np.asarray(gdata_adcp.vel_east)
vel_north = np.asarray(gdata_adcp.vel_north)
'''
##############################################################################
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
target_timeSfv3 = []
target_press_surf = np.empty((len(folder_exps),43))
target_press_surf[:] = np.nan
target_vflux = np.empty((len(folder_exps),43))
target_vflux[:] = np.nan
target_uflux = np.empty((len(folder_exps),43))
target_uflux[:] = np.nan
target_lhtflux = np.empty((len(folder_exps),43))
target_lhtflux[:] = np.nan
target_shtflux = np.empty((len(folder_exps),43))
target_shtflux[:] = np.nan

for i,folder in enumerate(folder_exps):

    file_track = folder + storm_id + '.' + cycle + '.hfsa.trak.atcfunix'

    #%% Get list files
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.storm.atm.*.grb2')))

    #%% Reading FV3 grid
    fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
    lon_hafs_fv3 = np.asarray(fv3.lon_0)
    lat_hafs_fv3 = np.asarray(fv3.lat_0)

    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
        
    #################################################################################
    #%% Retrieve Enthalpu fluxes following saildrone trajectory

    files_model = files_hafs_fv3
    time_name = 'time'
    lat_name = 'latitude'
    lon_name = 'longitude'
    #depth_level = 0
    if np.min(lon_hafs_fv3) < 0:
        lon_obss = lon
    else:
        lon_obss = lon + 360
    lat_obss = lat
    timestamp_obss = timestampS

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    target_timeSfv, target_press_surf[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)
     
    _, target_uflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'UFLX_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_vflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'VFLX_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_lhtflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'LHTFL_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_shtflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'SHTFL_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv3.append(target_timeSfv)

#########################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon, lat,'.-',color='blue',label='Saildrone sd'+dataset_id)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend(loc='upper right',bbox_to_anchor=[1.4,1.0])
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
#plt.savefig('track1.png')    

#########################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon, lat,'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.legend(loc='upper right',bbox_to_anchor=[1.6,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lon)-3,np.nanmax(lon)+3])
plt.ylim([np.nanmin(lat)-3,np.nanmax(lat)+3])
#plt.savefig('track2.png')    

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,P,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_press_surf[i,0:len(target_timeSfv3[i])]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.ylim([950,1050])
plt.ylim([1000,1020])
#plt.savefig('sea_surf_press.png')

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,shtfluxS,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_shtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Sensible Heat Flux '+ cycle,fontsize=18)
plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.savefig('sea_surf_press.png')

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,lhtfluxS,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Latent Heat Flux '+ cycle,fontsize=18)
plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.savefig('sea_surf_press.png')

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,CdS,'.-',color='blue',label='Saildrone sd'+dataset_id)
#for i in np.arange(len(exp_names)):
#    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Cd '+ cycle,fontsize=18)
#plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,ChS,'.-',color='blue',label='Saildrone sd'+dataset_id)
#for i in np.arange(len(exp_names)):
#    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Ch '+ cycle,fontsize=18)
#plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

