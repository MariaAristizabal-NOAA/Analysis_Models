#%% User input
# forecasting cycle to be used

# Lee
'''
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_update26feb2024.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/sd1064_hurricane_2023.nc'
#url_saildrone_adcp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_hurricane_2023_adcp_cbe4_c409_047b_U1696433734196.nc'

exp_names = ['HAFSv2_h3a0']
exp_labels = ['HFSAv2a_baseline']
exp_colors = ['orange']
'''

'''
cycle = '2024081612'
storm_num = '05'
basin = 'al'
storm_id = '05l'
storm_name = 'Ernesto'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1031_hurricane_2024.nc'
'''

cycle = '2024092500'
storm_num = '09'
basin = 'al'
storm_id = '09l'
storm_name = 'Helene'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1057_hurricane_2024.nc'

exp_names = ['HFSA_oper','HFSB_oper','HAFSv2p0p1a_2024rt']
exp_labels = ['HFSA_oper','HFSB_oper','HFSA_para']
exp_colors = ['purple','lime','dodgerblue']
hafs_ab = ['hfsa','hfsb','hfsa']

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
Sw_dn = np.tile(200,len(temp_sb37_mean)) # default value inside COARE

# lat = latitude defined positive to north
latitude = np.asarray(gdata.latitude)

# lw_dn = downward (positive) longwave radiation (W/m^2)
Lw_dn = 400-1.6*np.abs(latitude) # default value inside COARE

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
ok = np.isfinite(u)

u = u[ok]
zu = zu
t = temp_air_mean[oktimeg][ok]
zt = zt
rh = rh_mean[oktimeg][ok]
zq = zq
P = baro_pres_mean[oktimeg][ok]
ts = temp_sb37_mean[oktimeg][ok]
sw_dn = Sw_dn[oktimeg][ok]
lw_dn = Lw_dn[oktimeg][ok]
lat = latitude[oktimeg][ok]
lon = longitude[oktimeg][ok]
jd = JD[oktimeg][ok]
zi = zi
rain = rain
Ss = sal_sb37_mean[oktimeg][ok]
cp = wave_dominant_period[oktimeg][ok]
sigH = wave_significant_height[oktimeg][ok]

A = coare36vn_zrf_et(u, zu, t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon, jd, zi, rain, Ss, cp, sigH, zrf_u=10.0, zrf_t=10.0, zrf_q=10.0)

timeS = times[oktimeg][ok]
timestampS = timestamps[oktimeg][ok]

usr= A[:,0]
tauS = A[:,1]
shtfluxS = A[:,2]
lhtfluxS = A[:,3]
CdS = A[:,12]
ChS = A[:,13]
CeS = A[:,14]
UrfN = A[:,20]
Cdn_10 = A[:,34]
Chn_10 = A[:,35]
Cen_10 = A[:,36]

##############################################################
UrfN[np.abs(UrfN)>100] = np.nan
tauS[tauS > 100] = np.nan
fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,tauS,'.-',markersize=1,color='orange',label='Wind Stress')
ax1.legend(loc='center right')
plt.ylabel('($N/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Wind Stress using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,P,'.-',markersize=1,label='Sea Level Pressure')
plt.ylabel('(hPa)',fontsize=14)
ax2.legend()

##############################################################
shtfluxS[np.abs(shtfluxS) > 2000] = np.nan
lhtfluxS[np.abs(lhtfluxS) > 2000] = np.nan
fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,shtfluxS,'.-',markersize=1,color='lightcoral',label='Sensible heat flux')
ax1.plot(timeS,lhtfluxS,'.-',markersize=1,color='lightseagreen',label='latent heat flux')
ax1.legend(loc='lower left')
plt.ylabel('($W/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Enthalpy Fluxes using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,P,'.-',markersize=1,label='Sea Level Pressure')
plt.ylabel('(hPa)',fontsize=14)
ax2.legend(loc='upper left')
#ax2.legend()

##############################################################
shtfluxS[np.abs(shtfluxS) > 2000] = np.nan
lhtfluxS[np.abs(lhtfluxS) > 2000] = np.nan
fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,shtfluxS,'.-',markersize=1,color='lightcoral',label='Sensible heat flux')
ax1.plot(timeS,lhtfluxS,'.-',markersize=1,color='lightseagreen',label='latent heat flux')
ax1.plot(timeS,tauS,'.-',markersize=1,color='orange',label='Wind Stress')
#ax1.legend(loc='center right')
#ax1.legend(loc='center right')
ax1.legend(loc='lower left')
plt.ylabel('($W/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Enthalpy Fluxes using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,u,'.-',markersize=1,label='Wind Speed at Installed Height')
plt.ylabel('(m/s)',fontsize=14)
ax2.legend()

##############################################################
fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,tauS,'.-',markersize=1,color='orange',label='Wind Stress')
#ax1.legend(loc='center right')
ax1.legend(loc='lower left')
plt.ylabel('($W/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Enthalpy Fluxes using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,u,'.-',markersize=1,label='Wind Speed at Installed Height')
plt.ylabel('(m/s)',fontsize=14)
ax2.legend()

##############################################################################
# Read adcp saildrone data
'''
#url_adcp = url_saildrone_adcp
#gdata_adcp = xr.open_dataset(url_adcp)#,decode_times=False)

#latitude_adcp = np.asarray(gdata_adcp.latitude)
#longitude_adcp = np.asarray(gdata_adcp.longitude)
#time_adcp = np.asarray(gdata_adcp.time)
#depth_adcp = np.asarray(gdata_adcp.depth)
#dataset_id_adcp = gdata_adcp.drone_id
#vel_east = np.asarray(gdata_adcp.vel_east)
#vel_north = np.asarray(gdata_adcp.vel_north)
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
target_vflux = np.empty((len(folder_exps),43))
target_vflux[:] = np.nan
target_uflux = np.empty((len(folder_exps),43))
target_uflux[:] = np.nan
target_lhtflux = np.empty((len(folder_exps),43))
target_lhtflux[:] = np.nan
target_shtflux = np.empty((len(folder_exps),43))
target_shtflux[:] = np.nan
target_surf_pres = np.empty((len(folder_exps),43))
target_surf_pres[:] = np.nan
target_surf_temp = np.empty((len(folder_exps),43))
target_surf_temp[:] = np.nan
target_fric_vel = np.empty((len(folder_exps),43))
target_fric_vel[:] = np.nan
target_10u = np.empty((len(folder_exps),43))
target_10u[:] = np.nan
target_10v = np.empty((len(folder_exps),43))
target_10v[:] = np.nan
target_2temp = np.empty((len(folder_exps),43))
target_2temp[:] = np.nan
target_10temp = np.empty((len(folder_exps),43))
target_10temp[:] = np.nan

for i,folder in enumerate(folder_exps):

    #%% Get storm track from trak atcf files
    if hafs_ab[i] == 'hfsa':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if hafs_ab[i] == 'hfsb':
        file_track = folder + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
    print(file_track)

    #%% Get list files
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.storm.atm.*.grb2')))
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))

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
    #fv3s = pygrib.open(files_hafs_fv3[0])
    #fv3 = fv3s.select(name='Temperature')[0]
    #lat_hafs_fv3, lon_hafs_fv3 = fv3.latlons()

    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
        
    #################################################################################
    #%% Retrieve Enthalpy fluxes following saildrone trajectory

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

    #fv3 = cfgrib.open_dataset(files_hafs_fv3[0],filter_by_keys={'typeOfLevel': 'meanSea'})
    #kwargs = dict(filter_by_keys={'typeOfLevel': 'meanSea'})

    #fv3 = cfgrib.open_dataset(files_hafs_fv3[0],filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})
   # kwargs = dict(filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})

    #Latent heat net flux
    #fv3s.select(name='Latent heat net flux')[0].values
    #grbindx.select(shortName='lhtfl',typeOfLevel='surface',level=0)[0].values
    target_timeSfv, target_lhtflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'lhtfl','surface','0')

    #Instantaneous surface sensible heat flux
    #grbindx.select(shortName='ishf',typeOfLevel='surface',level=0)[0].values
    _, target_shtflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'ishf','surface','0')

    #U-component of atmospheric surface momentum flux
    _, target_uflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'utaua','surface','0')

    #V-component of atmospheric surface momentum flux
    _, target_vflux[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'vtaua','surface','0')

    #Surface pressure
    _, target_surf_pres[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'sp','surface','0')

    #Temperature:K (instant):regular_ll:surface:level 0
    _, target_surf_temp[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'t','surface','0')
    
    #Frictional velocity
    _, target_fric_vel[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'fricv','surface','0')

    #10 metre U wind component
    #grbindx.select(shortName='10u',typeOfLevel='heightAboveGround',level=10)[0].values
    _, target_10u[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'10u','heightAboveGround','10')

    #10 metre V wind component
    #grbindx.select(shortName='10v',typeOfLevel='heightAboveGround',level=10)[0].values
    _, target_10v[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'10v','heightAboveGround','10')

    #2 metre temperature
    #grbindx.select(shortName='2t',typeOfLevel='heightAboveGround',level=2)[0].values
    _, target_2temp[i,0:len(files_model)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,'2t','heightAboveGround','2')

    target_timeSfv3.append(target_timeSfv)

#########################################################################
# Calculate CD
#rho_air = target_surf_pres/(target_surf_temp * 287.058)
rho_air = target_surf_pres/(target_2temp * 287.058)
target_10w = np.sqrt(target_10u**2 + target_10v**2)
target_ustar = np.sqrt(target_uflux**2 + target_vflux**2)/rho_air
CD = target_ustar**2/target_10w**2

#CH = target_shtflux/(rho_air * 1025 * target_10w * (target_surf_temp-target_2temp)) 
CH = target_shtflux/(rho_air * 1004 * target_10w * (target_surf_temp-target_2temp)) 

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
    plt.plot(target_timeSfv3[i],target_surf_pres[i,0:len(target_timeSfv3[i])]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
plt.ylim([970,1025])
#plt.ylim([1000,1020])
#plt.savefig('sea_surf_press.png')

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],shtfluxS[0:779],'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(timeS,shtfluxS,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_shtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend(loc='upper left')
plt.title('Sensible Heat Flux '+ cycle,fontsize=18)
plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.savefig('sea_surf_press.png')

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],lhtfluxS[0:779],'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(timeS,lhtfluxS,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend(loc='upper left')
plt.title('Latent Heat Flux '+ cycle,fontsize=18)
plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.savefig('sea_surf_press.png')

#########################################################################
target_momflux = np.sqrt(target_uflux**2 + target_vflux**2)
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],lhtfluxS[0:779],'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(timeS,tauS,'.-',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_momflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend(loc='upper left')
plt.title('Wind Stress '+ cycle,fontsize=18)
plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

#########################################################################
UrfN[np.abs(UrfN)>100] = np.nan
Cdn_10[Cdn_10>0.1] = np.nan
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],Cdn_10[0:779]*1000,'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(timeS,Cdn_10*1000,'.-',color='blue',label='Saildrone sd'+dataset_id)
#plt.legend(loc='lower right')
plt.legend(loc='lower left')
plt.title('Cdn_10 x 1000 '+ cycle,fontsize=18)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

#########################################################################
Chn_10[Chn_10*1000>20] = np.nan
Chn_10[Chn_10<0] = np.nan
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],Chn_10[0:779]*1000,'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(timeS,Chn_10*1000,'.-',color='blue',label='Saildrone sd'+dataset_id)
plt.legend(loc='upper left')
plt.title('Chn_10 x 1000 '+ cycle,fontsize=18)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

#########################################################################
UrfN[np.abs(UrfN)>100] = np.nan
Cdn_10[Cdn_10>0.1] = np.nan
fig,ax = plt.subplots(figsize=(10, 5))
#plt.plot(UrfN[0:779],Cdn_10[0:779]*1000,'.',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(UrfN,Cdn_10*1000,'.',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_10w[i],CD[i]*1000,'o',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend(loc='upper left')
plt.title('Drag Coefficient '+ cycle,fontsize=18)
plt.xlabel('10-m Wind Speed (m/s)',fontsize=16) 
plt.ylabel('Cdn_10 x 1000',fontsize=16) 

#########################################################################
Chn_10[Chn_10*1000>20] = np.nan
Chn_10[Chn_10<0] = np.nan
fig,ax = plt.subplots(figsize=(10, 5))
#plt.plot(UrfN[0:780],Chn_10[0:780]*1000,'.',color='blue',label='Saildrone sd'+dataset_id)
plt.plot(UrfN,Chn_10*1000,'.',color='blue',label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_10w[i,:],CH[i,:]*1000,'o',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend(loc='upper left')
plt.title('Sensible/Latent Heat Transfer Coefficient '+ cycle,fontsize=18)
plt.xlabel('10-m Wind Speed (m/s)',fontsize=16) 
plt.ylabel('Chn_10 x 1000',fontsize=16) 

#########################################################################
'''
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(UrfN,Cen_10*1000,'.',color='blue',label='Saildrone sd'+dataset_id)
#for i in np.arange(len(exp_names)):
#    plt.plot(target_10w[i],CE[i]*1000,'.',color=exp_colors[i],label=exp_labels[i])
plt.legend(loc='lower right')
plt.title('Latent Heat Transfer Coefficient '+ cycle,fontsize=18)
plt.xlabel('10-m Wind Speed (m/s)',fontsize=16) 
plt.ylabel('Cen_10 x 1000',fontsize=16) 
'''
#########################################################################
UrfN[np.abs(UrfN)>100] = np.nan
fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timeS[0:779],UrfN[0:779],'.-',markersize=1,label='Saildrone sd'+dataset_id)
plt.plot(timeS,UrfN,'.-',markersize=1,label='Saildrone sd'+dataset_id)
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3[i],target_10w[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='center right')
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
plt.ylabel('(m/s)',fontsize=14)
plt.title('10 m Wind Speed'+ cycle,fontsize=18)

#########################################################################
'''
# Calculate CD
#rho_air = target_surf_pres/(target_surf_temp * 287.058)
rho_air = target_surf_pres/(target_2temp * 287.058)
target_10w = np.sqrt(target_10u**2 + target_10v**2)
target_ustar = np.sqrt(target_uflux**2 + target_vflux**2)/rho_air
CD = target_ustar**2/target_10w**2

#CH = target_shtflux/(rho_air * 1025 * target_10w * (target_surf_temp-target_2temp))
CH = target_shtflux/(rho_air * 1004 * target_10w * (target_surf_temp-target_2temp))

i=0
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],target_10u[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='10u',markersize=7)
plt.plot(target_timeSfv3[i],target_10v[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='10v',markersize=7)
plt.plot(target_timeSfv3[i],target_10w[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='10w',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],target_uflux[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='uflx',markersize=7)
plt.plot(target_timeSfv3[i],target_vflux[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='vflx',markersize=7)
plt.plot(target_timeSfv3[i],target_ustar[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='ustar',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],CD[i,0:len(target_timeSfv3[i])]*1000,'o-',markeredgecolor='k',label='CD*1000',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],target_shtflux[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='shtflx',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],target_surf_temp[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='surf_temp',markersize=7)
plt.plot(target_timeSfv3[i],target_2temp[i,0:len(target_timeSfv3[i])],'o-',markeredgecolor='k',label='2temp',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(target_timeSfv3[i],CH[i,0:len(target_timeSfv3[i])]*1000,'o-',markeredgecolor='k',label='CH*1000',markersize=7)
plt.legend()
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
'''
