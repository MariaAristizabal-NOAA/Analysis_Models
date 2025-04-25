#%% User input
cycle = '2023090706'

url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_update26feb2024.nc'
url_saildrone_adcp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_hurricane_2023_adcp_cbe4_c409_047b_U1696433734196.nc'

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

sys.path.append(folder_COARE)
from coare36vn_zrf_et import coare36vn_zrf_et

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#############################################################################
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

#############################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

##############################################################################
# Obtain enthalpy fluxes using COARE algorithm

oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

'''
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
'''

'''
lenn = len(wind_speed_mean[oktimeg])

u = np.tile(np.nanmean(wind_speed_mean[oktimeg]),lenn)
zu = zu
t = np.tile(np.nanmean(temp_air_mean[oktimeg]),lenn) 
zt = zt
rh = np.tile(np.nanmean(rh_mean[oktimeg]),lenn) 
zq = zq
P = np.tile(np.nanmean(baro_pres_mean[oktimeg]),lenn) 
ts = np.tile(np.nanmean(temp_sb37_mean[oktimeg]),lenn) 
sw_dn = np.tile(np.nanmean(sw_dn[oktimeg]),lenn)
lw_dn = np.tile(np.nanmean(lw_dn[oktimeg]),lenn)
lat = np.tile(np.nanmean(latitude[oktimeg]),lenn)
lon = np.tile(np.nanmean(longitude[oktimeg]),lenn)
jd = np.tile(np.nanmean(JD[oktimeg]),lenn)
zi = zi
rain = rain
Ss = np.tile(np.nanmean(sal_sb37_mean[oktimeg]),lenn) 
cp = np.tile(np.nanmean(wave_dominant_period[oktimeg]),lenn) 
sigH = np.tile(np.nanmean(wave_significant_height[oktimeg]),lenn) 
'''

lenn = 1000

#u = np.tile(10,lenn)
u = np.linspace(0,60,lenn)
zu = zu
t = np.tile(28,lenn) 
zt = zt
rh = np.tile(82,lenn) 
zq = zq
P = np.tile(1012,lenn) 
ts = np.tile(29,lenn) 
sw_dn = np.tile(200,lenn)
lw_dn = np.tile(370,lenn)
lat = np.tile(18.9,lenn)
lon = np.tile(-53.6,lenn)
jd = np.tile(250,lenn)
zi = zi
rain = rain
Ss = np.tile(36.8,lenn) 
cp = np.tile(10,lenn) 
sigH = np.tile(2.8,lenn) 

A = coare36vn_zrf_et(u, zu, t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon, jd, zi, rain, Ss, cp, sigH, zrf_u=10.0, zrf_t=10.0, zrf_q=10.0)

#timeS = times[oktimeg]
#timestampS = timestamps[oktimeg]
timeS = times[0:lenn]
timestampS = timestamps[0:lenn]

tauS = A[:,1] # Wind stress (N/m^2)
shtfluxS = A[:,2] # Sensible heat flux (W/m^2)
lhtfluxS = A[:,3] # Latent heat flux (W/m^2)
CdS = A[:,12] # Drag coefficient at height zu (unitless)
ChS = A[:,13] # Sensible heat transfer coefficient at height zu (unitless)
CeS = A[:,14] # Latent heat transfer coefficient at height zu (unitless)

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

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(u,CdS,'.-',color='blue',label='Cd')
plt.legend(loc='lower right')
plt.title('  ',fontsize=18)
plt.xlabel('Wind speed (m/s)')

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(u,ChS,'.-',color='lightcoral',label='Ch')
plt.plot(u,CeS,'.-',color='lightseagreen',label='Cl')
plt.legend(loc='lower right')
plt.title('  ',fontsize=18)
plt.xlabel('Wind speed (m/s)')


