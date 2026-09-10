#!/usr/bin/env python3

"""This scrip plots the sea surface height time series for a specific location. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
import requests
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.interpolate import interp1d
from scipy.signal import find_peaks

#================================================================
def get_var_from_model_following_trajectory(files_model,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs):

    # depth_level is an optional argument
    # kwargs = dict(depth_level = 0)

    depth_level = kwargs.get('depth_level', None)

    target_var = np.empty((len(files_model)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        model = xr.open_dataset(file,engine="netcdf4")
        if file.split('.')[-1] == 'nc':
            t = model[time_name][:]
            timestamp = mdates.date2num(t)[0]
            target_time.append(mdates.num2date(timestamp))
            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])

        # Interpolating lat_obs and lon_obs into model grid
        sublon = np.interp(timestamp,timestamp_obs,lon_obs)
        sublat = np.interp(timestamp,timestamp_obs,lat_obs)
        if lon_model.ndim == 1:
            oksort_lon = np.argsort(lon_model)
            oksort_lat = np.argsort(lat_model)
            lonn = lon_model[oksort_lon]
            latt = lat_model[oksort_lat]
        if lon_model.ndim == 2:
            oksort_lon = np.argsort(lon_model[0,:])
            oksort_lat = np.argsort(lat_model[:,0])
            lonn = lon_model[0,oksort_lon]
            latt = lat_model[oksort_lat,0]
        oklon = int(np.round(np.interp(sublon,lonn,np.arange(len(lonn)))))
        oklat = int(np.round(np.interp(sublat,latt,np.arange(len(latt)))))
        if file.split('.')[-1] == 'nc':
            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name][0,depth_level,oksort_lat,oksort_lon])
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name][0,oksort_lat,oksort_lon])
            
            if np.isfinite(var[oklat,oklon]): 
                target_var[x] = var[oklat,oklon]
            else:
                target_var[x] = var[oklat-1,oklon]

    return target_time, target_var

#================================================================
# Parse the yaml config file
print('Parse the config file: Tide_gauge_ssh_vs_ARAFS_time_series.yml:')
with open('Tide_gauge_ssh_vs_ARAFS_time_series.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
#conf['fhour'] = int(conf['fhhh'][1:])
#conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
#conf['validTime'] = conf['initTime'] + conf['fcstTime']

experiments = conf['experiments']
folder_exp = conf['COMarafs']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

#================================================================
# Charleston, OR
station_id = "9432780"

url = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station_id}.json"

response = requests.get(url)
data = response.json()["stations"][0]

lat_station = data["lat"]
lon_station = data["lng"]
name_station = data["name"]

# Setup NOAA CO-OPS API URL
url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

params = {
    "begin_date": "20230101",
    "end_date": "20230228",
    "station": station_id,
    "product": "hourly_height",  # Options: hourly_height, water_level (6-min), high_low
    "datum": "MSL",  # Mean Lower Low Water (Standard tidal datum)
    "units": "english",  # 'english' (feet) or 'metric' (meters)
    "time_zone": "gmt",  # 'gmt', 'lst' (local standard), or 'lst_ldt' (local day light)
    "format": "json",
    "application": "PythonScript",
}

response = requests.get(url, params=params)
data = response.json()

if "data" in data:
    df = pd.DataFrame(data["data"])
else:
    print("Error or no data returned:", data)

time_obs = pd.to_datetime(df['t'])
timestamp_obs = time_obs
msl_obs = df['v'].astype(float)*0.3048

#================================================================
# Read ocean files
oceanf = glob.glob(os.path.join(folder_exp[0],'*f006.nc'))[0].split('/')[-1].split('.')
ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

if ocean == 'mom6':
    files_model = sorted(glob.glob(os.path.join(folder_exp[0],'*mom6*.nc')))

if ocean == 'hycom':
    files_model = sorted(glob.glob(os.path.join(folder_exp[0],'*hycom.2d*.nc')))

target_ssh = np.empty((len(experiments),len(files_model)))
target_ssh[:] = np.nan

for f,folder in enumerate(folder_exp):
    print(folder)

    oceanf = glob.glob(os.path.join(folder,'*f006.nc'))[0].split('/')[-1].split('.')
    ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

    if ocean == 'mom6':
        files_model = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
        var_name = 'SSH'
        lon_name = 'xh'
        lat_name = 'yh'
        time_name = 'time'
    
    if ocean == 'hycom':
        files_model = sorted(glob.glob(os.path.join(folder,'*hycom.2d*.nc')))
        var_name = 'sea_surface_height'
        lon_name = 'Longitude'
        lat_name = 'Latitude'
        time_name = 'Time'
    
    depth_level = 0
    timestamp_obs = timestamp_obs
    lon_obs = np.tile(lon_station,len(timestamp_obs))
    lat_obs = np.tile(lat_station,len(timestamp_obs))
    kwargs = dict(depth_level = 0)
    
    #oklo = np.isfinite(lon_obss)
    #lon_obs = lonB[oklo]
    #lat_obs = latB[oklo]
    #timestamp_obs = timestamp_obss[oklo]
    
    target_time, target_ssh[f,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

#================================================================
# Find peak high tides and low tides
ssh_obs = np.asarray(msl_obs)
time_obs = np.asarray(time_obs)
timestamp = time_obs

high_peaks, _ = find_peaks(ssh_obs, distance=20) # ~10 hours apart for M2
low_peaks, _ = find_peaks(-ssh_obs, distance=20)

# Interpolate smooth curves through peaks
f_upper = interp1d(timestamp[high_peaks], ssh_obs[high_peaks], kind='cubic', fill_value="extrapolate")
f_lower = interp1d(timestamp[low_peaks], ssh_obs[low_peaks], kind='cubic', fill_value="extrapolate")

spline_upper_env = f_upper(timestamp)
spline_lower_env = f_lower(timestamp)

#================================================================
plt.figure(figsize=(10,5))
plt.plot(time_obs,ssh_obs,'.-',label='Tidal Gauge')
plt.plot(time_obs,spline_upper_env,'.-',label='Upper Envelope')
plt.legend()
plt.title(name_station + ', lon = ' + str(lon_station) + ', lat = ' + str(lat_station))
plt.ylabel('Water Level (m) (with respect to MSL)')

#================================================================
plt.figure(figsize=(10,5))
for ex,exp in enumerate(experiments):
    plt.plot(target_time,target_ssh[ex,:],'o-',markeredgecolor='k',label=exp)
plt.ylabel('SSH (m)')
plt.title(name_station + ', lon = ' + str(lon_station) + ', lat = ' + str(lat_station))
plt.legend()

#================================================================
fig, ax1 = plt.subplots(figsize=(10, 5))

# Plot 1st line on the primary y-axis (left)
color1 = 'tab:orange'
ax1.set_xlabel('')
ax1.set_ylabel('Water Level (m) (with respect to MSL)', color=color1)
line1 = ax1.plot(time_obs,spline_upper_env,color=color1, linewidth=2, label='Upper Envelope')
ax1.tick_params(axis='y', labelcolor=color1)
#ax1.set_ylim([0.8,1.4])

# Create twin y-axis on the right
ax2 = ax1.twinx()  

# Plot 2nd line on the secondary y-axis (right)
color2 = 'tab:purple'
ax2.set_ylabel('SSH (m)', color=color2)
for ex,exp in enumerate(experiments):
    line2 = ax2.plot(target_time,target_ssh[ex,:],color=color2, linewidth=2, linestyle='-.',label=exp)
ax2.tick_params(axis='y', labelcolor=color2)
#ax2.set_ylim([0.8,1.4])

# Combined Legend for both lines
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left')
plt.title(name_station + ', lon = ' + str(lon_station) + ', lat = ' + str(lat_station))
plt.grid(True, alpha=0.3)
plt.show()

start_date = datetime(2023, 1, 9)
end_date = datetime(2023, 1, 17)
ax2.set_xlim(start_date, end_date)

#================================================================
# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.axis('scaled')
ax.plot(lon_station,lat_station,'*r',markersize=5)
# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
plt.xlim([-130,-110])
plt.ylim([10,60])
plt.title(name_station,fontsize=16)
