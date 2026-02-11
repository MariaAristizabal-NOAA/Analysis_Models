#!/usr/bin/env python3

"""This scrip plots the sea surface height time series for a specific location. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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
print('Parse the config file: plot_ocean1.yml:')
with open('plot_ocean1.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
#conf['fhour'] = int(conf['fhhh'][1:])
#conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
#conf['validTime'] = conf['initTime'] + conf['fcstTime']

experiments = conf['experiments']
folder_exp = conf['COMarafs']
file_obs = conf['FileData']
lon_station = conf['lon_station']
lat_station = conf['lat_station']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

#================================================================
# Read data from tide gauge
with open(file_obs, "r") as f:
    lines = f.readlines()

time_obs = []
timestamp_obs = []
ssh_obs = []
for line in lines[1:]:
    print(line)
    dd = line.split('\t')[0]
    tt = line.split('\t')[1]
    date_str = dd + ' ' + tt
    ttime = datetime.strptime(date_str, "%m/%d/%Y %H:%M")
    time_obs.append(ttime)
    timestamp_obs.append(mdates.date2num(ttime))
    ssh_obs.append(float(line.split('\t')[4].split('\n')[0]))
        
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

plt.figure(figsize=(10,5))
plt.plot(time_obs,ssh_obs,'.-',label='Tidal Gauge')
plt.legend()
plt.title('Arena Cove CA, lon = ' + str(lon_station) + ', lat = ' + str(lat_station))
plt.ylabel('Water Level (m) (with respect to MSL)')

plt.figure(figsize=(10,5))
for ex,exp in enumerate(experiments):
    plt.plot(target_time,target_ssh[ex,:],'.-',label=exp)
plt.legend()
plt.title('Arena Cove CA, lon = ' + str(lon_station) + ', lat = ' + str(lat_station))
plt.ylabel('SSH (m)')


#================================================================
# create figure and axes instances
fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.axis('scaled')
ax.plot(lon_station,lat_station,'*r',markersize=5)
# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
plt.xlim([-130,-110])
plt.ylim([10,60])

#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
#plt.close("all")
