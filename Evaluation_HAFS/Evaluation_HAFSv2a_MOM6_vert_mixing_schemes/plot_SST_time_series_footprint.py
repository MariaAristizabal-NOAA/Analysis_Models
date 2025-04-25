#!/usr/bin/env python3

"""This script is to plot SST footprint time series."""

import os
import sys
import glob

import yaml
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath

from geo4HYCOM import haversine

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

#================================================================
def latlon_str2num(string):
    """Convert lat/lon string into numbers."""
    value = pd.to_numeric(string[:-1], errors='coerce') / 10
    if string.endswith(('N', 'E')):
        return value
    else:
        return -value

#================================================================
def get_adeck_track(adeck_file):

    cols = ['basin', 'number', 'ymdh', 'technum', 'tech', 'tau', 'lat', 'lon', 'vmax', 'mslp', 'type','rad', 'windcode', 'rad1', 'rad2', 'rad3', 'rad4', 'pouter', 'router', 'rmw', 'gusts', 'eye','subregion', 'maxseas', 'initials', 'dir', 'speed', 'stormname', 'depth','seas', 'seascode', 'seas1', 'seas2', 'seas3', 'seas4', 'userdefined','userdata1', 'userdata2', 'userdata3', 'userdata4', 'userdata5', 'userdata6', 'userdata7', 'userdata8','userdata9', 'userdata10', 'userdata11', 'userdata12', 'userdata13', 'userdata14', 'userdata15', 'userdata16']

    # Read in adeckFile as pandas dataFrame
    print('Read in adeckFile ...')
    adeck = pd.read_csv(adeck_file, index_col=False, names=cols, dtype=str, header=None, skipinitialspace=True)

    adeck['lat'] = adeck['lat'].apply(latlon_str2num)
    adeck['lon'] = adeck['lon'].apply(latlon_str2num)
    #adeck['vmax'] = adeck['vmax'].apply(lambda x: x if x>0 else np.nan)
    #adeck['mslp'] = adeck['mslp'].apply(lambda x: x if x>800 else np.nan)
    adeck['init_time'] = pd.to_datetime(adeck['ymdh'], format='%Y%m%d%H', errors='coerce')
    adeck['valid_time'] = pd.to_timedelta(adeck['tau'].apply(pd.to_numeric, errors='coerce'), unit='h') + adeck['init_time']

    fhour, ind = np.unique(adeck['tau'],return_index=True)
    lat_adeck = np.asarray(adeck['lat'][ind])
    lon_adeck = np.asarray(adeck['lon'][ind])
    init_time = np.asarray(adeck['init_time'][ind])
    valid_time = np.asarray(adeck['valid_time'][ind])

    return fhour,lat_adeck,lon_adeck,init_time,valid_time

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_ocean3.yml:')
with open('plot_ocean3.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

#================================================================
# Cycle the different experiments

ff = np.arange(0,127,3)
sst_mean = np.empty((len(conf['COMhafs']),len(ff)))
sst_mean[:] = np.nan
sst_max = np.empty((len(conf['COMhafs']),len(ff)))
sst_max[:] = np.nan
sst_min = np.empty((len(conf['COMhafs']),len(ff)))
sst_min[:] = np.nan
dsst_mean = np.empty((len(conf['COMhafs']),len(ff)))
dsst_mean[:] = np.nan
dsst_max = np.empty((len(conf['COMhafs']),len(ff)))
dsst_max[:] = np.nan
dsst_min = np.empty((len(conf['COMhafs']),len(ff)))
dsst_min[:] = np.nan

for exp in np.arange(len(conf['COMhafs'])): 

    # Get lat and lon from adeck file
    adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
    adeck_file = os.path.join(conf['COMhafs'][exp],adeck_name)

    fhour,lat_adeck,lon_adeck,init_time,valid_time = get_adeck_track(adeck_file)

    print('lon_adeck = ',lon_adeck)
    print('lat_adeck = ',lat_adeck)

    #================================================================
    # Read SST
    for f,fh in enumerate(ff):
    
        # Get lat and lon from adeck file
        adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
        adeck_file = os.path.join(conf['COMhafs'][exp],adeck_name)
    
        fhour,lat_adeck,lon_adeck,init_time,valid_time = get_adeck_track(adeck_file)
    
        print('lon_adeck = ',lon_adeck)
        print('lat_adeck = ',lat_adeck)

        oceanf = glob.glob(os.path.join(conf['COMhafs'][exp],'*f006.nc'))[0].split('/')[-1].split('.')
        ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

        if f < len(fhour):
            if ocean == 'mom6':
                fname0 =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+'f003.nc'
                fname = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+'f'+fhour[f]+'.nc'
                ncfile0 = os.path.join(conf['COMhafs'][exp], fname0)
                ncfile = os.path.join(conf['COMhafs'][exp], fname)
                if os.path.exists(ncfile):
                    print(f'ncfile: {ncfile}')
                    nc0 = xr.open_dataset(ncfile0)
                    nc = xr.open_dataset(ncfile)
                    sst0 = np.asarray(nc0['SST'][0,:,:])
                    sst = np.asarray(nc['SST'][0,:,:])
                    lon = np.asarray(nc.xh)
                    lat = np.asarray(nc.yh)
                    lon_mesh,lat_mesh = np.meshgrid(lon,lat)
                else:
                    continue
    
            if ocean == 'hycom':
                fname = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+'f'+fhour[f]+'.nc'
                ncfile = os.path.join(conf['COMhafs'][exp], fname)
                if os.path.exists(ncfile):
                    print(f'ncfile: {ncfile}')
                    nc = xr.open_dataset(ncfile)
                    sst = np.asarray(nc['temperature'][0,0,:,:])
                    lon = np.asarray(nc.Longitude)
                    lat = np.asarray(nc.Latitude)
                    lon_mesh,lat_mesh = np.meshgrid(lon,lat)
                else:
                    continue

            # Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
            lon[lon>180] = lon[lon>180] - 360
            #sort_lon = np.argsort(lon)
            #lon = lon[sort_lon]
        
            # define grid boundaries
            lonmin_new = np.min(lon)
            lonmax_new = np.max(lon)
            latmin = np.min(lat)
            latmax = np.max(lat)
            print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
    
            Rkm=200    # search radius [km]
            #lns,lts = np.meshgrid(lon,lat)
            dummy = np.ones(lat_mesh.shape)
        
            dR = haversine(lon_mesh,lat_mesh,lon_adeck[f],lat_adeck[f])/1000.
            dumb = dummy.copy()
            dumb[dR>Rkm] = np.nan

            dsst = sst - sst0

            sstt = sst*dumb
            sst_mean[exp,f] = np.nanmean(sstt)
            sst_max[exp,f] = np.nanmax(sstt)
            sst_min[exp,f] = np.nanmin(sstt)

            dsstt = dsst*dumb
            dsst_mean[exp,f] = np.nanmean(dsstt)
            dsst_max[exp,f] = np.nanmax(dsstt)
            dsst_min[exp,f] = np.nanmin(dsstt)

fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMhafs'])): 
    plt.plot(np.arange(0,127,3),sst_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(np.arange(0,127,3),sst_min[exp,:],sst_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(np.arange(0,127,3),sst_max[exp,:],'.-',color=conf['exp_colors'][exp],markeredgecolor='k',markersize=7)

ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Sea Surface Temperature Footprint (200 km around storm eye))',fontsize=14)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMhafs'])):
    plt.plot(np.arange(0,127,3),dsst_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(np.arange(0,127,3),dsst_min[exp,:],dsst_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(np.arange(0,127,3),dsst_max[exp,:],'.-',color=conf['exp_colors'][exp],markeredgecolor='k',markersize=7)

ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Change in SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Sea Surface Temperature Change Footprint (200 km around storm eye))',fontsize=14)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))


