#!/usr/bin/env python3

"""This script is to plot out HAFS vertical wind shear foot print time series."""

import os
import sys
import logging
import math
import datetime

import yaml
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import grib2io
from netCDF4 import Dataset

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

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
print('Parse the config file: plot_atmos2.yml:')
with open('plot_atmos2.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

#================================================================
# Cycle the different experiments

ff = np.arange(0,127,3)
vws_mean = np.empty((len(conf['COMhafs']),len(ff)))
vws_mean[:] = np.nan
vws_max = np.empty((len(conf['COMhafs']),len(ff)))
vws_max[:] = np.nan
vws_min = np.empty((len(conf['COMhafs']),len(ff)))
vws_min[:] = np.nan

for exp in np.arange(len(conf['COMhafs'])): 

    # Get lat and lon from adeck file
    adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
    adeck_file = os.path.join(conf['COMhafs'][exp],adeck_name)

    fhour,lat_adeck,lon_adeck,init_time,valid_time = get_adeck_track(adeck_file)

    print('lon_adeck = ',lon_adeck)
    print('lat_adeck = ',lat_adeck)

    #================================================================
    # Read fluxes
    for f,fh in enumerate(ff):
    
        # Get lat and lon from adeck file
        adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
        adeck_file = os.path.join(conf['COMhafs'][exp],adeck_name)
    
        fhour,lat_adeck,lon_adeck,init_time,valid_time = get_adeck_track(adeck_file)
    
        print('lon_adeck = ',lon_adeck)
        print('lat_adeck = ',lat_adeck)
       
        fname = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+'f'+fhour[f]+'.grb2'
        grib2file = os.path.join(conf['COMhafs'][exp], fname)
        print(f'grib2file: {grib2file}')
        grb = grib2io.open(grib2file,mode='r')
    
        print('Extracting lat, lon')
        lat = np.asarray(grb.select(shortName='NLAT')[0].data)
        lon = np.asarray(grb.select(shortName='ELON')[0].data)
     
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
    
        #==========================================================
        print('Extracting UGRD, VGRD at 200 hPa')
        levstr='200 mb'
        u200 = grb.select(shortName='UGRD', level=levstr)[0].data
        u200 = np.asarray(u200) * 1.94384 # convert m/s to kt
        
        v200 = grb.select(shortName='VGRD', level=levstr)[0].data
        v200 = np.asarray(v200) * 1.94384 # convert m/s to kt
        
        #==========================================================
        print('Extracting UGRD, VGRD at 850 hPa')
        levstr='850 mb'
        u850 = grb.select(shortName='UGRD', level=levstr)[0].data
        u850 = np.asarray(u850) * 1.94384 # convert m/s to kt
        
        v850 = grb.select(shortName='VGRD', level=levstr)[0].data
        v850 = np.asarray(v850) * 1.94384 # convert m/s to kt
        
        # Calculate VWS
        uvws = u200-u850
        vvws = v200-v850
        mag = (uvws**2+vvws**2)**.5

        #==========================================================
        Rkm=200    # search radius [km]
        #lns,lts = np.meshgrid(lon,lat)
        dummy = np.ones(lat.shape)
    
        dR = haversine(lon,lat,lon_adeck[f],lat_adeck[f])/1000.
        dumb = dummy.copy()
        dumb[dR>Rkm] = np.nan
    
        vws = mag*dumb
        vws_mean[exp,f] = np.nanmean(vws)
        vws_max[exp,f] = np.nanmax(vws)
        vws_min[exp,f] = np.nanmin(vws)

fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMhafs'])): 
    plt.plot(fhour.astype('int'),vws_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(fhour.astype('int'),vws_min[exp,:],vws_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(fhour.astype('int'),vws_max[exp,:],'.-',color=conf['exp_colors'][exp],markeredgecolor='k',markersize=7)

ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Vertical Wind Shear ($Kt.$)',fontsize=14,labelpad=10)
plt.title('Vertical Wind Shear Footprint (200 km around storm eye))',fontsize=14)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

#plt.savefig(fig_name, bbox_inches='tight')

