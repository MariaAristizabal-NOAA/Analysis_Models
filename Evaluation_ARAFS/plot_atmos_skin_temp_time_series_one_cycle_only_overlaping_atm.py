#!/usr/bin/env python3

#exp_names = ['HFSAv2a_baseline_latest','HKPP','HKP2']
#exp_labels = ['HFSAv2a_baseline','HKPP','HKP2']
#exp_colors = ['orange','green','dodgerblue']

"""This script is to plot out HAFS atmos skin temp time series"""

import os
import sys
import logging
import math
import datetime

import yaml
import numpy as np
import pandas as pd
#from scipy.ndimage import gaussian_filter

import grib2io
from netCDF4 import Dataset
import xarray as xr

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

from scipy.interpolate import griddata

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

xlim = conf['xlim']
ylim = conf['ylim']

####################################################################
# Need to find overlapping areas for atm.
COMmodel_ocean = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/'
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(COMmodel_ocean+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

# Read ocean lat and lon and land mask
ocean_name = conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'
ocean_file = os.path.join(COMmodel_ocean+conf['ymdh']+'/00E/', ocean_name)
ocean = xr.open_dataset(ocean_file)
lon_ocean = np.asarray(ocean['geolon'])
lat_ocean = np.asarray(ocean['geolat'])
wet_ocean = np.asarray(ocean['wet_c'])[1:,1:]

# interpolate lon_ocean onto lon atmosphere
# --- Original grid ---
lon_orig = lon_ocean + 180
lat_orig = lat_ocean

# Example data on original grid
#data = np.ones((lon_ocean.shape[0],lon_ocean.shape[1]))
wet_ocean[wet_ocean==0] = np.nan
data = wet_ocean

# --- Target grid ---
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
lon_atm = np.copy(lon)
lat_atm = np.copy(lat)
if abs(np.max(lon_atm) - 360.) < 10.:
    lon_atm[lon_atm>180] = lon_atm[lon_atm>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon_atm = lon_atm - lon_offset

lon_targ = lon_atm
lat_targ = lat_atm

# --- Interpolate ---
points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
values = data.ravel()
data_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,3)
skt_mean = np.empty((len(conf['COMmodels']),len(ff)))
skt_mean[:] = np.nan
skt_max = np.empty((len(conf['COMmodels']),len(ff)))
skt_max[:] = np.nan
skt_min = np.empty((len(conf['COMmodels']),len(ff)))
skt_min[:] = np.nan
skt_std = np.empty((len(conf['COMmodels']),len(ff)))
skt_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    #================================================================
    # Read fluxes
    for f,fh in enumerate(ff):

        if len(str(fh))==3:
            fhour = str(fh)
        if len(str(fh))==2:
            fhour = '0'+str(fh)
        if len(str(fh))==1:
            fhour = '00'+str(fh)

        fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.f'+fhour+'.grb2'
        grib2file = os.path.join(conf['COMmodels'][exp]+conf['ymdh']+'/00E', fname)

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
    
        print('Extracting SHTFL at surface')
        level = 'surface'
        skt = grb.select(shortName='TMP',level=level)[0].data
        skt = skt - 273.15
        skt = skt * data_interp

        #oklon = np.logical_and(lon>=xlim[0],lon<xlim[1]) 
        #oklat = np.logical_and(lat>=ylim[0],lat<ylim[1]) 

        #lonsub = lon[oklon]
        #lonsub = lon[oklon]
    
        skt_mean[exp,f] = np.nanmean(skt)
        skt_max[exp,f] = np.nanmax(skt)
        skt_min[exp,f] = np.nanmin(skt)
        skt_std[exp,f] = np.nanstd(skt)

########
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,skt_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Skin Temp. ($^oC$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Skin Temperature',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,skt_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,skt_min[exp,:],skt_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Domain Mean Skin Temp. ($^oC$)',fontsize=14,labelpad=10)
plt.title('Skin Temperature',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))


