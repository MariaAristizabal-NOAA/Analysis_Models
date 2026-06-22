#!/usr/bin/env python3

"""This script is to plot out HAFS SST time series"""

import os
import sys
import logging
import math
import datetime

import yaml
import numpy as np
import pandas as pd

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
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

xlim = conf['xlim']
ylim = conf['ylim']

####################################################################
# Need to find overlapping areas for atm.
fname = conf['stormModel'].lower()+'.'+conf['ymdhs'][0]+'.'+conf['fhhhs'][0]+'.grb2'
grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdhs'][0]+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

levstr='10 m above ground'
ugrd = grb.select(shortName='UGRD', level=levstr)[0].data
atm_domain = np.copy(ugrd)
atm_domain[np.isfinite(atm_domain)] = 1

# Read ocean lat and lon and land mask
ocean_name = conf['ymdhs'][0]+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhhs'][0]+'.nc'
ocean_file = os.path.join(conf['COMmodels'][0]+conf['ymdhs'][0]+'/00E/', ocean_name)
ocean = xr.open_dataset(ocean_file)
lon_ocean = np.asarray(ocean['geolon'])
lat_ocean = np.asarray(ocean['geolat'])
#wet_ocean = np.asarray(ocean['wet_c'])[1:,1:]

# data on original grid
data = np.ones((lon.shape[0],lon.shape[1]))

# --- original grid ---
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

lon_orig = lon_atm
lat_orig = lat_atm

# interpolate lon_atm onto lon_ocean
# --- Original grid ---
lon_targ = lon_ocean + 180
lat_targ = lat_ocean

# --- Interpolate ---
points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
values = data.ravel()
data_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,3)
sst_mean = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
sst_mean[:] = np.nan
sst_max = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
sst_max[:] = np.nan
sst_min = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
sst_min[:] = np.nan
sst_std = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
sst_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 
    print(conf['COMmodels'][exp])
    for c,ymdh in enumerate(conf['ymdhs']):
        print(ymdh)
        for f,fh in enumerate(ff):
    
            if len(str(fh))==3:
                fhour = str(fh)
            if len(str(fh))==2:
                fhour = '0'+str(fh)
            if len(str(fh))==1:
                fhour = '00'+str(fh)
    
            fname = ymdh+'.'+conf['stormModel'].lower()+'.mom6.'+'f'+fhour+'.nc'
            ncfile = os.path.join(conf['COMmodels'][exp]+ymdh+'/00E/', fname)
            try: 
                nc = xr.open_dataset(ncfile)
                sst = np.asarray(nc['SST'][0,:,:])
                sst_atm_domain = sst * data_interp
                sst_mean[exp,c,f] = np.nanmean(sst_atm_domain)
                sst_max[exp,c,f] = np.nanmax(sst_atm_domain)
                sst_min[exp,c,f] = np.nanmin(sst_atm_domain)
                sst_std[exp,c,f] = np.nanstd(sst_atm_domain)
            except FileNotFoundError:
                print("The file does not exist.")
    
########
for exp in np.arange(len(conf['COMmodels'])):
    fig,ax = plt.subplots(figsize = (9,6))
    if np.isfinite(np.mean(sst_mean[exp,:,:])):
        for c,ymdh in enumerate(conf['ymdhs']):
            if c%5 == 0:
                plt.plot(ff,sst_mean[exp,c,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdhs'][c][0:-4],markeredgecolor='k',markersize=7,alpha=0.5)
            else:
                plt.plot(ff,sst_mean[exp,c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=7,alpha=0.5)
ax.legend(loc='upper left',bbox_to_anchor=(0.95,1.1))
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Sea Surface Temperature',fontsize=14)
plt.xlim([0,121])
plt.xticks(np.arange(0,120,12))

########
for exp in np.arange(len(conf['COMmodels'])):
    fig,ax = plt.subplots(figsize = (9,6))
    if np.isfinite(np.mean(sst_mean[exp,:,:])):
        for c,ymdh in enumerate(conf['ymdhs']):
            if c%5 == 0:
                plt.plot(ff,sst_mean[exp,c,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdhs'][c],markeredgecolor='k',markersize=7)

ax.legend(loc='upper left',bbox_to_anchor=(0.90,1.0))
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Sea Surface Temperature',fontsize=14)
plt.xlim([0,121])
plt.xticks(np.arange(0,120,12))


