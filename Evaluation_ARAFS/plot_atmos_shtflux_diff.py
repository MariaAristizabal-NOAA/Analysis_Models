#!/usr/bin/env python3

"""This script is to plot out HAFS atmospheric sensible heat flux and 10-m wind."""

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

import matplotlib.pyplot as plt

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.ndimage import gaussian_filter

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Read FV3 files
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdh']+'/00E/', fname)

print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat_atm = np.asarray(grb.select(shortName='NLAT')[0].data)
lon_atm = np.asarray(grb.select(shortName='ELON')[0].data)

print('Extracting SHTFL at surface')

grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')
shf1 = grb.select(shortName='SHTFL')[0].data

grib2file = os.path.join(conf['COMmodels'][1]+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')
shf2 = grb.select(shortName='SHTFL')[0].data

diff = gaussian_filter(shf2 - shf1,10)
#===================================================================================================
lon = lon_atm
lat = lat_atm

print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

lon_atm = lon

#===================================================================================================
print('Plotting Sensible Heat Flux and 10-m wind')

var_name= 'SHTFL diff'
units = '($W/m^2$)'

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-50,51,5)
cf = plt.contourf(lon_atm,lat_atm,diff,levels=cflevels,cmap='coolwarm',extend='both',alpha=1,transform=transform)
#cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.8, extendrect=True)

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

title = 'Sensible Heat Flux Differences (W/m$^2$) \n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center') #, x=1.2,y

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(diff),1)) + ' (W/m$^2$) \n' + 'Max = ' + str(np.round(np.nanmax(diff),1)) + ' (W/m$^2$)' + ' Min = ' + str(np.round(np.nanmin(diff),1)) + ' (W/m$^2$)',transform=ax.transAxes)

#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
