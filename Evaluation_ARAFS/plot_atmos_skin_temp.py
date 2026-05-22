#!/usr/bin/env python3

"""This script is to plot out HAFS atmospheric skin temperature"""

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
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)

print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat_atm = np.asarray(grb.select(shortName='NLAT')[0].data)
lon_atm = np.asarray(grb.select(shortName='ELON')[0].data)

print('Extracting LHTFL at surface')

#grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')
level = 'surface'
tmp = grb.select(shortName='TMP',level=level)[0].data
tmp = tmp - 273.15 # convert K to degC

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

print('Plotting skin temp.')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

#cflevels = np.arange(0,801,50)
#cflevels = np.arange(0,1001,50)
cflevels = np.arange(-2,36,2)
cf = plt.contourf(lon_atm,lat_atm,tmp,levels=cflevels,cmap='turbo',extend='both',alpha=1,transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.03, extendrect=True,shrink=0.85)

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

title = 'Skin Temperature ($^oC$) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(tmp),1)) + ' ($^oC$) ' + ' Max = ' + str(np.round(np.nanmax(tmp),1)) + ' ($^oC$)'+ ' Min = ' + str(np.round(np.nanmin(tmp),1)) + ' ($^oC$)',transform=ax.transAxes)

#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
