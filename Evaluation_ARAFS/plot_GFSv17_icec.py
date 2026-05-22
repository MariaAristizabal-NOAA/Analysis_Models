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
from scipy.ndimage import gaussian_filter

import grib2io
from netCDF4 import Dataset

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyproj
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import (LongitudeLocator, LongitudeFormatter, LatitudeLocator, LatitudeFormatter)

# Parse the yaml config file
print('Parse the config file: plot_GFSv17.yml:')
with open('plot_GFSv17.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

xlim = conf['xlim']
ylim = conf['ylim']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

grib2file = conf['grib2file']

print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

#print('Extracting lat, lon')
#lat = np.asarray(grb.select(shortName='NLAT')[0].data)
#lon = np.asarray(grb.select(shortName='ELON')[0].data)
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180

'''
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
[nlat, nlon] = np.shape(lon)
'''

print('Extracting SHTFL at surface')
icec = grb.select(shortName='ICEC')[0].data
#shf = gaussian_filter(shf, 5)

print('Extracting lat and lon')
lat, lon =  grb.select(shortName='ICEC')[0].latlons()
#===================================================================================================
print('Plotting ice coverage')
#fig_prefix = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']

# Set default figure parameters

mpl.rcParams['figure.figsize'] = [8, 8]
mpl.rcParams["figure.dpi"] = 150
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['legend.fontsize'] = 8

mpl.rcParams['figure.figsize'] = [8, 5.4]
#fig_name = fig_prefix+'.'+'shtflux_wind10m.'+conf['fhhh'].lower()+'.png'
cbshrink = 1.0
lonmin = np.min(lon)
lonmax = np.max(lon)
lonint = 10.0
latmin = np.min(lat)
latmax = np.max(lat)
latint = 10.0

'''
plt.figure()
plt.contourf(lon-180,lat,shf)
plt.colorbar()
'''

lon_offset = 180
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

cflevels = np.arange(0,1.1,0.1)
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')
cf = ax.contourf(lon-180, lat, icec, cflevels, cmap='Blues_r', transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=40, shrink=1, extendrect=True, ticks=cflevels[::2])

plt.xlim([-20,100])
plt.ylim([-20,85])

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Ice Coverage (shaded)'
title_left = """{0}\n{1}""".format(model_info,var_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)

#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
