#!/usr/bin/env python3

"""This script is to plot out HAFS atmospheric MSLP and 10-m wind."""

import os

import yaml
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import grib2io

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_ocean.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Read FV3 files
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E', fname)
print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat_atm = grb.select(shortName='NLAT')[0].data
lon_atm = grb.select(shortName='ELON')[0].data

print('Extracting UGRD, VGRD at 10 m above ground')
levstr='10 m above ground'
ugrd = grb.select(shortName='UGRD', level=levstr)[0].data
ugrd = ugrd * 1.94384 # convert m/s to kt

vgrd = grb.select(shortName='VGRD', level=levstr)[0].data
vgrd = vgrd * 1.94384 # convert m/s to kt

# Calculate wind speed
wspd = (ugrd**2+vgrd**2)**.5

print('Extracting MSLET')
slp = grb.select(shortName='MSLET')[0].data
slp = slp * 0.01 # convert Pa to hPa
#slp = gaussian_filter(slp, 5)

#================================================================
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
print('raw lonlat limit: ', np.min(lon_atm), np.max(lon_atm), np.min(lat_atm), np.max(lat_atm))
if abs(np.max(lon_atm) - 360.) < 10.:
    lon_atm[lon_atm>180] = lon_atm[lon_atm>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon_atm = lon_atm - lon_offset
print('new lonlat limit: ', np.min(lon_atm), np.max(lon_atm), np.min(lat_atm), np.max(lat_atm))

#================================================================
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

#cflevels = np.arange(0,51,5)
cslevels = np.arange(940,1040,5)
lblevels = np.arange(940,1040,10)
cmap = plt.get_cmap('rainbow')

cf = ax.contourf(lon_atm,lat_atm,slp,levels=cslevels,cmap=cmap,extend='both',alpha=0.8,transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, shrink=0.8,extendrect=True,ticks=lblevels)

cs = ax.contour(lon_atm,lat_atm,slp,levels=cslevels,cmap=cmap,extend='both',alpha=0.8,transform=transform)
#cs = ax.contour(lon_atm,lat_atm,slp,levels=cslevels,colors='k',extend='both',alpha=0.8,transform=transform)
#lb = plt.clabel(cs, levels=lblevels, inline_spacing=1, fmt='%d', fontsize=8)

ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

title = 'Sea Level Pressure (hPa) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')


