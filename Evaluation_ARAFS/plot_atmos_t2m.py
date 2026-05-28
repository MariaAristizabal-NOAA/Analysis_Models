#!/usr/bin/env python3

"""This script is to plot out HAFS atmospheric surface temperature, MSLP and 10-m wind."""

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

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)

print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
[nlat, nlon] = np.shape(lon)

print('Extracting Temperature at 2 m above ground')
levstr='2 m above ground'
tmp = grb.select(shortName='TMP', level=levstr)[0].data
tmp = tmp - 273.15 # convert K to degC
#tmp = gaussian_filter(tmp, 2)

#===================================================================================================
print('Plotting 2 m temperature, MSLET and 10 m wind')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
#ax.axis('equal')
ax.axis('scaled')

#cflevels = np.linspace(-20, 40, 121)
cflevels = np.arange(-24,26,4)
#ctmp = plt.get_cmap('nipy_spectral')
#cmap = mpl.colors.LinearSegmentedColormap.from_list('sub_'+ctmp.name,ctmp(np.linspace(0.04, 0.98, 201)))
cf = ax.contourf(lon, lat, tmp, levels=cflevels, cmap='rainbow', extend='both', transform=transform)
ax.contour(lon, lat, tmp, levels=[0], colors='k', transform=transform)
ax.contour(lon, lat, tmp, levels=[-2], colors='grey', transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True, ticks=cflevels)

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

# Add borders and coastlines
#ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='whitesmoke')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

title = '2 m Temperature (${^o}$C) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
