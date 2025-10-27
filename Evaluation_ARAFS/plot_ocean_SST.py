#!/usr/bin/env python3

"""This script is to plot ARAFS ocean SST"""

import os

import yaml
import numpy as np
import pandas as pd
import xarray as xr

import grib2io

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Parse the yaml config file
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

################################################################
# SST from ocean
fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'
ncfile = os.path.join(conf['COMarafs'], fname_ocn)
nc = xr.open_dataset(ncfile)

var_ocn = np.asarray(nc['SST'][0,:,:])
lon_ocn = np.asarray(nc.xh)
lat_ocn = np.asarray(nc.yh)

lonmin_raw = np.min(lon_ocn)
lonmax_raw = np.max(lon_ocn)
print('raw lonlat limit: ', np.min(lon_ocn), np.max(lon_ocn), np.min(lat_ocn), np.max(lat_ocn))

#================================================================
# Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
lon_ocn[lon_ocn>180] = lon_ocn[lon_ocn>180] - 360
lon_ocn[lon_ocn<-180] = lon_ocn[lon_ocn<-180] + 360
sort_lon_ocn = np.argsort(lon_ocn)
lon_ocn = lon_ocn[sort_lon_ocn]

# define grid boundaries
lonmin_ocn = np.min(lon_ocn)
lonmax_ocn = np.max(lon_ocn)
latmin_ocn = np.min(lat_ocn)
latmax_ocn = np.max(lat_ocn)
print('new lonlat limit: ', np.min(lon_ocn), np.max(lon_ocn), np.min(lat_ocn), np.max(lat_ocn))

# Shift central longitude so the Southern Hemisphere and North Indin Ocean domains are plotted continuously
if np.logical_and(lonmax_ocn >= 90, lonmax_ocn <=180):
    central_longitude = 90
else:
    central_longitude = -90
print('central longitude: ',central_longitude)
#central_longitude = -90

# sort var according to the new longitude
var_ocn = var_ocn[:,sort_lon_ocn]

#################################################################
var_name= 'sst'
units = '($^oC$)'

# create figure and axes instances
fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
ax.axis('scaled')

cflevels = np.arange(8,29,2)
cmap = plt.get_cmap('turbo')
cf = ax.contourf(lon_ocn, lat_ocn, var_ocn, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
lb = ax.contour(lon_ocn, lat_ocn, var_ocn, levels=[26], colors='grey', alpha=0.7,transform=ccrs.PlateCarree(),linewidths=0.5)
ax.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=20, shrink=0.6, extendrect=True, ticks=cflevels)
cb.ax.tick_params(labelsize=8)

ax.set_extent([lonmin_raw, lonmax_raw, latmin_ocn, latmax_ocn], crs=ccrs.PlateCarree())

# Add gridlines and labels
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 20))
gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 10))
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Sea Surface Temperature (${^o}$C)'
title_left = """{0}\n{1}""".format(model_info,var_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
footer = os.environ.get('FOOTERgraph','Experimental ARAFS Product').strip()
ax.text(1.0,-0.1, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)


