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

xlim = conf['xlim']
ylim = conf['ylim']

temp_lim = conf['temp_lim']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

################################################################
# SST from ocean
fname_ocn = conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'
ncfile = os.path.join(conf['COMarafs'], fname_ocn)
nc = xr.open_dataset(ncfile)

var_ocn = np.asarray(nc['SST'][0,:,:])
lon_ocn = np.asarray(nc.geolon)
lat_ocn = np.asarray(nc.geolat)

#================================================================

central_longitude = 180
#################################################################
var_name= 'sst'
units = '($^oC$)'

# create figure and axes instances
#fig = plt.figure(figsize=(8,4))
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
ax.axis('scaled')

cflevels = np.arange(temp_lim[0],temp_lim[1]+1,1)
cmap = plt.get_cmap('turbo')
cf = ax.contourf(lon_ocn, lat_ocn, var_ocn, levels=cflevels[::2], cmap=cmap, extend='both', transform=ccrs.PlateCarree())
lb = ax.contour(lon_ocn, lat_ocn, var_ocn, levels=[26], colors='grey', alpha=0.7,transform=ccrs.PlateCarree(),linewidths=0.5)
lb = ax.contour(lon_ocn, lat_ocn, var_ocn, levels=[0], colors='grey', alpha=1,transform=ccrs.PlateCarree(),linewidths=1)
lb = ax.contour(lon_ocn, lat_ocn, var_ocn, levels=[-1], colors='k', alpha=1,transform=ccrs.PlateCarree(),linewidths=1)
ax.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, shrink=0.8, extendrect=True, ticks=cflevels[::2])
cb.ax.tick_params(labelsize=8)

ax.set_extent([xlim[0], xlim[1], ylim[0], ylim[1]], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

title = 'Sea Surface Temperature (${^o}$C) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

