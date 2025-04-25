#!/usr/bin/env python3

"""This scrip plots fields from the coupler output. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Set Cartopy data_dir location
#cartopy.config['data_dir'] = conf['cartopyDataDir']
#print(conf)

#================================================================
# Read MOM6 file

fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc' 

ncfile = os.path.join(conf['COMhafs'], fname) 
nc = xr.open_dataset(ncfile)

var = np.asarray(nc['SST'][0,:,:])

lon = np.asarray(nc.xh)
lat = np.asarray(nc.yh)
lonmin_raw = np.min(lon)
lonmax_raw = np.max(lon)
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

#================================================================
# Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
lon[lon>180] = lon[lon>180] - 360
sort_lon = np.argsort(lon)
lon = lon[sort_lon]

# define grid boundaries
lonmin = np.min(lon)
lonmax = np.max(lon)
latmin = np.min(lat)
latmax = np.max(lat)
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

# Shift central longitude so the Southern Hemisphere and North Indin Ocean domains are plotted continuously
if np.logical_and(lonmax >= 90, lonmax <=180):
    central_longitude = 90
else:
    central_longitude = -90
print('central longitude: ',central_longitude)

# sort var according to the new longitude
var = var[:,sort_lon]

#================================================================
var_name= 'sst'
units = '($^oC$)'

# create figure and axes instances
fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
ax.axis('scaled')

cflevels = np.linspace(16, 32, 17)
#cmap = plt.get_cmap('turbo')
cmap = plt.get_cmap('jet')
cf = ax.contourf(lon, lat, var, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
lb = ax.contour(lon, lat, var, levels=[26], colors='grey', alpha=0.7,transform=ccrs.PlateCarree(),linewidths=0.5)
ax.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=20, shrink=0.6, extendrect=True, ticks=cflevels[::4])
cb.ax.tick_params(labelsize=8)

if conf['trackon']=='yes':
    lon_adeck[np.logical_or(lon_adeck<lonmin,lon_adeck>lonmax)] = np.nan
    ax.plot(lon_adeck,lat_adeck,'-ok',markersize=2,alpha=0.4,transform=ccrs.PlateCarree(central_longitude=0))
    nhour = int((int(conf['fhhh'][1:])/3))
    if nhour <= len(fhour):
        ax.plot(lon_adeck[nhour],lat_adeck[nhour],'ok',markersize=6,alpha=0.4,markerfacecolor='None',transform=ccrs.PlateCarree(central_longitude=0))
   
ax.set_extent([lonmin_raw, lonmax_raw, latmin, latmax], crs=ccrs.PlateCarree())

# Add gridlines and labels
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 20))
gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 10))
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

# Add borders and coastlines
#ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='whitesmoke')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
title_center = 'Sea Surface Temperature (${^o}$C)'
ax.set_title(title_center, loc='center', y=1.05, fontsize=8)
title_left = conf['stormModel']+' '+conf['stormName']+conf['stormID']
ax.set_title(title_left, loc='left', fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title_right, loc='right', fontsize=8)

pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
plt.savefig(pngFile,bbox_inches='tight',dpi=150)
#plt.close("all")

