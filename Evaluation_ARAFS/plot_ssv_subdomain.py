#!/usr/bin/env python3

"""This scrip plots the sea surface temperature for an area 500 km around the storm eye. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
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

#================================================================
# Read ocean files

fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc' 

ncfile = os.path.join(conf['COMarafs'], fname) 
nc = xr.open_dataset(ncfile)

ssu = np.asarray(nc['SSU'][0,:,:])
ssv = np.asarray(nc['SSV'][0,:,:])
lon = np.asarray(nc.xh)
lat = np.asarray(nc.yh)
lonq = np.asarray(nc.xq) 
latq = np.asarray(nc.yq)

#================================================================
# Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
lon[lon>180] = lon[lon>180] - 360
lon[lon<-180] = lon[lon<-180] + 360
sort_lon = np.argsort(lon)
lon = lon[sort_lon]

lonq[lonq>180] = lonq[lonq>180] - 360
lonq[lonq<-180] = lonq[lonq<-180] + 360
sort_lonq = np.argsort(lonq)
lonq = lonq[sort_lonq]

ssu = ssu[:,sort_lonq]
ssv = ssv[:,sort_lon]

skip=6
lnh = lon[::skip]
lth = lat[::skip]
lnq = lonq[::skip]
ltq = latq[::skip]
ssu = ssu[::skip,::skip]
ssv = ssv[::skip,::skip]

lnnq,ltth = np.meshgrid(lnq,lth)
lnnh,lttq = np.meshgrid(lnh,ltq)
lnnh,ltth = np.meshgrid(lnh,lth)

#interpolator = RegularGridInterpolator((lth,lnq),ssu,bounds_error=False, fill_value=None)
#ssu_interp = interpolator((lttq,lnnh))

interpolator = RegularGridInterpolator((lth,lnq),ssu,bounds_error=False, fill_value=None)
ssu_interp = interpolator((ltth,lnnh))

interpolator = RegularGridInterpolator((ltq,lnh),ssv,bounds_error=False, fill_value=None)
ssv_interp = interpolator((ltth,lnnh))

ssV_interp = np.sqrt(ssu_interp**2 + ssv_interp**2)

#================================================================
# create figure and axes instances
var_name= 'ssu'
units = '($m/s$)'
cflevels = np.arange(-1,1.1,0.1)
cmap = plt.get_cmap('turbo')

fig = plt.figure(figsize=(6,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.axis('scaled')
cf = ax.contourf(lnh, lth, ssV_interp, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
ax.contour(lnh, lth, ssV_interp, cflevels, colors='grey',alpha=0.5, linewidths=0.5, transform=ccrs.PlateCarree())
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=20, shrink=0.6, extendrect=True, ticks=cflevels[::4])
cb.ax.tick_params(labelsize=8)

q = plt.quiver(lnnh,ltth,ssu_interp,ssv_interp,scale=20,transform=ccrs.PlateCarree())

plt.xlim(xlim)
plt.ylim(ylim)

# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Sea Surface Velocity ($^oC$)'
storm_info = conf['stormID']
title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
#footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
#ax.text(1.0,-0.08, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)

#pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
#plt.close("all")
        

