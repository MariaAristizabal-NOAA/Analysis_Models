#!/usr/bin/env python3

"""This scrip plots the sea surface height for the entire ocean domain. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
import grib2io
import datetime as datetime

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.ndimage import gaussian_filter

#================================================================
def colors_iwv():
    colors=[(255,255,255,255),(0, 0, 204, 255), (0, 31, 255, 255), (0, 112, 255, 255), (0, 194, 255, 255), (1, 249, 236, 255),\
            (3, 225, 159, 255), (6, 200, 82, 255), (18, 182, 14, 255), (97, 206, 16, 255), (176, 231, 5, 255), (255, 255,0,255),\
            (225, 225, 0,255), (255,195, 0, 255), (255, 165, 0, 255), (255, 112, 0, 255), (255, 59, 0, 255), (255, 7, 0, 255),\
            (208, 0, 38, 255), (154, 0, 82, 255), (100, 0, 125, 255), (127, 61, 165, 255)]
    colors = [tuple(ti/255.0 for ti in element) for element in colors]
    return colors

#================================================================
def compute_iwv(grb, nlat, nlon, g=9.80665):
    """
    IWV (a.k.a. precipitable water) = (1/g) ∫ q dp
    Sign-safe: integrates with ascending pressure so dp > 0.
    Returns IWV
    """

    Ps_pa = np.array([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925, 950, 975, 1000])

    Qs = np.empty((len(Ps_pa),nlat,nlon))
    Qs[:] = np.nan

    for lv,Ps in enumerate(Ps_pa):
        print(lv,' ',Ps)
        Qs[lv,:,:] = grb.select(fullName="Specific Humidity",level=str(Ps)+' mb')[0].data * 100

    # trapezoidal layer means
    dp   = np.diff(Ps_pa)                                    # (L-1,), all > 0
    qbar = 0.5 * (Qs[:-1] + Qs[1:])
    dp3  = dp[:, None, None]

    IWV_kgm2 = np.sum(qbar * dp3, axis=0) / g           # kg m^-2
    IWV_mm   = IWV_kgm2                                     # 1 kg m^-2 == 1 mm

    # Clean tiny negative numerical noise
    IWV_mm = np.where(IWV_mm < 0, np.clip(IWV_mm, 0, None), IWV_mm)
    IWV = IWV_mm

    return IWV

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

cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Read FV3 files
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'

grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdh']+'/00E/', fname)
print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat_atm = grb.select(shortName='NLAT')[0].data
lon_atm = grb.select(shortName='ELON')[0].data

nlat = lat_atm.shape[0]
nlon = lat_atm.shape[1]

IWV1 = compute_iwv(grb,nlat, nlon)

grib2file = os.path.join(conf['COMmodels'][1]+conf['ymdh']+'/00E/', fname)
print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

IWV2 = compute_iwv(grb,nlat, nlon)

diff = gaussian_filter(IWV2 - IWV1,10)

#================================================================
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
#================================================================
var_name= 'Diff IWV'
units = '($Kg m^{-1} s^{-1}$)'

bounds = np.arange(-5,5.1,1)
ticks = np.arange(-5,5.1,1)

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

#fig = plt.figure(figsize=(8,4))
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')
cf = ax.contourf(lon_atm, lat_atm, diff, levels=bounds, cmap='coolwarm', extend='both', transform=transform)

cb = plt.colorbar(cf, orientation='vertical', pad=0.02, extendrect=True,ticks=ticks,shrink=0.8)

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

title = 'Diff IWV ($mm$) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(diff),2)) + ' (mm) \n' + 'Max = ' + str(np.round(np.nanmax(diff),1)) + ' (mm)' + ' Min = ' + str(np.round(np.nanmin(diff),1)) + ' (mm)',transform=ax.transAxes)

#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
