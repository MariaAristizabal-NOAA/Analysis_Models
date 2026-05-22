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

#================================================================
def colors_ivt():
    colors=[(255,255,255,255),(255, 255, 0, 255),(255,229, 0, 255), (255,201, 0, 255), (255, 173, 0, 255), \
            (255, 130, 0, 255), (255, 80, 0, 255), (255, 30, 0, 255), (235, 0, 16, 255), (184, 0, 58, 255),\
            (133, 0, 99, 255), (87, 0, 136, 255)]

    colors = [tuple(ti/255.0 for ti in element) for element in colors]
    return colors

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
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
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

#================================================================
# Read ocean files
fname =  conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'

ncfile = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
if os.path.exists(ncfile):
    nc = xr.open_dataset(ncfile)
    var = np.asarray(nc['SST'][0,:,:])
    lon_ocn = np.asarray(nc.geolon)
    lat_ocn = np.asarray(nc.geolat)
else:
    print('file ',ncfile,' does not exist.')

#================================================================
# Read FV3 files
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat_atm = grb.select(shortName='NLAT')[0].data
lon_atm = grb.select(shortName='ELON')[0].data

nlat = lat_atm.shape[0]
nlon = lat_atm.shape[1]

IWV = compute_iwv(grb,nlat, nlon)

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
bounds = [20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
ticks = [20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60]

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-3,34,2)
cmap = 'turbo'
if os.path.exists(ncfile):
    ax.contourf(lon_ocn+180, lat_ocn, var, levels=cflevels, cmap=cmap, extend='both',alpha=0.1,transform=transform)

cf = ax.contourf(lon_atm,lat_atm,IWV,levels=bounds,colors=colors_iwv(),extend='both',alpha=0.8,transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, shrink=0.8,extendrect=True,ticks=ticks)

ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

title = 'IWV ($mm$) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(IWV),1)) + ' (mm) ' + 'Max = ' + str(np.round(np.nanmax(IWV),1)) + ' (mm)',transform=ax.transAxes)


#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
