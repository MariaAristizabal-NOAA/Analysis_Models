#!/usr/bin/env python3

"""This scrip plots integrated vapor transport (IVT). """ 

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
def compute_ivt(grb, nlat, nlon, g=9.80665):

    Ps_pa = np.array([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925, 950, 975, 1000])

    Qs = np.empty((len(Ps_pa),nlat,nlon))
    Qs[:] = np.nan
    Us = np.empty((len(Ps_pa),nlat,nlon))
    Us[:] = np.nan
    Vs = np.empty((len(Ps_pa),nlat,nlon))
    Vs[:] = np.nan

    for lv,Ps in enumerate(Ps_pa):
        print(lv,' ',Ps)
        Qs[lv,:,:] = grb.select(fullName="Specific Humidity",level=str(Ps)+' mb')[0].data * 100
        Us[lv,:,:] = grb.select(fullName="U-Component of Wind",level=str(Ps)+' mb')[0].data
        Vs[lv,:,:] = grb.select(fullName="V-Component of Wind",level=str(Ps)+' mb')[0].data

    # trapezoidal layer means
    dp   = np.diff(Ps_pa)                                    # (L-1,), all > 0
    qbar = 0.5 * (Qs[:-1] + Qs[1:])
    ubar = 0.5 * (Us[:-1] + Us[1:])
    vbar = 0.5 * (Vs[:-1] + Vs[1:])
    dp3  = dp[:, None, None]

    IVT_u = np.sum(qbar * ubar * dp3, axis=0) / g
    IVT_v = np.sum(qbar * vbar * dp3, axis=0) / g
    IVT  = np.hypot(IVT_u, IVT_v)

    return IVT_u, IVT_v, IVT

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
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
IVT_u, IVT_v, IVT = compute_ivt(grb,nlat, nlon)

#================================================================
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

bounds = [250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600]

# create figure and axes instances
#fig,ax = plt.subplots(figsize=(7,4))
#ax.set_facecolor("bisque")

fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-3,34,2)
cmap = 'turbo'

if os.path.exists(ncfile):
    plt.contourf(lon_ocn+180, lat_ocn, var, levels=cflevels, cmap=cmap, extend='both',alpha=0.1,transform=transform)

cf = plt.contourf(lon_atm,lat_atm,IVT,levels=bounds,colors=colors_ivt(),extend='both',alpha=0.8,transform=transform)
plt.contour(lon_atm,lat_atm,IVT,levels=bounds,colors='grey',linewidths=0.2,extend='both',alpha=0.8,transform=transform)
#cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, extendrect=True,shrink=0.8)

skip = 100
#plt.quiver(lon_atm[::skip,::skip]+180, lat_atm[::skip,::skip], IVT_u[::skip,::skip], IVT_v[::skip,::skip],scale=122, scale_units='xy', width=0.003, headlength=6, headwidth=4,alpha=0.5,transform=transform)
Q = ax.quiver(lon_atm[::skip,::skip], lat_atm[::skip,::skip], IVT_u[::skip,::skip], IVT_v[::skip,::skip],scale=122, scale_units='xy', width=0.003, headlength=6, headwidth=4,alpha=1.0,transform=transform)
ax.quiverkey(Q, 1.05, 1.05, 500, r'500 $kg m^{-1} s^{-1}$', labelpos='E', coordinates='axes')

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

title = 'IVT ($Kg m^{-1} s^{-1}$) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title, loc='center')

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(IVT),1)) + ' ($Kg m^{-1} s^{-1}$) ' + 'Max = ' + str(np.round(np.nanmax(IVT),1)) + ' ($Kg m^{-1} s^{-1}$)',transform=ax.transAxes)

#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
