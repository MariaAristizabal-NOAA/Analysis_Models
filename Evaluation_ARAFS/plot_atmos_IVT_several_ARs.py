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
print('Parse the config file: plot_atmos_IVT_several_ARs.yml:')
with open('plot_atmos_IVT_several_ARs.yml', 'rt') as f:
    conf = yaml.safe_load(f)

ymdhs = conf['ymdhs']
fhhhs = conf['fhhhs']

xlim = conf['xlim']
ylim = conf['ylim']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']

#================================================================
# Read lat and lot from FV3 file
fname = conf['stormModel'].lower()+'.'+ymdhs[0]+'.'+fhhhs[0]+'.grb2'
grib2file = os.path.join(conf['COMarafs']+ymdhs[0]+'/00E/', fname)
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
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

bounds = [250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600]
#bounds = [600,800, 1000, 1200, 1400, 1600]

fig = plt.figure(figsize=(10,7))
ax = plt.axes(projection=myproj)
ax.axis('scaled')

ax.add_feature(cfeature.LAND, facecolor='black')
ax.add_feature(cfeature.OCEAN, facecolor='midnightblue')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='dimgrey')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='white', edgecolor='dimgrey')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

for f,ymdh in enumerate(ymdhs):
    fname = conf['stormModel'].lower()+'.'+ymdh+'.'+fhhhs[f]+'.grb2'
    grib2file = os.path.join(conf['COMarafs']+ymdh+'/00E/', fname)
    grb = grib2io.open(grib2file,mode='r')
    lon = grb.select(shortName='ELON')[0].data

    IVT_u, IVT_v, IVT = compute_ivt(grb,nlat, nlon)

    IVT[IVT<500] = np.nan

    cf = plt.contourf(lon_atm,lat_atm,IVT,levels=bounds,colors=colors_ivt(),extend='both',alpha=0.4,transform=transform)
    #plt.contour(lon_atm,lat_atm,IVT,levels=[510],colors='grey',linewidths=1,extend='both',alpha=1,transform=transform)

for f,ymdh in enumerate(ymdhs):
    fname = conf['stormModel'].lower()+'.'+ymdh+'.'+fhhhs[f]+'.grb2'
    grib2file = os.path.join(conf['COMarafs']+ymdh+'/00E/', fname)
    grb = grib2io.open(grib2file,mode='r')
    lon = grb.select(shortName='ELON')[0].data

    IVT_u, IVT_v, IVT = compute_ivt(grb,nlat, nlon)

    IVTt = np.copy(IVT)
    IVTt[lon_atm < -150+lon_offset] = np.nan
    IVTt[lon_atm > xlim[1]+lon_offset] = np.nan

    latmax = np.where(IVTt == np.nanmax(IVTt))[0][0]
    lonmax = np.where(IVTt == np.nanmax(IVTt))[1][0]
    lat_max = lat_atm[latmax,lonmax]
    lon_max = lon_atm[latmax,lonmax]

    ax.text(lon_max,lat_max,str(f+1),color='limegreen',fontweight='bold',fontsize=14)

    latt = [lat_atm[latmax,lonmax]]
    lonn = [lon_atm[latmax,lonmax]]
    ivtu = [IVT_u[latmax,lonmax]]
    ivtv = [IVT_v[latmax,lonmax]]
    Q = ax.quiver(lonn, latt, ivtu, ivtv,scale=150, scale_units='xy', width=0.003, headlength=6, headwidth=4,alpha=1.0,transform=transform)

#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
