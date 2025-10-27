#!/usr/bin/env python3

"""This script is to plot out ARAFS atmospheric surface temperature and ATAFS ocean SST"""

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import grib2io
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

################################################################
# Surface temp from atmosphere
fname = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMarafs'], fname)
print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

print('Extracting surface temperature')
levstr='surface'
tmp = grb.select(shortName='TMP', level=levstr)[0].data
tmp = tmp - 273.15 # convert K to degC

xlim = conf['xlim']
ylim = conf['ylim']

# Transform to geographic coordinated
lon[lon>180] = lon[lon>180] - 360

# subset
oklon = np.logical_and(lon[0,:]>=xlim[0],lon[0,:]<xlim[1]) 
oklat = np.logical_and(lat[:,0]>=ylim[0],lat[:,0]<ylim[1]) 

lon_sub = lon[oklat,:][:,oklon]
lat_sub = lat[oklat,:][:,oklon]
tmp_sub = tmp[oklat,:][:,oklon]

################################################################
# SST from ocean
fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'
ncfile = os.path.join(conf['COMarafs'], fname_ocn)
nc = xr.open_dataset(ncfile)

sst = np.asarray(nc['SST'][0,:,:])
lono = np.asarray(nc.xh)
lato = np.asarray(nc.yh)

# subset
oklono = np.logical_and(lono>=xlim[0],lono<xlim[1])
oklato = np.logical_and(lato>=ylim[0],lato<ylim[1])

lono_sub = lono[oklono]
lato_sub = lato[oklato]
sst_sub = sst[oklato,:][:,oklono]

#################################################################
# interpolation from fv3 to MOM6
# --- Original grid ---
lon_orig = lon_sub
lat_orig = lat_sub
# Example data on original grid
data = tmp_sub

# --- Target grid ---
lon_targ, lat_targ =  np.meshgrid(lono_sub, lato_sub)

# --- Interpolate ---
points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
values = data.ravel()
tmp_sub_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')
#tmp_sub_interp = griddata(points, values, (lon_targ, lat_targ), method='nearest')

#################################################################
print('Plotting surface temperature atmosphere')

plt.figure(figsize=(12,8))
cflevels = np.arange(8,29,2)
cmap = plt.get_cmap('turbo')
cf = plt.contourf(lon_sub,lat_sub,tmp_sub, levels=cflevels, cmap=cmap, extend='both')
plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
lb = plt.contour(lon_sub,lat_sub,tmp_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
plt.axis('scaled')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Surface Temperature (${^o}$C, shaded)'
title_left = """{0}{1}""".format(model_info,var_info)
plt.title(title_left, loc='left', y=0.99)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
plt.title(title_right, loc='right', y=0.99)

#################################################################
print('Plotting SST')

plt.figure(figsize=(12,8))
cflevels = np.arange(8,29,2)
cmap = plt.get_cmap('turbo')
cf = plt.contourf(lono_sub,lato_sub,sst_sub, levels=cflevels, cmap=cmap, extend='both')
plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
lb = plt.contour(lono_sub,lato_sub,sst_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
plt.axis('scaled')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'SST (${^o}$C, shaded)'
title_left = """{0}{1}""".format(model_info,var_info)
plt.title(title_left, loc='left', y=0.99)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
plt.title(title_right, loc='right', y=0.99)

#################################################################
print('Plotting surface temperature atmosphere interpolated to MOM6 grid')

plt.figure(figsize=(12,8))
cflevels = np.arange(8,29,2)
cmap = plt.get_cmap('turbo')
cf = plt.contourf(lono_sub,lato_sub,tmp_sub_interp, levels=cflevels, cmap=cmap, extend='both')
plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
lb = plt.contour(lon_sub,lat_sub,tmp_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
plt.axis('scaled')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Surface Temperature Interpolated(${^o}$C, shaded)'
title_left = """{0}{1}""".format(model_info,var_info)
plt.title(title_left, loc='left', y=0.99)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
plt.title(title_right, loc='right', y=0.99)

#################################################################
print('Plotting surf temp FV3 - SST MOM6')

plt.figure(figsize=(12,8))
cflevels = np.arange(-0.25,0.26,0.01)
cmap = plt.get_cmap('bwr')
cf = plt.contourf(lono_sub,lato_sub,tmp_sub_interp-sst_sub, levels=cflevels, cmap=cmap, extend='both')
plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
plt.axis('scaled')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Surf. temp FV3 - MOM6 (${^o}$C, shaded)'
title_left = """{0}{1}""".format(model_info,var_info)
plt.title(title_left, loc='left', y=0.99)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
plt.title(title_right, loc='right', y=0.99)

#################################################################
print('Surf. temp FV3 - MOM6 bar plot')
diff = np.ravel(tmp_sub_interp-sst_sub)

plt.figure()
plt.hist(diff,bins=1000)
plt.xlim([-0.35,0.35])
