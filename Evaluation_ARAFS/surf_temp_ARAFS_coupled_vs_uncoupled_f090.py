#!/usr/bin/env python3

"""This script is to compare ARAFS surface temperature with OSTIA SST"""

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import grib2io
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Parse the yaml config file
print('Parse the config file: plot_atmos3.yml:')
with open('plot_atmos3.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
experiments = conf['Experiments']
data_prod = conf['Data_prod']
fhhhs = conf['fhhhs']
commodels = conf['COMmodels']
comobs = conf['COMobs']
fileobs = conf['Fileobs']

xlim = conf['xlim']
ylim = conf['ylim']

################################################################
# Read Surface temp from experiment
fhhh = 'f090'
for ex,exp in enumerate(experiments):
    print(exp)
    if exp == 'ARAFS' or exp == 'ARAFS_OCN':
        conf['fhours'] = int(fhhh[1:])
        conf['fcstTime'] = pd.to_timedelta(conf['fhours'], unit='h')
        conf['validTime'] = conf['initTime'] + conf['fcstTime']
        fname = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+fhhh+'.grb2'
        grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname)
        if os.path.exists(grib2file): 
            print(f'grib2file: {grib2file}')
        else:
            fname = conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+fhhh+'.grb2'
            grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname)
            print(f'grib2file: {grib2file}')
        grb = grib2io.open(grib2file,mode='r')
        
        print('Extracting lat, lon')
        lat = grb.select(shortName='NLAT')[0].data
        lon = grb.select(shortName='ELON')[0].data
        
        print('Extracting surface temperature')
        levstr='surface'
        tmp = grb.select(shortName='TMP', level=levstr)[0].data
        tmp = tmp - 273.15 # convert K to degC
        
        # Transform to geographic coordinated
        lon[lon>180] = lon[lon>180] - 360
        
        # subset
        oklon = np.logical_and(lon[0,:]>=xlim[0],lon[0,:]<xlim[1]) 
        oklat = np.logical_and(lat[:,0]>=ylim[0],lat[:,0]<ylim[1]) 
        
        lon_sub = lon[oklat,:][:,oklon]
        lat_sub = lat[oklat,:][:,oklon]
        tmp_sub = tmp[oklat,:][:,oklon]

        if exp == 'ARAFS':
            lon_arafs = lon_sub
            lat_arafs = lat_sub
            tmp_arafs = tmp_sub
        if exp == 'ARAFS_OCN':
            tmp_arafs_ocn = tmp_sub
    
    ################################################################
    if exp == 'MOM6':
        # SST from ocean
        fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+fhhh+'.nc'
        ncfile = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname_ocn)
        nc = xr.open_dataset(ncfile)
        
        tmp = np.asarray(nc['SST'][0,:,:])
        uo = np.asarray(nc['uo'][0,0,:,:])
        lon = np.asarray(nc.xh)
        lonq = np.asarray(nc.xq)
        lat = np.asarray(nc.yh)
        
        # subset
        oklon = np.logical_and(lon>=xlim[0],lon<xlim[1])
        oklonq = np.logical_and(lonq>=xlim[0],lonq<xlim[1])
        oklat = np.logical_and(lat>=ylim[0],lat<ylim[1])
        
        lon_sub = lon[oklon]
        lonq_sub = lonq[oklonq]
        lat_sub = lat[oklat]
        tmp_sub = tmp[oklat,:][:,oklon]
        uo_sub = uo[oklat,:][:,oklonq]

        plt.figure(figsize=(12,8))
        cflevels = np.arange(-1,1.1,0.1)
        cmap = plt.get_cmap('turbo')
        cf = plt.contourf(lonq_sub,lat_sub,uo_sub, levels=cflevels, cmap=cmap, extend='both')
        plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
        plt.axis('scaled')

        #################################################################
    plt.figure(figsize=(12,8))
    cflevels = np.arange(8,29,2)
    cmap = plt.get_cmap('turbo')
    cf = plt.contourf(lon_sub,lat_sub,tmp_sub, levels=cflevels, cmap=cmap, extend='both')
    plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
    lb = plt.contour(lon_sub,lat_sub,tmp_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
    plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
    plt.axis('scaled')
 
    model_info = exp
    var_info = ' Surface Temperature (${^o}$C, shaded)'
    title_left = """{0}{1}""".format(model_info,var_info)
    plt.title(title_left, loc='left', y=0.99)
    title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    plt.title(title_right, loc='right', y=0.99)
 
    #################################################################
plt.figure(figsize=(12,8))
cflevels = np.arange(-2,2.1,0.2)
cmap = plt.get_cmap('bwr')
cf = plt.contourf(lon_arafs,lat_arafs,tmp_arafs-tmp_arafs_ocn, levels=cflevels, cmap=cmap, extend='both')
plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
plt.axis('scaled')

model_info = 'ARAFS - ARAFS_OCN \n'
var_info = 'Surface Temperature (${^o}$C, shaded)'
title_left = """{0}{1}""".format(model_info,var_info)
plt.title(title_left, loc='left', y=0.99)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
plt.title(title_right, loc='right', y=0.99)



