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

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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

cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

#################################################################
# Obtain time series for TMP and WTMP
ff = np.arange(0,121,3)
xpos = [1000,1000,1000,1000,1000,1000]
ypos = [500,1000,1500,2000,2200,2300]

tmp = np.empty((len(conf['COMmodels']),len(ff),len(xpos)))
tmp[:] = np.nan
wtmp = np.empty((len(conf['COMmodels']),len(ff),len(xpos)))
wtmp[:] = np.nan

for ex,exp in enumerate(conf['COMmodels']):
    print(exp)

    for f,fh in enumerate(ff):
        if len(str(fh))==3:
            fhour = str(fh)
        if len(str(fh))==2:
            fhour = '0'+str(fh)
        if len(str(fh))==1:
            fhour = '00'+str(fh)
        print(fhour)

        fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.f'+fhour+'.grb2'
        grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/00E/',fname)
        grb = grib2io.open(grib2file,mode='r')
    
        print('Extracting lat, lon')
        lat = grb.select(shortName='NLAT')[0].data
        lon = grb.select(shortName='ELON')[0].data
    
        for i in np.arange(len(xpos)):
            print('Extracting surface temperature')
            levstr='surface'
            tmp[ex,f,i] = grb.select(shortName='TMP', level=levstr)[0].data[ypos[i],xpos[i]] - 273.15 
        
            print('Extracting water temperature')
            wtmp[ex,f,i] = grb.select(shortName='WTMP', level=levstr)[0].data[ypos[i],xpos[i]] - 273.15

color = ['red','green','orange','lime']
label = ['uncoupled','coupled']
for i in np.arange(len(xpos)):
    plt.figure(figsize=(10,6))
    plt.plot(ff,tmp[0,:,i],'.-',color=color[0],label=label[0]+' TMP')
    plt.plot(ff,tmp[1,:,i],'.-',color=color[1],label=label[1]+ ' TMP')
    plt.plot(ff,wtmp[0,:,i],'.-',color=color[2],label=label[0]+' WTMP')
    plt.plot(ff,wtmp[1,:,i],'.-',color=color[3],label=label[1]+' WTMP')
    plt.title('TMP Lat='+str(lat[ypos[i],xpos[i]])+' Lon='+str(lon[ypos[i],xpos[i]]-360))
    plt.legend()

################################################################
# Read Surface temp from experiment
for ex,exp in enumerate(conf['COMmodels']):
    print(exp)
    fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
    grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/00E/',fname)
    grb = grib2io.open(grib2file,mode='r')
    
    print('Extracting lat, lon')
    lat = grb.select(shortName='NLAT')[0].data
    lon = grb.select(shortName='ELON')[0].data
    
    print('Extracting surface temperature')
    levstr='surface'
    tmp = grb.select(shortName='TMP', level=levstr)[0].data
    tmp = tmp - 273.15 # convert K to degC

    print('Extracting water temperature')
    wtmp = grb.select(shortName='WTMP', level=levstr)[0].data
    wtmp = wtmp - 273.15 # convert K to degC
     
    if conf['exp_labels'][ex]=='Uncoupled':
        tmp_uncoupled = tmp
        wtmp_uncoupled = wtmp
    if conf['exp_labels'][ex]=='Coupled':
        tmp_coupled = tmp
        wtmp_coupled = wtmp

    # For cartopy plotting
    print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
    if abs(np.max(lon) - 360.) < 10.:
        lon[lon>180] = lon[lon>180] - 360.
        lon_offset = 0.
    else:
        lon_offset = 180.
    lon = lon - lon_offset
    print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
        
    #################################################################
    print('Plotting TMP surface')
    
    myproj = ccrs.PlateCarree(lon_offset)
    transform = ccrs.PlateCarree(lon_offset)
    
    fig = plt.figure()
    ax = plt.axes(projection=myproj)
    ax.axis('scaled')
    
    cflevels = np.arange(-3,31,2)
    cf = plt.contourf(lon,lat,tmp,levels=cflevels,cmap='turbo',extend='both',alpha=1,transform=transform)
    #cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
    cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True)
    
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    
    plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
    plt.ylim([ylim[0], ylim[1]])
    
    title = 'TMP Surface (${^o}$C) ' + conf['exp_labels'][ex] +'\n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    
    ax.set_title(title, loc='center')

    #################################################################
    print('Plotting WTMP')
    
    myproj = ccrs.PlateCarree(lon_offset)
    transform = ccrs.PlateCarree(lon_offset)
    
    fig = plt.figure()
    ax = plt.axes(projection=myproj)
    ax.axis('scaled')
    
    cflevels = np.arange(-3,31,2)
    cf = plt.contourf(lon,lat,wtmp,levels=cflevels,cmap='turbo',extend='both',alpha=1,transform=transform)
    #cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
    cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True)
    
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    
    plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
    plt.ylim([ylim[0], ylim[1]])
    
    title = 'WTMP (${^o}$C) ' + conf['exp_labels'][ex] +'\n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    ax.set_title(title, loc='center')

    #################################################################
    print('Plotting TMP - WTMP')
    
    myproj = ccrs.PlateCarree(lon_offset)
    transform = ccrs.PlateCarree(lon_offset)
    
    fig = plt.figure()
    ax = plt.axes(projection=myproj)
    ax.axis('scaled')
    
    cflevels = np.arange(-2,2.1,0.1)
    cf = plt.contourf(lon,lat,tmp-wtmp,levels=cflevels,cmap='bwr',extend='both',alpha=1,transform=transform)
    #cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
    cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True)
    
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    
    plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
    plt.ylim([ylim[0], ylim[1]])
    
    title = 'TMP - WTMP (${^o}$C) ' + conf['exp_labels'][ex] +'\n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    ax.set_title(title, loc='center')

#################################################################
print('Plotting TMP surface difference')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-2,2.1,0.1)
cf = plt.contourf(lon,lat,tmp_coupled-tmp_uncoupled,levels=cflevels,cmap='bwr',extend='both',alpha=1,transform=transform)
#cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True)

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

title = 'TMP Surface (${^o}$C) ' + '(Coupled - Uncoupled) \n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

#################################################################
print('Plotting WTMP difference')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-2,2.1,0.1)
cf = plt.contourf(lon,lat,wtmp_coupled-wtmp_uncoupled,levels=cflevels,cmap='bwr',extend='both',alpha=1,transform=transform)
#cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=50, shrink=cbshrink, extendrect=True, ticks=cflevels[::2])
cb = plt.colorbar(cf, orientation='vertical', pad=0.03, shrink=0.85, extendrect=True)

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

title = 'WTMP Surface (${^o}$C) ' + '(Coupled - Uncoupled) \n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title, loc='center')
