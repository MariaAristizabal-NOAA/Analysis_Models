#!/usr/bin/env python3

"""This script plots out HAFS the 3-hours accumulated precipitation, the mean sea level pressure and the 1000-500 geopotential thickness."""

import os

import yaml
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import grib2io

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
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

fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)

grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
[nlat, nlon] = np.shape(lon)

print('Extracting 24h accumulate precipitation at surface')

if int(conf['fhhh'][1:]) - 24 >= 0:
    fhhs = np.arange(int(conf['fhhh'][1:])-24,int(conf['fhhh'][1:])+1,3)

    fhhhs = []
    for f in fhhs:
        if f<10: fhhhs.append('f00'+str(f)) 
        if f>10 and f<100: fhhhs.append('f0'+str(f)) 
        if f>100: fhhhs.append('f'+str(f)) 
        
    apcp = np.empty((len(fhhhs),lat.shape[0],lat.shape[1]))
    apcp[:] = np.nan
    for f,fhhh in enumerate(fhhhs):
        fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.'+fhhh+'.grb2'
        grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
        print(grib2file)
        grb = grib2io.open(grib2file,mode='r')
        #apcp[f,:,:] = grb.select(shortName='APCP')[0].data*0.0393701  # convert kg/m^2 to in
        apcp[f,:,:] = grb.select(shortName='APCP')[0].data  # convert kg/m^2 to in

apcp_24 = np.sum(apcp,axis=0)
#===================================================================================================
print('APCP')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# Inches
'''
cflevels = [0.01,0.1,0.25,  # green
            0.5,0.75,1,1.25, # blue
            1.5,1.75,2,     # purple
            2.5,3,4,5, # red
            7,10,15,20] # yellow 

cfcolors = ['lawngreen','mediumseagreen', 'green',              # green
            'darkblue','dodgerblue','cornflowerblue','cyan', # blue
            'mediumslateblue','darkorchid','darkmagenta',       # purple
            'darkred','firebrick','orangered','coral',   # red
            'darkgoldenrod','gold','yellow']       # yellow
'''

# mm
cflevels = [0,   #white
            0.1,1,2,  # green
            4,6,8,10, # blue
            12,14,16,     # purple
            18,20,25,30, # red
            35,40,45,50, # yello
            55]    #pink

cfcolors = ['white',   
            'lawngreen','mediumseagreen', 'green',              # green
            'darkblue','dodgerblue','cornflowerblue','cyan', # blue
            'mediumslateblue','darkorchid','darkmagenta',       # purple
            'darkred','firebrick','orangered','coral',   # red
            'darkgoldenrod','gold','yellow',       # yellow
            'lightsalmon']   # pink

cm = matplotlib.colors.ListedColormap(cfcolors)
norm = matplotlib.colors.BoundaryNorm(cflevels, cm.N)
#cm = matplotlib.colors.ListedColormap(colors_precip)
#norm = matplotlib.colors.BoundaryNorm(bounds_precip, cm.N)
cbshrink = 1.0

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

#cf = ax.contourf(lon, lat, apcp, cflevels, cmap=cm, norm=norm, transform=transform)
cf = ax.contourf(lon, lat, apcp_24, cflevels, cmap=cm, norm=norm, transform=transform,extend='both')
#cb = plt.colorbar(cf, orientation='vertical',shrink=0.8,pad=0.03,extendrect=True,ticks=cflevels)
cb = plt.colorbar(cf, orientation='vertical',shrink=0.8,pad=0.03,ticks=cflevels[1:-1],extend='both',extendfrac=(0.0, 0.0))
cb.ax.set_yticklabels(cflevels[1:-1])

# Add borders and coastlines
#ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='whitesmoke')
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

dt1 = pd.to_timedelta(fhhs[0], unit='h')
time1 = conf['initTime'] + dt1
dt2 = pd.to_timedelta(fhhs[-1], unit='h')
time2 = conf['initTime'] + dt2

title = '24hr Accumulated Precipitation (mm) \n ' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+' Valid: '+time1.strftime('%Y%m%d%HZ ')+'to '+time2.strftime('%Y%m%d%HZ ')

ax.set_title(title, loc='center') #, x=1.2,y=1.01)

ax.text(0.1,0.1,'Mean = ' + str(np.round(np.nanmean(apcp_24),1)) + ' (mm) ' + 'Max = ' + str(np.round(np.nanmax(apcp_24),1)) + ' (mm)',transform=ax.transAxes)


#plt.show()
#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
