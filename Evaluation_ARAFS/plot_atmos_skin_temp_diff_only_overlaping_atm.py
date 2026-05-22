#!/usr/bin/env python3

"""This script is to plot out HAFS atmospheric surface temperature, MSLP and 10-m wind."""

import os

import yaml
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import grib2io
import xarray as xr

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.interpolate import griddata

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

xlim = conf['xlim']
ylim = conf['ylim']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

####################################################################
# Need to find non overlapping areas for atm.
COMmodel_ocean = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/'
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(COMmodel_ocean+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

# Read ocean lat and lon and land mask
ocean_name = conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'
ocean_file = os.path.join(COMmodel_ocean+conf['ymdh']+'/00E/', ocean_name)
ocean = xr.open_dataset(ocean_file)
lon_ocean = np.asarray(ocean['geolon'])
lat_ocean = np.asarray(ocean['geolat'])
wet_ocean = np.asarray(ocean['wet_c'])[1:,1:]

# interpolate lon_ocean onto lon atmosphere
# --- Original grid ---
lon_orig = lon_ocean + 180
lat_orig = lat_ocean

# Example data on original grid
#data = np.ones((lon_ocean.shape[0],lon_ocean.shape[1]))
wet_ocean[wet_ocean==0] = np.nan
data = wet_ocean

# --- Target grid ---
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
lon_atm = np.copy(lon)
lat_atm = np.copy(lat)
if abs(np.max(lon_atm) - 360.) < 10.:
    lon_atm[lon_atm>180] = lon_atm[lon_atm>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon_atm = lon_atm - lon_offset

lon_targ = lon_atm
lat_targ = lat_atm

# --- Interpolate ---
points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
values = data.ravel()
data_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')

#================================================================
# Read FV3 files
fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+conf['fhhh']+'.grb2'
grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdh']+'/00E/', fname)

print(f'grib2file: {grib2file}')
grb = grib2io.open(grib2file,mode='r')

print('Extracting lat, lon')
lat = np.asarray(grb.select(shortName='NLAT')[0].data)
lon = np.asarray(grb.select(shortName='ELON')[0].data)

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

print('Extracting skin temperature')
levstr='surface'

grib2file = os.path.join(conf['COMmodels'][0]+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')
tmp1 = grb.select(shortName='TMP',level=levstr)[0].data
tmp1 = tmp1 - 273.15 # convert K to degC

grib2file = os.path.join(conf['COMmodels'][1]+conf['ymdh']+'/00E/', fname)
grb = grib2io.open(grib2file,mode='r')
tmp2 = grb.select(shortName='TMP',level=levstr)[0].data
tmp2 = tmp2 - 273.15 # convert K to degC

diff = (tmp2 - tmp1) * data_interp

#===================================================================================================
print('Plotting skin temp. diff')

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure()
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cflevels = np.arange(-2,2.1,0.2)
cf = ax.contourf(lon, lat, diff, levels=cflevels, cmap='coolwarm', extend='both', transform=transform)
cb = plt.colorbar(cf, orientation='vertical', pad=0.02, shrink=0.85, extendrect=True, ticks=cflevels)

plt.xlim([xlim[0]+lon_offset, xlim[1]+lon_offset])
plt.ylim([ylim[0], ylim[1]])

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

title = 'Skin Temperature Difference (${^o}$C) \n' + conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')

ax.set_title(title, loc='center')

ax.text(0.05,0.1,'Mean = ' + str(np.round(np.nanmean(diff),1)) + ' ($^oC$)' + ' Max = ' + str(np.round(np.nanmax(diff),1)) + ' ($^oC$)' + ' Min = ' + str(np.round(np.nanmin(diff),1)) + ' ($^oC$)',transform=ax.transAxes)

'''
model_info = os.environ.get('TITLEgraph','').strip()
var_info = '2 m Temperature Difference(${^o}$C, shaded)'
storm_info = conf['stormID']
title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
'''

#plt.savefig(fig_name, bbox_inches='tight')
#plt.close(fig)
