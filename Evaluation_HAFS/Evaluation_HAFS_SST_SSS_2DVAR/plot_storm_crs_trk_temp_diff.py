#!/usr/bin/env python3

"""This scrip plots an along forecast track transect of temperature. """

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

#================================================================
def latlon_str2num(string):
    """Convert lat/lon string into numbers."""
    value = pd.to_numeric(string[:-1], errors='coerce') / 10
    if string.endswith(('N', 'E')):
        return value
    else:
        return -value

#================================================================
def get_adeck_track(adeck_file):

    cols = ['basin', 'number', 'ymdh', 'technum', 'tech', 'tau', 'lat', 'lon', 'vmax', 'mslp', 'type','rad', 'windcode', 'rad1', 'rad2', 'rad3', 'rad4', 'pouter', 'router', 'rmw', 'gusts', 'eye','subregion', 'maxseas', 'initials', 'dir', 'speed', 'stormname', 'depth','seas', 'seascode', 'seas1', 'seas2', 'seas3', 'seas4', 'userdefined','userdata1', 'userdata2', 'userdata3', 'userdata4', 'userdata5', 'userdata6', 'userdata7', 'userdata8','userdata9', 'userdata10', 'userdata11', 'userdata12', 'userdata13', 'userdata14', 'userdata15', 'userdata16']

    # Read in adeckFile as pandas dataFrame
    print('Read in adeckFile ...')
    adeck = pd.read_csv(adeck_file, index_col=False, names=cols, dtype=str, header=None, skipinitialspace=True)

    adeck['lat'] = adeck['lat'].apply(latlon_str2num)
    adeck['lon'] = adeck['lon'].apply(latlon_str2num)
    #adeck['vmax'] = adeck['vmax'].apply(lambda x: x if x>0 else np.nan)
    #adeck['mslp'] = adeck['mslp'].apply(lambda x: x if x>800 else np.nan)
    adeck['init_time'] = pd.to_datetime(adeck['ymdh'], format='%Y%m%d%H', errors='coerce')
    adeck['valid_time'] = pd.to_timedelta(adeck['tau'].apply(pd.to_numeric, errors='coerce'), unit='h') + adeck['init_time']

    fhour, ind = np.unique(adeck['tau'],return_index=True)
    lat_adeck = np.asarray(adeck['lat'][ind])
    lon_adeck = np.asarray(adeck['lon'][ind])
    init_time = np.asarray(adeck['init_time'][ind])
    valid_time = np.asarray(adeck['valid_time'][ind])

    return fhour,lat_adeck,lon_adeck,init_time,valid_time

#================================================================
# Conversion from geographic longitude and latitude to HYCOM convention
def HYCOM_coord_to_geo_coord(lonh,lath):
    lonh = np.asarray(lonh)
    if np.ndim(lonh) > 0:
        long = [ln-360 if ln>=180 else ln for ln in lonh]
    else:
        long = [lonh-360 if lonh>=180 else lonh][0]
    latg = lath
    return long, latg

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_ocean2.yml:')
with open('plot_ocean2.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

#================================================================
# Get lat and lon from adeck file

if conf['trackon']=='yes':
    adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
    adeck_file1 = os.path.join(conf['COMhafs1'],adeck_name)
    adeck_file2 = os.path.join(conf['COMhafs2'],adeck_name)

    fhour1,lat_adeck1,lon_adeck1,init_time1,valid_time1 = get_adeck_track(adeck_file1)
    #fhour2,lat_adeck2,lon_adeck2,init_time2,valid_time2 = get_adeck_track(adeck_file2)

    print('lon_adeck1 = ',lon_adeck1)
    print('lat_adeck1 = ',lat_adeck1)
    #print('lon_adeck2 = ',lon_adeck2)
    #print('lat_adeck2 = ',lat_adeck2)
#================================================================
# Read ocean files

oceanf = glob.glob(os.path.join(conf['COMhafs1'],'*f006.nc'))[0].split('/')[-1].split('.')

ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

if ocean == 'mom6':
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'

if ocean == 'hycom':
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+conf['fhhh']+'.nc'

ncfile1 = os.path.join(conf['COMhafs1'], fname)
nc1 = xr.open_dataset(ncfile1)
ncfile2 = os.path.join(conf['COMhafs2'], fname)
nc2 = xr.open_dataset(ncfile2)

if ocean == 'mom6':
    varr1 = np.asarray(nc1['temp'][0,:,:,:])
    mld1 = np.asarray(nc1['MLD_0125'][0,:,:])
    zl1 = np.asarray(nc1['z_l'])
    lon1 = np.asarray(nc1.xh)
    lat1 = np.asarray(nc1.yh)
    varr2 = np.asarray(nc2['temp'][0,:,:,:])
    mld2 = np.asarray(nc2['MLD_0125'][0,:,:])
    zl2 = np.asarray(nc2['z_l'])
    lon2 = np.asarray(nc2.xh)
    lat2 = np.asarray(nc2.yh)

if ocean == 'hycom':
    mld = np.asarray(nc['mixed_layer_thickness'][0,:,:])
    zl = np.asarray(nc['Z'])
    lonh = np.asarray(nc.Longitude)
    lath = np.asarray(nc.Latitude)
    lon, lat = HYCOM_coord_to_geo_coord(lonh,lath)

lonmin_raw1 = np.min(lon1)
lonmax_raw1 = np.max(lon1)
print('raw lonlat1 limit: ', np.min(lon1), np.max(lon1), np.min(lat1), np.max(lat1))
lonmin_raw2 = np.min(lon2)
lonmax_raw2 = np.max(lon2)
print('raw lonlat2 limit: ', np.min(lon2), np.max(lon2), np.min(lat2), np.max(lat2))

#================================================================
#%% temp along storm path
# file1
lon_adeck_interp1 = np.interp(lat1,lat_adeck1,lon_adeck1,left=np.nan,right=np.nan)
if len(np.where((np.isfinite(lon_adeck_interp1)))[0]) == 0:
    lon_adeck_interp1 = np.interp(lat1,lat_adeck1,lon_adeck1)

lat_adeck_interp1 = np.copy(lat1)
lat_adeck_interp1[np.isnan(lon_adeck_interp1)] = np.nan

lon_adeck_int1 = lon_adeck_interp1[np.isfinite(lon_adeck_interp1)]
lat_adeck_int1 = lat_adeck_interp1[np.isfinite(lat_adeck_interp1)]
lon_adeck_int1[0:len(lon_adeck_int1)] = lon_adeck_int1
lat_adeck_int1[0:len(lat_adeck_int1)] = lat_adeck_int1

oklon1 = np.round(np.interp(lon_adeck_int1,lon1,np.arange(len(lon1)))).astype(int)
oklat1 = np.round(np.interp(lat_adeck_int1,lat1,np.arange(len(lat1)))).astype(int)

'''
# file2
lon_adeck_interp2 = np.interp(lat2,lat_adeck,lon_adeck,left=np.nan,right=np.nan)
if len(np.where((np.isfinite(lon_adeck_interp2)))[0]) == 0:
    lon_adeck_interp2 = np.interp(lat2,lat_adeck,lon_adeck)

lat_adeck_interp2 = np.copy(lat2)
lat_adeck_interp2[np.isnan(lon_adeck_interp2)] = np.nan

lon_adeck_int2 = lon_adeck_interp2[np.isfinite(lon_adeck_interp2)]
lat_adeck_int2 = lat_adeck_interp2[np.isfinite(lat_adeck_interp2)]
lon_adeck_int2[0:len(lon_adeck_int2)] = lon_adeck_int2
lat_adeck_int2[0:len(lat_adeck_int2)] = lat_adeck_int2

oklon2 = np.round(np.interp(lon_adeck_int2,lon2,np.arange(len(lon2)))).astype(int)
oklat2 = np.round(np.interp(lat_adeck_int2,lat2,np.arange(len(lat2)))).astype(int)
'''

var1 = np.empty((len(zl1),len(lon_adeck_int1)))
var1[:] = np.nan 
mldd1 = np.empty((len(lon_adeck_int1)))
mldd1[:] = np.nan
var2 = np.empty((len(zl1),len(lon_adeck_int1)))
var2[:] = np.nan 
mldd2 = np.empty((len(lon_adeck_int1)))
mldd2[:] = np.nan
for x in np.arange(len(lon_adeck_int1)):
    var1[:,x] = varr1[:,oklat1[x],oklon1[x]]
    mldd1[x] = np.asarray(mld1[oklat1[x],oklon1[x]])
    var2[:,x] = varr2[:,oklat1[x],oklon1[x]]
    mldd2[x] = np.asarray(mld2[oklat1[x],oklon1[x]])

diff = var1 - var2

#================================================================
# Temp
'''
okfhour = conf['fhhh'][1:] == fhour
if len(lat_adeck[okfhour])!=0:
    lat_eye = lat_adeck[okfhour][0]
else:
    lat_eye = np.nan

kw = dict(levels=np.arange(15,31.1,0.5))
fig,ax = plt.subplots(figsize=(8,4))
ctr = ax.contourf(lat[oklat],-zl,var,cmap='Spectral_r',**kw,extend='both')
cbar = fig.colorbar(ctr,extendrect=True)
#cbar.set_label('$^oC$',fontsize=14)
cs = ax.contour(lat[oklat],-zl,var,[26],colors='k')
ax.plot(lat[oklat],-mldd,'-',color='green')
ax.plot(np.tile(lat_eye,len(zl)),-zl,'-k')
ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
ax.set_ylim([-300,0])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Latitude')

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Temperature ($^oC$) X-section along track'
storm_info = conf['stormName']+conf['stormID']
title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
ax.text(1.0,-0.2, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)

pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.crs_trk_temp'+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
#plt.close()
'''
#================================================================
# Temp1 - Temp2
kw = dict(levels=np.arange(-0.3,0.31,0.05))
fig,ax = plt.subplots(figsize=(8,4))
ctr = ax.contourf(lat1[oklat1],-zl1,diff,cmap='seismic',**kw,extend='both')
cbar = fig.colorbar(ctr,extendrect=True)
#cbar.set_label('$^oC$',fontsize=14)
#cs = ax.contour(lat1[oklat1],-zl1,diff,[0],colors='grey',alpha=0.3)
#ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
ax.plot(lat1[oklat1],-mldd1,'-',color='dodgerblue',label='original')
ax.plot(lat1[oklat1],-mldd2,'-',color='orange',label='modified')
ax.set_ylim([-300,0])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Latitude')
ax.legend()

model_info = os.environ.get('TITLEgraph','').strip()
var_info = 'Temperature Difference ($^oC$) X-section along track. Min dT = ' + str(np.round(np.nanmin(diff),2))
storm_info = conf['stormName']+conf['stormID']
title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
ax.text(1.0,-0.2, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)

pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.crs_trk_temp'+'.change.'+conf['fhhh'].lower()+'.png'
plt.savefig(pngFile,bbox_inches='tight',dpi=150)
#plt.close()

fig,ax = plt.subplots(figsize=(8,4))
ax.plot(lat1[oklat1],-mldd1,'o-',color='dodgerblue',label='original')
ax.plot(lat1[oklat1],-mldd2,'o-',color='orange',label='modified')
ax.legend()
ax.set_title('Mixed Layer Depth',fontsize=14)
ax.set_ylabel('(m)',fontsize=14)
ax.set_xlabel('Latitude',fontsize=14)
