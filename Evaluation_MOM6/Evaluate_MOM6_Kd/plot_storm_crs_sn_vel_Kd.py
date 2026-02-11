#!/usr/bin/env python3

"""This scrip profile if Eddy Diffusivity Kd. """

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import matplotlib
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
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
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
    adeck_file = os.path.join(conf['COMhafs'],adeck_name)

    fhour,lat_adeck,lon_adeck,init_time,valid_time = get_adeck_track(adeck_file)

    print('lon_adeck = ',lon_adeck)
    print('lat_adeck = ',lat_adeck)

#================================================================
# Read ocean files

oceanf = glob.glob(os.path.join(conf['COMhafs'],'*f006.nc'))[0].split('/')[-1].split('.')

ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

if ocean == 'mom6':
    fname000 =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+'f000.nc'
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc'

if ocean == 'hycom':
    fname000 =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+'f000.nc'
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+conf['fhhh']+'.nc'

ncfile000 = os.path.join(conf['COMhafs'], fname000)
nc000 = xr.open_dataset(ncfile000)
ncfile = os.path.join(conf['COMhafs'], fname)
nc = xr.open_dataset(ncfile)

if ocean == 'mom6':
    ssuu = np.asarray(nc['uo'][0,:,:,:])
    ssvv = np.asarray(nc['vo'][0,:,:,:])
    temp = np.asarray(nc['temp'][0,:,:,:])
    kdd_ePBL = np.asarray(nc['Kd_ePBL'][0,:,:,:])
    kdd_shear = np.asarray(nc['Kd_shear'][0,:,:,:])
    lon = np.asarray(nc.xh)
    lat = np.asarray(nc.yh)
    lonq = np.asarray(nc.xq)
    latq = np.asarray(nc.yq)
    mldd = np.asarray(nc['MLD_0125'][0,:,:])
    zl = np.asarray(nc['z_l'])
    zi = np.asarray(nc['z_i'])

if ocean == 'hycom':
    ssu = np.asarray(nc['u_velocity'][0,:,:,:])/100
    ssv = np.asarray(nc['v_velocity'][0,:,:,:])/100
    varr = np.sqrt(ssu**2 + ssv**2)
    mldd = np.asarray(nc['mixed_layer_thickness'][0,:,:])
    zl = np.asarray(nc['Z'])
    lonh = np.asarray(nc.Longitude)
    lath = np.asarray(nc.Latitude)
    lon, lat = HYCOM_coord_to_geo_coord(lonh,lath)

lonmin_raw = np.min(lon)
lonmax_raw = np.max(lon)
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

#================================================================
#%% Latitudinal transect
okfhour = conf['fhhh'][1:] == fhour
if len(lon_adeck[okfhour])==0 or len(lat_adeck[okfhour])==0:
    print('There is not latitude or longitude for the center of the storm at this forecast hour. Exiting plot_storm_crs_sn_temp.py')
    sys.exit()
else:
    lon_eye = lon_adeck[okfhour][0]
    xpos = int(np.round(np.interp(lon_eye,lon,np.arange(len(lon)))))
    lat_eye = lat_adeck[okfhour][0]
    ypos = int(np.round(np.interp(lat_eye,lat,np.arange(len(lat)))))
    
    ssu = ssuu[:,:,xpos]
    ssv = ssvv[:,:,xpos]
    mld = mldd[:,xpos]
    kd_ePBL = kdd_ePBL[:,:,xpos]
    kd_shear = kdd_shear[:,:,xpos]

    #zzl,ltth = np.meshgrid(zl,lonq)
    zll,lttq = np.meshgrid(zl,latq)

    interpolator = RegularGridInterpolator((zl,lat),ssu,bounds_error=False, fill_value=None)
    ssu_interp = interpolator((zll,lttq)).T

    var = np.sqrt(ssu_interp**2 + ssv**2)
    
    #================================================================
    # Temp
    lat_eye = lat_adeck[okfhour][0]

    fig,ax = plt.subplots(figsize=(8,4))
    cmap = plt.colormaps['YlOrRd']
    ncolors = cmap.N
    levels = np.arange(0,2.1,0.1)
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N,extend='max')
    ctr = ax.pcolor(latq,-zl,var,cmap='YlOrRd',norm=norm)
    cbar = fig.colorbar(ctr,extendrect=True)
    #ctr = ax.contourf(latq,-zl,var,cmap='YlOrRd',**kw,extend='both')
    #cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(latq,-zl,var,[26],colors='k')
    ax.plot(lat,-mld,'-',color='green')
    ax.plot(np.tile(lat_eye,len(zl)),-zl,'-k')
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.set_xlim([lat_eye-4,lat_eye+4])
    ax.set_ylim([-150,0])
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Latitude')
    
    if lon_eye >= 0:
        hemis = 'E'
    else:
        hemis = 'W'
        
    model_info = os.environ.get('TITLEgraph','').strip()
    var_info = 'Speed ($m/s$) X-section at  ' + str(np.round(lon_eye,2)) + ' ' + hemis
    storm_info = conf['stormName']+conf['stormID']
    title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
    ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
    title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
    footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
    ax.text(1.0,-0.2, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)
    
    pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.crs_sn_temp'+'.'+conf['fhhh'].lower()+'.png'
    #plt.savefig(pngFile,bbox_inches='tight',dpi=150)
    #plt.close()
    
    ##################
    fig,ax = plt.subplots(figsize=(8,4))
    #ctr = ax.contourf(lat,-zi,kd_ePBL,cmap='plasma',**kw,extend='both')
    cmap = plt.colormaps['plasma']
    ncolors = cmap.N
    levels = np.arange(0,0.21,0.01)
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N,extend='max')
    ctr = ax.pcolor(lat,-zi,kd_ePBL,cmap='plasma',norm=norm)
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$m^2/s$',fontsize=14)
    ax.plot(lat,-mld,'-',color='green')
    ax.plot(np.tile(lat_eye,len(zl)),-zl,'-k')
    ax.set_ylim([-150,0])
    ax.set_xlim([lat_eye-4,lat_eye+4])
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Latitude')

    if lon_eye >= 0:
        hemis = 'E'
    else:
        hemis = 'W'

    model_info = os.environ.get('TITLEgraph','').strip()
    var_info = 'Kd ePBL ($m^2/s$) X-section at  ' + str(np.round(lon_eye,2)) + ' ' + hemis
    storm_info = conf['stormName']+conf['stormID']
    title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
    ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
    title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
    footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
    ax.text(1.0,-0.2, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)

    ##################
    fig,ax = plt.subplots(figsize=(8,4))
    #ctr = ax.contourf(lat,-zi,kd_ePBL,cmap='plasma',**kw,extend='both')
    cmap = plt.colormaps['plasma']
    ncolors = cmap.N
    levels = np.arange(0,0.21,0.01)
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N,extend='max')
    ctr = ax.pcolor(lat,-zi,kd_shear,cmap='plasma',norm=norm)
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$m^2/s$',fontsize=14)
    ax.plot(lat,-mld,'-',color='green')
    ax.plot(np.tile(lat_eye,len(zl)),-zl,'-k')
    ax.set_ylim([-150,0])
    ax.set_xlim([lat_eye-4,lat_eye+4])
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Latitude')

    if lon_eye >= 0:
        hemis = 'E'
    else:
        hemis = 'W'

    model_info = os.environ.get('TITLEgraph','').strip()
    var_info = 'Kd shear ($m^2/s$) X-section at  ' + str(np.round(lon_eye,2)) + ' ' + hemis
    storm_info = conf['stormName']+conf['stormID']
    title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
    ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
    title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
    ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
    footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
    ax.text(1.0,-0.2, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)

    #================================================================
