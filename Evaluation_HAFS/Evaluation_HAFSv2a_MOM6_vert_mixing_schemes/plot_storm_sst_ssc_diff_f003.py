#!/usr/bin/env python3

"""This scrip plots the sea surface temperature for an area 500 km around the storm eye. """ 

import os
import sys
import glob
import yaml

import xarray as xr
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
  
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from geo4HYCOM import haversine

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
# Parse the yaml config file
print('Parse the config file: plot_ocean2.yml:')
with open('plot_ocean2.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
conf['fhour'] = int(conf['fhhh'][1:])
conf['fcstTime'] = pd.to_timedelta(conf['fhour'], unit='h')
conf['validTime'] = conf['initTime'] + conf['fcstTime']

# Set Cartopy data_dir location
cartopy.config['data_dir'] = conf['cartopyDataDir']
print(conf)

#================================================================
# Get lat and lon from adeck file

if conf['trackon']=='yes':
    adeck_name = conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.trak.atcfunix'
    adeck_file1 = os.path.join(conf['COMhafs1'],adeck_name)
    fhour1,lat_adeck1,lon_adeck1,init_time1,valid_time1 = get_adeck_track(adeck_file1)

    adeck_file2 = os.path.join(conf['COMhafs2'],adeck_name)
    fhour2,lat_adeck2,lon_adeck2,init_time2,valid_time2 = get_adeck_track(adeck_file2)

    print('lon_adeck = ',lon_adeck1)
    print('lat_adeck = ',lat_adeck1)

#================================================================
# Read ocean files

oceanf = glob.glob(os.path.join(conf['COMhafs1'],'*f006.nc'))[0].split('/')[-1].split('.')

ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

if ocean == 'mom6':
    fname000 =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+'f003.nc' 
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+conf['fhhh']+'.nc' 

if ocean == 'hycom':
    fname000 =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+'f000.nc'
    fname =  conf['stormID'].lower()+'.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.hycom.3z.'+conf['fhhh']+'.nc'

ncfile0001 = os.path.join(conf['COMhafs1'], fname000)
nc0001 = xr.open_dataset(ncfile0001)
ncfile1 = os.path.join(conf['COMhafs1'], fname) 
nc1 = xr.open_dataset(ncfile1)

ncfile0002 = os.path.join(conf['COMhafs2'], fname000)
nc0002 = xr.open_dataset(ncfile0002)
ncfile2 = os.path.join(conf['COMhafs2'], fname) 
nc2 = xr.open_dataset(ncfile2)

if ocean == 'mom6':
    varr0001 = np.asarray(nc0001['SST'][0,:,:])
    varr1 = np.asarray(nc1['SST'][0,:,:])
    ssu1 = np.asarray(nc1['SSU'][0,:,:])
    ssv1 = np.asarray(nc1['SSV'][0,:,:])
    lon = np.asarray(nc1.xh)
    lat = np.asarray(nc1.yh)
    lonq = np.asarray(nc1.xq) 
    latq = np.asarray(nc1.yq)

    varr0002 = np.asarray(nc0002['SST'][0,:,:])
    varr2 = np.asarray(nc2['SST'][0,:,:])
    ssu2 = np.asarray(nc2['SSU'][0,:,:])
    ssv2 = np.asarray(nc2['SSV'][0,:,:])

if ocean == 'hycom':
    varr000 = np.asarray(nc000['temperature'][0,0,:,:])
    varr = np.asarray(nc['temperature'][0,0,:,:])
    ssu = np.asarray(nc['u_velocity'][0,0,:,:])/100
    ssv = np.asarray(nc['v_velocity'][0,0,:,:])/100
    lon = np.asarray(nc.Longitude)
    lat = np.asarray(nc.Latitude)

lonmin_raw = np.min(lon)
lonmax_raw = np.max(lon)
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

#================================================================
# Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
lon[lon>180] = lon[lon>180] - 360
lon[lon<-180] = lon[lon<-180] + 360
sort_lon = np.argsort(lon)
lon = lon[sort_lon]

# define grid boundaries
lonmin_new = np.min(lon)
lonmax_new = np.max(lon)
latmin = np.min(lat)
latmax = np.max(lat)
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

# Shift central longitude so the Southern Hemisphere and North Indin Ocean domains are plotted continuously
if np.logical_and(lonmax_new >= 90, lonmax_new <=180):
    central_longitude = 90
else:
    central_longitude = -90
print('central longitude: ',central_longitude)

# sort var according to the new longitude
varr0001 = varr0001[:,sort_lon]
varr1 = varr1[:,sort_lon]

varr0002 = varr0002[:,sort_lon]
varr2 = varr2[:,sort_lon]

#================================================================
var_name= 'sst'
units = '($^oC$)'

Rkm=600    # search radius [km]
lns,lts = np.meshgrid(lon,lat)
dummy = np.ones(lns.shape)

lnsq,ltsq = np.meshgrid(lon,latq)
dummy2 = np.ones(lnsq.shape)

#skip=6
skip=4
ln = lns[::skip,::skip]
lt = lts[::skip,::skip]

nhour = int((int(conf['fhhh'][1:])/3))
okfhour = conf['fhhh'][1:] == fhour1
if len(lon_adeck1[okfhour])!=0 and len(lat_adeck1[okfhour])!=0:
    if lat_adeck1[nhour] < (latmax+5.0):
        dR=haversine(lns,lts,lon_adeck1[nhour],lat_adeck1[nhour])/1000.
        dumb=dummy.copy()
        dumb[dR>Rkm]=np.nan
    
        var1 = varr1*dumb
        dvar1 = np.asarray(varr1 - varr0001)*dumb

        var2 = varr2*dumb
        dvar2 = np.asarray(varr2 - varr0002)*dumb
    
        # create figure and axes instances
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
        ax.axis('scaled')
        
        cflevels = np.arange(-1,1.1,0.1)
        cmap = plt.get_cmap('seismic')
        cf = ax.contourf(lon, lat, var2-var1, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
        ax.contour(lon, lat, var2-var1, cflevels, colors='grey',alpha=0.5, linewidths=0.5, transform=ccrs.PlateCarree())
        cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=20, shrink=0.6, extendrect=True, ticks=cflevels[::4])
        cb.ax.tick_params(labelsize=8)
        
        if conf['trackon']=='yes':
            lon_adeck1[np.logical_or(lon_adeck1<lonmin_new,lon_adeck1>lonmax_new)] = np.nan
            ax.plot(lon_adeck1,lat_adeck1,'-ok',markersize=2,alpha=0.4,transform=ccrs.PlateCarree(central_longitude=0))
            ax.plot(lon_adeck2,lat_adeck2,'-ok',markersize=2,alpha=0.4,transform=ccrs.PlateCarree(central_longitude=0))

            ax.plot(lon_adeck1[nhour],lat_adeck1[nhour],'ok',markersize=10,alpha=0.6,markerfacecolor='None',transform=ccrs.PlateCarree(central_longitude=0))
            ax.plot(lon_adeck2[nhour],lat_adeck2[nhour],'ok',markersize=10,alpha=0.6,markerfacecolor='None',transform=ccrs.PlateCarree(central_longitude=0))
       
            if np.logical_and(lon_adeck1[nhour]-5.5 > lonmin_new,lon_adeck1[nhour]+5.5 < lonmax_new):
                mnmx="(min,max)="+"(%6.1f"%np.nanmin(var2-var1)+","+"%6.1f)"%np.nanmax(var2-var1)
                plt.text(lon_adeck1[nhour]-2.15-central_longitude,lat_adeck1[nhour]-4.75,mnmx,fontsize=8,color='DarkOliveGreen',fontweight='bold',bbox=dict(boxstyle="round",color='w',alpha=0.5))
                ax.set_extent([lon_adeck1[nhour]-5.5,lon_adeck1[nhour]+5.5,lat_adeck1[nhour]-5,lat_adeck1[nhour]+5],crs=ccrs.PlateCarree())
            else:
                print('Longitude track limits are out of the ocean domain')
    
        # Add gridlines and labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
        gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}
        
        # Add borders and coastlines
        ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
        model_info = os.environ.get('TITLEgraph','').strip()
        var_info = 'Sea Surface Temperature Difference($^oC$)'
        storm_info = conf['stormName']+conf['stormID']
        title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
        ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
        footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
        ax.text(1.0,-0.08, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)
    
        pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.'+var_name+'.'+conf['fhhh'].lower()+'.png'
        #plt.savefig(pngFile,bbox_inches='tight',dpi=150)
        #plt.close("all")
        
        if ocean == "hycom":
            ssu = ssu[:,sort_lon]
            ssv = ssv[:,sort_lon]
            ssu = ssu[::skip,::skip]
            ssv = ssv[::skip,::skip]

            #q = plt.quiver(ln,lt,ssu,ssv,scale=20,transform=ccrs.PlateCarree())

        if ocean == "mom6":
            lonq[lonq>180] = lonq[lonq>180] - 360
            lonq[lonq<-180] = lonq[lonq<-180] + 360
            sort_lonq = np.argsort(lonq)
            lonq = lonq[sort_lonq]

            ssu1 = ssu1[:,sort_lonq]
            ssv1 = ssv1[:,sort_lon]
            ssu2 = ssu2[:,sort_lonq]
            ssv2 = ssv2[:,sort_lon]

            lnnq,ltth = np.meshgrid(lonq,lat)
            lnnh,lttq = np.meshgrid(lon,latq)

            interpolator1 = RegularGridInterpolator((lat,lonq),ssu1,bounds_error=False, fill_value=None)
            ssu1_interp = interpolator1((lttq,lnnh))

            interpolator2 = RegularGridInterpolator((lat,lonq),ssu2,bounds_error=False, fill_value=None)
            ssu2_interp = interpolator2((lttq,lnnh))

            dvel = np.sqrt((ssu2_interp-ssu1_interp)**2 + (ssv2-ssv1)**2)

            lon_skip = lon[::skip]
            lat_skip = lat[::skip]
            lonq_skip = lonq[::skip]
            latq_skip = latq[::skip]
            ssu1_interp_skip = ssu1_interp[::skip,::skip]
            ssv1_skip = ssv1[::skip,::skip]
            ssu2_interp_skip = ssu2_interp[::skip,::skip]
            ssv2_skip = ssv2[::skip,::skip]

            #q = plt.quiver(lnnh,lttq,ssu_interp,ssv,scale=20,transform=ccrs.PlateCarree())

        dR=haversine(lnsq,ltsq,lon_adeck1[nhour],lat_adeck1[nhour])/1000.
        dumb=dummy2.copy()
        dumb[dR>Rkm]=np.nan

        # create figure and axes instances
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
        ax.axis('scaled')
        
        cflevels = np.arange(-1,1.1,0.1)
        cmap = plt.get_cmap('seismic')
        cf = ax.contourf(lon, latq, dvel*dumb, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
        lb = ax.contour(lon, latq, dvel*dumb, cflevels, colors='grey', alpha=0.5, linewidths=0.5, transform=ccrs.PlateCarree())
        cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=20, shrink=0.6, extendrect=True, ticks=cflevels[::4])
        cb.ax.tick_params(labelsize=8)
            
        q = plt.quiver(lon_skip,latq_skip,ssu2_interp_skip-ssu1_interp_skip,ssv2_skip-ssv1_skip,scale=10,transform=ccrs.PlateCarree())
        
        if conf['trackon']=='yes':
            lon_adeck1[np.logical_or(lon_adeck1<lonmin_new,lon_adeck1>lonmax_new)] = np.nan
            ax.plot(lon_adeck1,lat_adeck1,'-ok',markersize=2,alpha=0.4,transform=ccrs.PlateCarree(central_longitude=0))

            ax.plot(lon_adeck1[nhour],lat_adeck1[nhour],'ok',markersize=10,alpha=0.6,markerfacecolor='None',transform=ccrs.PlateCarree(central_longitude=0))
            ax.plot(lon_adeck2[nhour],lat_adeck2[nhour],'ok',markersize=10,alpha=0.6,markerfacecolor='None',transform=ccrs.PlateCarree(central_longitude=0))
       
            if np.logical_and(lon_adeck1[nhour]-5.5 > lonmin_new,lon_adeck1[nhour]+5.5 < lonmax_new):
                mnmx="(min,max)="+"(%6.1f"%np.nanmin(dvel*dumb)+","+"%6.1f)"%np.nanmax(dvel*dumb)
                plt.text(lon_adeck1[nhour]-2.15-central_longitude,lat_adeck1[nhour]-4.75,mnmx,fontsize=8,color='DarkOliveGreen',fontweight='bold',bbox=dict(boxstyle="round",color='w',alpha=0.5))
                ax.set_extent([lon_adeck1[nhour]-5.5,lon_adeck1[nhour]+5.5,lat_adeck1[nhour]-5,lat_adeck1[nhour]+5],crs=ccrs.PlateCarree())
            else:
                print('Longitude track limits are out of the ocean domain')
    
        # Add gridlines and labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
        gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}
        
        # Add borders and coastlines
        ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
        model_info = os.environ.get('TITLEgraph','').strip()
        var_info = 'Sea Surface Velocity Difference ($m/s$)'
        storm_info = conf['stormName']+conf['stormID']
        title_left = """{0}\n{1}\n{2}""".format(model_info,var_info,storm_info)
        ax.set_title(title_left, loc='left', y=0.99,fontsize=8)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+conf['fhhh'].upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        ax.set_title(title_right, loc='right', y=0.99,fontsize=8)
        footer = os.environ.get('FOOTERgraph','Experimental HAFS Product').strip()
        ax.text(1.0,-0.08, footer, fontsize=8, va="top", ha="right", transform=ax.transAxes)
        
        pngFile = conf['stormName'].upper()+conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.storm.'+var_name+'.change.'+conf['fhhh'].lower()+'.png'
        #plt.savefig(pngFile,bbox_inches='tight',dpi=150)
        #plt.close("all")

else:
    print('There is not latitude or longitude for the center of the storm at this forecast hour. Exiting plotting script')
