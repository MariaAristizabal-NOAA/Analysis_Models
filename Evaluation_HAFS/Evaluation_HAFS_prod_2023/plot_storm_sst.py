"""

 plot_storm_sst.py
 -------------
    read a HYCOM 3z .nc file,
    extract footprint SST and plot in time series (R<=500km)


 ************************************************************************
 usage: python plot_storm_sst.py stormModel stormName stormID YMDH trackon COMhafs graphdir
 -----
 ************************************************************************


 HISTORY
 -------
    modified to implement new filenames and hycom domains, as well as
        improve graphics -JS & MA 06/2022
    modified to comply the convention of number of input argument and 
       graphic filename. -hsk 8/2020
    modified to take global varibles from kick_graphics.py -hsk 9/20/2018
    modified to fit for RT run by Hyun-Sook Kim 5/17/2017
    edited by Hyun-Sook Kim 9/18/2015
    modified by Hyun-Sook Kim 11/18/2016
---------------------------------------------------------------
"""

from utils4HWRF import readTrack6hrly
from geo4HYCOM import haversine

import os
import sys
import glob
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#plt.switch_backend('agg')

#================================================================
model = "hfsa"
storm = "lee"
tcid = "13l"
cycle = "2023090706"
trackon = "yes"
COMOUT = "/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HFSAv1p1_HYCOM/2023090706/13l/"
graphdir = "./"
fhour = 'f120'
cartopyDataDir = '/scratch1/NCEPDEV/hwrf/noscrub/local/share/cartopy'

#model = sys.argv[1]
#storm = sys.argv[2]
#tcid = sys.argv[3]
#cycle = sys.argv[4]
#trackon = sys.argv[5]
#COMOUT = sys.argv[6]
#graphdir = sys.argv[7]

if not os.path.isdir(graphdir):
      p=Path(graphdir)
      p.mkdir(parents=True)

print("code:   plot_storm_sst.py")

atcf = COMOUT + tcid + '.' + cycle + '.' + model + '.trak.atcfunix'

adt,aln,alt,pmn,vmx = readTrack6hrly(atcf)

# Set Cartopy data_dir location
cartopy.config['data_dir'] = cartopyDataDir

#   ------------------------------------------------------------------------------------
Rkm=500    # search radius [km]

# - get SST  *_3z_*.[nc] files
#afiles = sorted(glob.glob(os.path.join(COMOUT,'*3z*.nc')))
afile = COMOUT + tcid + '.' + cycle + '.' + model + '.hycom.3z.' + fhour + '.nc'

afile0 = COMOUT + tcid + '.' + cycle + '.' + model + '.hycom.3z.f000' + '.nc'
ncfile0 = xr.open_dataset(afile0)

var0 = ncfile0['temperature'].isel(Z=0)
lon = np.asarray(ncfile0.Longitude)
lat = np.asarray(ncfile0.Latitude)

var_name = 'sst'
units = '($^oC$)'

lns,lts = np.meshgrid(lon,lat)
skip = 6
ln = lns[::skip,::skip]
lt = lts[::skip,::skip]
dummy = np.ones(lns.shape)

k = int((int(fhour[1:])/6))
if alt[k] < (np.max(lat)+5.0):
    dR=haversine(lns,lts,aln[k],alt[k])/1000.
    dumb=dummy.copy()
    dumb[dR>Rkm]=np.nan
    ncfile = xr.open_dataset(afile)
    varr = ncfile['temperature'].isel(Z=0)
    var = np.asarray(varr[0])*dumb
    dvar = np.asarray(varr[0]-np.squeeze(var0))*dumb
    u0=ncfile['u_velocity'].isel(Z=0)
    v0=ncfile['v_velocity'].isel(Z=0)
    u0=np.squeeze(u0)[::skip,::skip]
    v0=np.squeeze(v0)[::skip,::skip]
    #define forecast hour
    fhr=k*6
    
    # create figure and axes instances
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.axis('scaled')
    
    cflevels = np.linspace(19, 32, 27)
    cmap = plt.get_cmap('turbo')
    cf = ax.contourf(lon, lat, var, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    ax.contour(lon, lat, var, cflevels, colors='grey',alpha=0.5, transform=ccrs.PlateCarree())
    cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.75, extendrect=True, ticks=cflevels[::4])
    cb.ax.tick_params(labelsize=8)
    if trackon[0].lower()=='y':
        plt.plot(aln,alt,'-ok',linewidth=3,alpha=0.6,markersize=2)
        plt.plot(aln[k],alt[k],'ok',markerfacecolor='none',markersize=10,alpha=0.6)
    mnmx="(min,max)="+"(%6.1f"%np.nanmin(var)+","+"%6.1f)"%np.nanmax(var)
    plt.text(aln[k]-2.15,alt[k]-4.75,mnmx,fontsize=8,color='DarkOliveGreen',fontweight='bold',bbox=dict(boxstyle="round",color='w',alpha=0.5))
    plt.axis([aln[k]-5.5,aln[k]+5.5,alt[k]-5,alt[k]+5])
    q=plt.quiver(ln,lt,u0,v0,scale=2000)
           # Add gridlines and labels
    #gl = ax.gridlines(crs=transform, draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
    gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    
    # Add borders and coastlines
    #ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='whitesmoke')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    title_center = 'Sea Surface Temperature (${^o}$C), Currents'
    ax.set_title(title_center, loc='center', y=1.05, fontsize=8)
    title_left = model.upper()+' '+storm.upper()+tcid.upper()
    ax.set_title(title_left, loc='left', fontsize=8)
    title_right = 'Init: '+cycle+'Z '+'F'+"%03d"%(fhr)
    ax.set_title(title_right, loc='right', fontsize=8)
    
    pngFile=os.path.join(graphdir,storm.upper()+tcid.upper()+'.'+cycle+'.'+model.upper()+'.ocean.storm.'+var_name+'.f'+"%03d"%(fhr)+'.png')
    plt.savefig(pngFile,bbox_inches='tight',dpi=150)
    #plt.close("all")
           
    # create figure and axes instances for change plot
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.axis('scaled')
    cflevels = np.linspace(-4, 4, 17)
    cmap = plt.get_cmap('RdBu_r')
    cf = ax.contourf(lon, lat, dvar, levels=cflevels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    ax.contour(lon, lat, dvar, cflevels,colors='grey',alpha=0.5, transform=ccrs.PlateCarree())
    cb = plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.75, extendrect=True, ticks=cflevels[::4])
    cb.ax.tick_params(labelsize=8)
    if trackon[0].lower()=='y':
        plt.plot(aln,alt,'-ok',linewidth=3,alpha=0.6,markersize=2)
        plt.plot(aln[k],alt[k],'ok',markerfacecolor='none',markersize=10,alpha=0.6)
    mnmx="(min,max)="+"(%6.1f"%np.nanmin(dvar)+","+"%6.1f)"%np.nanmax(dvar)
    plt.text(aln[k]-2.15,alt[k]-4.75,mnmx,fontsize=8,color='DarkOliveGreen',fontweight='bold',bbox=dict(boxstyle="round",color='w',alpha=0.5))
    plt.axis([aln[k]-5.5,aln[k]+5.5,alt[k]-5,alt[k]+5])
           # Add gridlines and labels
    #gl = ax.gridlines(crs=transform, draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
    gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
           # Add borders and coastlines
    #ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='whitesmoke')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    title_center = 'Sea Surface Temperature Change (${^o}$C)'
    ax.set_title(title_center, loc='center', y=1.05, fontsize=8)
    title_left = model.upper()+' '+storm.upper()+tcid.upper()
    ax.set_title(title_left, loc='left', fontsize=8)
    title_right = 'Init: '+cycle+'Z '+'F'+"%03d"%(fhr)
    ax.set_title(title_right, loc='right', fontsize=8)
    pngFile=os.path.join(graphdir,storm.upper()+tcid.upper()+'.'+cycle+'.'+model.upper()+'.ocean.storm.'+var_name+'.change.f'+"%03d"%(fhr)+'.png')
    plt.savefig(pngFile,bbox_inches='tight',dpi=150)
    #plt.close("all")
    
