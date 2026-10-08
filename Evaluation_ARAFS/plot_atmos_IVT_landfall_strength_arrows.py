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
from matplotlib.patches import FancyArrowPatch
  
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
with open('plot_atmos_IVT_landfall_strength_arrows.yml', 'rt') as f:
    conf = yaml.safe_load(f)

ymdhs = conf['ymdhs']
fhhhs = conf['fhhhs']
curvatures = conf['curvatures']

# Limits for plotting
xlim = conf['xlim']
ylim = conf['ylim']

# Limits to find max IVT before landfall
xlim_landfall = conf['xlim_landfall']
ylim_landfall = conf['ylim_landfall']

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
# Obtain coastline data from Cartopy
coastline_geoms = cfeature.COASTLINE.geometries()

all_lons = []
all_lats = []

for geom in coastline_geoms:
    # Check if geometry is a MultiLineString or LineString
    geoms = geom.geoms if hasattr(geom, 'geoms') else [geom]

    for g in geoms:
        # Extract lon/lat coordinate arrays
        coords = np.array(g.coords)
        lon = coords[:, 0]
        lat = coords[:, 1]

        all_lons.append(lon)
        all_lats.append(lat)

all_lons_arr = np.concatenate(all_lons)
all_lats_arr = np.concatenate(all_lats)

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

ax.plot(xlim[0]+lon_offset,ylim[0],'.',color='k',label='Ralph/CW3E AR Strength Scale')
ax.plot(xlim[0]+lon_offset,ylim[0],'-',color='steelblue',label='Weak: IVT=250-500 $kg m^{-1} s^{-1}$')
ax.plot(xlim[0]+lon_offset,ylim[0],'-',color='orange',label='Moderate: IVT=500-750 $kg m^{-1} s^{-1}$')
ax.plot(xlim[0]+lon_offset,ylim[0],'-',color='crimson',label='Strong: IVT=750-1000 $kg m^{-1} s^{-1}$')
ax.plot(xlim[0]+lon_offset,ylim[0],'-',color='magenta',label='Extreme: IVT=1000-1250 $kg m^{-1} s^{-1}$')
ax.plot(xlim[0]+lon_offset,ylim[0],'-',color='black',label='Exceptional: IVT>1250 $kg m^{-1} s^{-1}$')

ax.legend(loc='lower left',fontsize='x-small')

for f,ymdh in enumerate(ymdhs):
    fname = conf['stormModel'].lower()+'.'+ymdh+'.'+fhhhs[f]+'.grb2'
    grib2file = os.path.join(conf['COMarafs']+ymdh+'/00E/', fname)
    grb = grib2io.open(grib2file,mode='r')

    IVT_u, IVT_v, IVT = compute_ivt(grb,nlat, nlon)

    #IVT[IVT<500] = np.nan
    IVTt = np.copy(IVT)
    IVTt[lon_atm < xlim_landfall[0]+lon_offset] = np.nan
    IVTt[lon_atm > xlim_landfall[1]+lon_offset] = np.nan
    IVTt[lat_atm < ylim_landfall[0]] = np.nan
    IVTt[lat_atm > ylim_landfall[1]] = np.nan

    IVTT = IVTt[np.isfinite(IVTt)]
    lon_ivt = lon_atm[np.isfinite(IVTt)]-180
    lat_ivt = lat_atm[np.isfinite(IVTt)]

    # Find points that are along the coastline
    lon_ivt_max = np.max(lon_ivt)
    lon_ivt_min = np.min(lon_ivt)
    lat_ivt_max = np.max(lat_ivt)
    lat_ivt_min = np.min(lat_ivt)
    oklon = np.logical_and(all_lons_arr < lon_ivt_max,all_lons_arr >= lon_ivt_min)
    all_lons_arr_su = all_lons_arr[oklon]
    all_lats_arr_su = all_lats_arr[oklon]
    oklat = np.logical_and(all_lats_arr_su < lat_ivt_max,all_lats_arr_su >= lat_ivt_min)
    all_lons_arr_sub = all_lons_arr_su[oklat]
    all_lats_arr_sub = all_lats_arr_su[oklat]

    lon_ivt_Sub = []
    lat_ivt_Sub = []
    IVT_Sub = []
    for x in np.arange(len(all_lons_arr_sub)):
        print(x)
        oklon_coast = np.logical_and(np.abs(lon_ivt) >= np.abs(all_lons_arr_sub[x]) - 0.1,np.abs(lon_ivt) < np.abs(all_lons_arr_sub[x]) + 0.1) 
        lon_ivt_su = lon_ivt[oklon_coast]
        lat_ivt_su = lat_ivt[oklon_coast]
        IVT_su = IVTT[oklon_coast]
        #oklat_coast = np.logical_and(np.abs(lat_ivt_su) <= np.abs(all_lats_arr_sub[x]) + 0.1,lat_ivt_su > np.abs(all_lats_arr_sub[x]) - 0.1) 
        oklat_coast = np.logical_and(np.abs(lat_ivt_su) <= np.abs(all_lats_arr_sub[x]) + 1,np.abs(lat_ivt_su) > np.abs(all_lats_arr_sub[x]) - 1) 
        lon_ivt_sub = lon_ivt_su[oklat_coast]
        lat_ivt_sub = lat_ivt_su[oklat_coast]
        IVT_sub = IVT_su[oklat_coast]
        if len(lon_ivt_sub)!=0 and len(lat_ivt_sub)!=0:
            lon_ivt_Sub.append(lon_ivt_sub)
            lat_ivt_Sub.append(lat_ivt_sub)
            IVT_Sub.append(IVT_sub)

    lon_ivt_SUB = np.asarray([item for sublist in lon_ivt_Sub for item in sublist])
    lat_ivt_SUB = np.asarray([item for sublist in lat_ivt_Sub for item in sublist])
    IVT_SUB = np.asarray([item for sublist in IVT_Sub for item in sublist])

    okmax = IVT_SUB == np.max(IVT_SUB)
    IVT_max = np.unique(IVT_SUB[okmax])[0]
    lon_max = np.unique(lon_ivt_SUB[okmax])[0]
    lat_max = np.unique(lat_ivt_SUB[okmax])[0]

    #cf = plt.contourf(lon_atm,lat_atm,IVTt,levels=bounds,colors=colors_ivt(),extend='both',alpha=0.4,transform=transform)

    end_point = (lon_max, lat_max)

    ###############
    # Find start point
    IVT[IVT<500] = np.nan
    IVTtt = np.copy(IVT)
    IVTtt[lon_atm < xlim_landfall[0]+lon_offset] = np.nan
    IVTtt[lon_atm > xlim_landfall[1]+lon_offset] = np.nan
    IVTtt[lat_atm < ylim_landfall[0]] = np.nan
    IVTtt[lat_atm > ylim_landfall[1]] = np.nan

    #IVT_tt = IVTtt[np.isfinite(IVTtt)]
    lon_ivtt = lon_atm[np.isfinite(IVTtt)]-180
    lat_ivtt = lat_atm[np.isfinite(IVTtt)]

    lon_ivt_max = np.max(lon_ivtt)
    lon_ivt_min = np.min(lon_ivtt)
    lat_ivt_max = np.max(lat_ivtt)
    lat_ivt_min = np.min(lat_ivtt)
    ##############

    if lon_ivt_min > xlim_landfall[0]:
        lon_start =  np.mean(lon_ivtt[lat_ivtt<ylim_landfall[0]+1])
        lat_start = lat_ivt_min 
    else:
        lon_start = xlim_landfall[0]
        lat_start = np.mean(lat_ivtt[lon_ivtt<xlim_landfall[0]+1])

    if IVT_max >= 250 and IVT_max< 500:
        color = 'steelblue'
        label = 'Weak'
        lon_start = lon_start + 2
        lat_start = lat_start + 2
    if IVT_max >= 500 and IVT_max< 750:
        color = 'orange'
        label = 'Moderate'
        lon_start = lon_start + 1
        lat_start = lat_start + 1
    if IVT_max >= 750 and IVT_max< 1000:
        color = 'crimson'
        label = 'Strong'
    if IVT_max >= 1000 and IVT_max< 1250:
        color = 'magenta'
        label = 'Extreme'
    if IVT_max >= 1250:
        color = 'black'
        label = 'Exceptional'

    start_point = (lon_start, lat_start)

    ax.text(lon_start-1,lat_start-1,str(f+1),color=color,fontweight='bold',fontsize=14,transform=ccrs.PlateCarree())

    if curvatures[f] == 'up':
        connectionstyle="arc3,rad=0.2"
    if curvatures[f] == 'down':
        connectionstyle="arc3,rad=-0.2"

    # Create the curved arrow patch
    arrow = FancyArrowPatch(
    posA=start_point,
    posB=end_point,
    connectionstyle=connectionstyle,  # Positive rad curves left/up; negative curves right/down
    arrowstyle="-|>",               # Arrowhead style
    mutation_scale=20,              # Arrowhead size
    color=color,
    linewidth=2.5,
    transform=ccrs.PlateCarree()    # Tells Cartopy coordinates are in Longitude/Latitude
    )

    # Add arrow to plot
    ax.add_patch(arrow)

#pngFile = conf['stormID'].upper()+'.'+conf['ymdh']+'.'+conf['stormModel']+'.ocean.'+var_name+'.'+conf['fhhh'].lower()+'.png'
#plt.savefig(pngFile,bbox_inches='tight',dpi=150)
