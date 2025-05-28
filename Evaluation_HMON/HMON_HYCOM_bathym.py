#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 13:40:33 2020

@author: aristizabal
"""

#%%
storm_id = '94l'
cycle = '2020080312'

# MAB
lon_lim = [-77,-52]
lat_lim = [35,46]

#%% 
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import os.path
import glob
import cmocean
import xarray as xr

import sys
sys.path.append('/home/Maria.Aristizabal/NCEP_scripts/')
from utils4HYCOM import readgrids, readdepth
#from utils4HYCOM import readdepth, readVar
#from utils4HYCOM2 import readBinz

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%%

ti = datetime.today()
ffig = '/home/Maria.Aristizabal/Figures/'+ str(ti.year) + '/' + ti.strftime('%b-%d') 
folder_fig =  ffig + '/' + storm_id + '_' + cycle + '/'

os.system('mkdir ' +  ffig)
os.system('mkdir ' +  folder_fig)

#%% Bathymetry file

bath_file = '/scratch2/NOS/nosofs/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])
bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#%% folder and file names    
ti = datetime.today() - timedelta(1)

folder_hmon_hycom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HMON_HYCOM_' + storm_id + '_' + str(ti.year) + '/' + 'HMON_HYCOM_' + storm_id + '_' + cycle + '/'

#%% Reading RTOFS grid    
grid_file = sorted(glob.glob(os.path.join(folder_hmon_hycom,'*regional.grid.*')))[0][:-2]

depth_file = sorted(glob.glob(os.path.join(folder_hmon_hycom,'*regional.depth.*')))[0][:-2]                                                

#%% Reading RTOFS grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
#lines_grid = [line.rstrip() for line in open(grid_file+'.b')]
lon_hycom = np.array(readgrids(grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(grid_file,'plat:',[0]))

bathy_hycom = np.asarray(readdepth(depth_file,'depth'))

bathy_hycom[bathy_hycom >= 10**30] = np.nan

#%% bathymetry HYCOM whole domain
kw = dict(levels = np.arange(-8000,1,500))

plt.figure(figsize = (10,5))
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,-bathy_hycom,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(np.min(lon_hycom-360),np.max(lon_hycom-360))
plt.ylim(np.min(lat_hycom),np.max(lat_hycom))
plt.title('HYCOM Bathymetry')

#%% Bathymetry GEBCO HYCOM domain
kw = dict(levels = np.arange(-8000,1,500))

plt.figure(figsize = (10,5))
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.title('GEBCO Bathymetry')

#%% Bathymetry GEBCO MAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-76,-70)
plt.ylim(35,42)
plt.title('GEBCO Bathymetry')

#%% Bathymetry HYCOM MAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,-bathy_hycom,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-76,-70)
plt.ylim(35,42)
plt.title('HYCOM Bathymetry')

#%% Bathymetry GEGCO SAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-85,-76)
plt.ylim(25,35)
plt.title('GEBCO Bathymetry')

#%% Bathymetry HYCOM SAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,-bathy_hycom,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-85,-76)
plt.ylim(25,35)
plt.title('HYCOM Bathymetry')

#%% Bathymetry GEGCO GoMex
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-98,-80)
plt.ylim(18,31)
plt.title('GEBCO Bathymetry')

#%% Bathymetry HYCOM GoMex
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,-bathy_hycom,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-98,-80)
plt.ylim(18,31)
plt.title('HYCOM Bathymetry')

#%% Bathymetry GEGCO Caribbean
kw = dict(levels = np.arange(-8000,1,500))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-90,-60)
plt.ylim(10,30)
plt.title('GEBCO Bathymetry')

#%% Bathymetry POM Caribbean
kw = dict(levels = np.arange(-8000,1,500))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,-bathy_hycom,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-90,-60)
plt.ylim(10,30)
plt.title('HYCOM Bathymetry')

