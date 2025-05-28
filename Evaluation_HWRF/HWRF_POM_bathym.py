#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:52:43 2020

@author: aristizabal
"""

#%%
storm_id = '95l'
cycle = '2020081018'

# MAB
lon_lim = [-77,-52]
lat_lim = [35,46]


#%%
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import os
import os.path
import glob
import cmocean

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

folder_hwrf_pom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF_POM_' + storm_id + '_' + str(ti.year) + '/' + 'HWRF_POM_' + storm_id + '_' + cycle + '/'

#%% Reading POM grid file
grid_file = sorted(glob.glob(os.path.join(folder_hwrf_pom,'*grid*.nc')))[0]
pom_grid = xr.open_dataset(grid_file)
lon_pom = np.asarray(pom_grid['east_e'][:])
lat_pom = np.asarray( pom_grid['north_e'][:])
zlevc = np.asarray(pom_grid['zz'][:])
topoz = np.asarray(pom_grid['h'][:])

#%% bathymetry POM whole domain
kw = dict(levels = np.arange(-8000,1,500))

plt.figure(figsize = (10,5))
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_pom,lat_pom,-topoz,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(np.min(lon_pom),np.max(lon_pom))
plt.ylim(np.min(lat_pom),np.max(lat_pom))
plt.title('POM Bathymetry')

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

#%% Bathymetry POM MAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_pom,lat_pom,-topoz,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-76,-70)
plt.ylim(35,42)
plt.title('POM Bathymetry')

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

#%% Bathymetry POM SAB
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_pom,lat_pom,-topoz,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-85,-76)
plt.ylim(25,35)
plt.title('POM Bathymetry')

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

#%% Bathymetry POM GoMex
kw = dict(levels = np.arange(-200,1,20))

plt.figure()
plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
plt.contourf(lon_pom,lat_pom,-topoz,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-98,-80)
plt.ylim(18,31)
plt.title('POM Bathymetry')

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
plt.contourf(lon_pom,lat_pom,-topoz,cmap=cmocean.cm.topo,**kw)
plt.colorbar()
plt.axis('scaled')
plt.xlim(-90,-60)
plt.ylim(10,30)
plt.title('POM Bathymetry')
