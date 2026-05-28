#!/usr/bin/env python3

#exp_names = ['HFSAv2a_baseline_latest','HKPP','HKP2']
#exp_labels = ['HFSAv2a_baseline','HKPP','HKP2']
#exp_colors = ['orange','green','dodgerblue']

"""This script is to plot out HAFS atmospheric sensible heat flux and 10-m wind."""

import os
import sys
import logging
import math
import datetime

import yaml
import numpy as np
import pandas as pd
#from scipy.ndimage import gaussian_filter

import grib2io
from netCDF4 import Dataset

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

#from geo4HYCOM import haversine

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,3)
shfl_mean = np.empty((len(conf['COMmodels']),len(ff)))
shfl_mean[:] = np.nan
shfl_max = np.empty((len(conf['COMmodels']),len(ff)))
shfl_max[:] = np.nan
shfl_min = np.empty((len(conf['COMmodels']),len(ff)))
shfl_min[:] = np.nan
shfl_std = np.empty((len(conf['COMmodels']),len(ff)))
shfl_std[:] = np.nan
lhfl_mean = np.empty((len(conf['COMmodels']),len(ff)))
lhfl_mean[:] = np.nan
lhfl_max = np.empty((len(conf['COMmodels']),len(ff)))
lhfl_max[:] = np.nan
lhfl_min = np.empty((len(conf['COMmodels']),len(ff)))
lhfl_min[:] = np.nan
lhfl_std = np.empty((len(conf['COMmodels']),len(ff)))
lhfl_std[:] = np.nan
totfl_mean = np.empty((len(conf['COMmodels']),len(ff)))
totfl_mean[:] = np.nan
totfl_max = np.empty((len(conf['COMmodels']),len(ff)))
totfl_max[:] = np.nan
totfl_min = np.empty((len(conf['COMmodels']),len(ff)))
totfl_min[:] = np.nan
totfl_std = np.empty((len(conf['COMmodels']),len(ff)))
totfl_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    #================================================================
    # Read fluxes
    for f,fh in enumerate(ff):

        if len(str(fh))==3:
            fhour = str(fh)
        if len(str(fh))==2:
            fhour = '0'+str(fh)
        if len(str(fh))==1:
            fhour = '00'+str(fh)

        fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.f'+fhour+'.grb2'
        grib2file = os.path.join(conf['COMmodels'][exp]+conf['ymdh']+'/00E', fname)

        print(f'grib2file: {grib2file}')
        grb = grib2io.open(grib2file,mode='r')
    
        print('Extracting lat, lon')
        lat = np.asarray(grb.select(shortName='NLAT')[0].data)
        lon = np.asarray(grb.select(shortName='ELON')[0].data)
     
        # Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
        lon[lon>180] = lon[lon>180] - 360
        #sort_lon = np.argsort(lon)
        #lon = lon[sort_lon]
    
        # define grid boundaries
        lonmin_new = np.min(lon)
        lonmax_new = np.max(lon)
        latmin = np.min(lat)
        latmax = np.max(lat)
        print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
    
        print('Extracting SHTFL at surface')
        shfl = grb.select(shortName='SHTFL')[0].data
        #shf_sort = shf[:,sort_lon]
        #shf = gaussian_filter(shf, 5)

        print('Extracting LHTFL at surface')
        lhfl = grb.select(shortName='LHTFL')[0].data

        #oklon = np.logical_and(lon>=xlim[0],lon<xlim[1]) 
        #oklat = np.logical_and(lat>=ylim[0],lat<ylim[1]) 

        #lonsub = lon[oklon]
        #lonsub = lon[oklon]
    
        shfl_mean[exp,f] = np.nanmean(shfl)
        shfl_max[exp,f] = np.nanmax(shfl)
        shfl_min[exp,f] = np.nanmin(shfl)
        shfl_std[exp,f] = np.nanstd(shfl)

        lhfl_mean[exp,f] = np.nanmean(lhfl)
        lhfl_max[exp,f] = np.nanmax(lhfl)
        lhfl_min[exp,f] = np.nanmin(lhfl)
        lhfl_std[exp,f] = np.nanstd(lhfl)

        totfl = shfl + lhfl
        totfl_mean[exp,f] = np.nanmean(totfl)
        totfl_max[exp,f] = np.nanmax(totfl)
        totfl_min[exp,f] = np.nanmin(totfl)
        totfl_std[exp,f] = np.nanstd(totfl)

'''
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])): 
    plt.plot(ff,shfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,shfl_min[exp,:],shfl_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(ff,shfl_max[exp,:],'.-',color=conf['exp_colors'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Sensible Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Sensible Heat Flux Footprint (200 km around storm eye))',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))
'''

########
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,shfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Sensible Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Sensible Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,shfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    #ax.fill_between(ff,shfl_mean[exp,:]-shfl_std[exp,:],shfl_mean[exp,:]+shfl_std[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    ax.fill_between(ff,shfl_min[exp,:],shfl_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Sensible Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Sensible Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,lhfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Latent Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Latent Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

#########

fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,lhfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    #ax.fill_between(ff,lhfl_mean[exp,:]-lhfl_std[exp,:],lhfl_mean[exp,:]+lhfl_std[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    ax.fill_between(ff,lhfl_min[exp,:],lhfl_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Latent Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Latent Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

#############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,totfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Enthalpy Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Enthalpy Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))


#############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,totfl_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    #ax.fill_between(ff,totfl_mean[exp,:]-totfl_std[exp,:],totfl_mean[exp,:]+totfl_std[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    ax.fill_between(ff,totfl_min[exp,:],totfl_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('Enthalpy Heat Flux ($W/m^2$)',fontsize=14,labelpad=10)
plt.title('Enthalpy Heat Flux',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

