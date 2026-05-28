#!/usr/bin/env python3

"""This script is to plot the IVT time series."""

import os
import sys
import logging
import math
import datetime

import yaml
import numpy as np
import pandas as pd

import grib2io
from netCDF4 import Dataset

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

import time

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,6)

"""
shfl_mean = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
shfl_mean[:] = np.nan
shfl_max = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
shfl_max[:] = np.nan
shfl_min = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
shfl_min[:] = np.nan
shfl_std = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
shfl_std[:] = np.nan
lhfl_mean = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
lhfl_mean[:] = np.nan
lhfl_max = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
lhfl_max[:] = np.nan
lhfl_min = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
lhfl_min[:] = np.nan
lhfl_std = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
lhfl_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 
    print(conf['COMmodels'][exp])
    for c,ymdh in enumerate(conf['ymdhs']):
        print(ymdh)
        for f,fh in enumerate(ff):
            start = time.perf_counter()
    
            if len(str(fh))==3:
                fhour = str(fh)
            if len(str(fh))==2:
                fhour = '0'+str(fh)
            if len(str(fh))==1:
                fhour = '00'+str(fh)
    
            print(fhour)
            fname = conf['stormModel'].lower()+'.'+ymdh+'.f'+fhour+'.grb2'
     
            grib2file = os.path.join(conf['COMmodels'][exp]+ymdh+'/00E/', fname)
            print(f'grib2file: {grib2file}')
            grb = grib2io.open(grib2file,mode='r')
        
            lat = grb.select(shortName='NLAT')[0].data
            lon = grb.select(shortName='ELON')[0].data
         
            nlat = lat.shape[0]
            nlon = lat.shape[1]
    
            print('Extracting SHTFL at surface')
            shfl = grb.select(shortName='SHTFL')[0].data

            print('Extracting LHTFL at surface')
            lhfl = grb.select(shortName='LHTFL')[0].data
      
            shfl_mean[exp,c,f] = np.nanmean(shfl)
            shfl_max[exp,c,f] = np.nanmax(shfl)
            shfl_min[exp,c,f] = np.nanmin(shfl)
            shfl_std[exp,c,f] = np.nanstd(shfl)

            lhfl_mean[exp,c,f] = np.nanmean(lhfl)
            lhfl_max[exp,c,f] = np.nanmax(lhfl)
            lhfl_min[exp,c,f] = np.nanmin(lhfl)
            lhfl_std[exp,c,f] = np.nanstd(lhfl)

            end = time.perf_counter()
            print(f"Elapsed time: {end - start:.6f} seconds")

np.savez('shfl_lhfl.npz', shfl_mean=shfl_mean, shfl_max=shfl_max, shfl_min=shfl_min, shfl_std=shfl_std,lhfl_mean=lhfl_mean, lhfl_max=lhfl_max, lhfl_min=lhfl_min, lhfl_std=lhfl_std)
"""

loaded_data = np.load('shfl_lhfl.npz')
shfl_mean = loaded_data['shfl_mean']
shfl_max = loaded_data['shfl_max']
shfl_min = loaded_data['shfl_min']
shfl_std = loaded_data['shfl_std']
lhfl_mean = loaded_data['lhfl_mean']
lhfl_max = loaded_data['lhfl_max']
lhfl_min = loaded_data['lhfl_min']
lhfl_std = loaded_data['lhfl_std']

########
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]

fig,ax = plt.subplots(figsize = (13,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,shfl_mean[exp,c,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=10)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('shfl ($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Sensible Heat Flux ' + conf['ymdhs'][c],fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

fig,ax = plt.subplots(figsize = (13,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,lhfl_mean[exp,c,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=10)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('lhfl ($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Latent Heat Flux ' + conf['ymdhs'][c],fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

##############
shfl_diff = np.diff(shfl_mean,axis=0)[0,:,:]

fig,ax = plt.subplots(figsize = (13,7))
for c,ymdh in enumerate(conf['ymdhs']):
    plt.plot(ff,shfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('shfl Diff($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Sensible Heat Flux (Coupled - Uncoupled)',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.savefig('shfl_couple_minus_uncoupled.png')

##############
lhfl_diff = np.diff(lhfl_mean,axis=0)[0,:,:]

fig,ax = plt.subplots(figsize = (13,7))
for c,ymdh in enumerate(conf['ymdhs']):
    plt.plot(ff,lhfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('lhfl Diff($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Latent Heat Flux (Coupled - Uncoupled)',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.savefig('lhfl_couple_minus_uncoupled.png')

##############
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]

fig,ax = plt.subplots(figsize = (13,7))
plt.plot(ff,shfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('shfl Diff($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Sensible (Coupled - Uncoupled)',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

##############
fig,ax = plt.subplots(figsize = (13,7))
plt.plot(ff,lhfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('lhfl Diff($W/m^2$)',fontsize=18,labelpad=10)
plt.title('Mean Latent (Coupled - Uncoupled)',fontsize=18)
plt.xlim([0,126])


