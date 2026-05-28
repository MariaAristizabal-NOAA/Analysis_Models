#!/usr/bin/env python3

"""This script is to plot the IWV time series."""

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

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,3)
apcp24_mean = np.empty((len(conf['COMmodels']),len(ff)))
apcp24_mean[:] = np.nan
apcp24_max = np.empty((len(conf['COMmodels']),len(ff)))
apcp24_max[:] = np.nan
apcp24_min = np.empty((len(conf['COMmodels']),len(ff)))
apcp24_min[:] = np.nan
apcp24_std = np.empty((len(conf['COMmodels']),len(ff)))
apcp24_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    print('Extracting 24h accumulate precipitation at surface')
    #================================================================
    for F,fh in enumerate(ff):

        if len(str(fh))==3:
            fhour = 'f'+str(fh)
        if len(str(fh))==2:
            fhour = 'f0'+str(fh)
        if len(str(fh))==1:
            fhour = 'f00'+str(fh)
        
        if int(fhour[1:]) - 24 >= 0:
            print('Extracting 24h accumulate precipitation at surface, forecast hour '+fhour)
            fhhs = np.arange(int(fhour[1:])-24,int(fhour[1:])+1,3)
        
            fhhhs = []
            for f in fhhs:
                if f<10: fhhhs.append('f00'+str(f))
                if f>10 and f<100: fhhhs.append('f0'+str(f))
                if f>100: fhhhs.append('f'+str(f))
                
            fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+fhhhs[0]+'.grb2'
            grib2file = os.path.join(conf['COMmodels'][exp]+conf['ymdh']+'/00E/', fname)
            grb = grib2io.open(grib2file,mode='r')
            lat = grb.select(shortName='NLAT')[0].data
        
            apcp = np.empty((len(fhhhs),lat.shape[0],lat.shape[1]))
            apcp[:] = np.nan
            for f,fhhh in enumerate(fhhhs):
                fname = conf['stormModel'].lower()+'.'+ conf['ymdh']+'.'+fhhh+'.grb2'
                grib2file = os.path.join(conf['COMmodels'][exp]+conf['ymdh']+'/00E/', fname)
                print(grib2file)
                grb = grib2io.open(grib2file,mode='r')
                apcp[f,:,:] = grb.select(shortName='APCP')[0].data*0.0393701  # convert kg/m^2 to in
        
            apcp_24 = np.sum(apcp,axis=0)

            print('Obtaining apcp24_mean, apcp24_std, apcp24_max, apcp24_min')
            apcp24_mean[exp,F] = np.nanmean(apcp_24) 
            print(apcp24_mean[exp,F])
            apcp24_std[exp,F] = np.nanstd(apcp_24) 
            print(apcp24_std[exp,F])
            apcp24_max[exp,F] = np.nanmax(apcp_24) 
            print(apcp24_max[exp,F])
            apcp24_min[exp,F] = np.nanmin(apcp_24) 
            print(apcp24_min[exp,F])

########
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,apcp24_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('24h Acc. Precip ($inches$)',fontsize=14,labelpad=10)
plt.title('Mean 24h Accumulated Precipitation ' + conf['ymdh'],fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,apcp24_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,apcp24_min[exp,:],apcp24_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(ff,apcp24_min[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
    plt.plot(ff,apcp24_max[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('24h Acc. Precip ($inches$)',fontsize=14,labelpad=10)
plt.title('Mean, Max and Min 24h Accumulated Precipitation ' + conf['ymdh'],fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,apcp24_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,apcp24_mean[exp,:]-apcp24_std[exp,:],apcp24_mean[exp,:]+apcp24_std[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(ff,apcp24_mean[exp,:]-apcp24_std[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
    plt.plot(ff,apcp24_mean[exp,:]+apcp24_std[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('24h Acc. Precip ($inches$)',fontsize=14,labelpad=10)
plt.title('Mean and Std 24h Accumulated Precipitation ' + conf['ymdh'],fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

