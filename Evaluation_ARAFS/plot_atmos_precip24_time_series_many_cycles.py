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

'''
apcp24_mean = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
apcp24_mean[:] = np.nan
apcp24_max = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
apcp24_max[:] = np.nan
apcp24_min = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
apcp24_min[:] = np.nan
apcp24_std = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
apcp24_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 
    print(conf['COMmodels'][exp])
    for c,ymdh in enumerate(conf['ymdhs']):
        print(ymdh)
        start = time.perf_counter()
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
    
                fname = conf['stormModel'].lower()+'.'+ ymdh+'.'+fhhhs[0]+'.grb2'
                grib2file = os.path.join(conf['COMmodels'][exp]+ymdh+'/00E/',fname)
                grb = grib2io.open(grib2file,mode='r')
                lat = grb.select(shortName='NLAT')[0].data
    
                apcp = np.empty((len(fhhhs),lat.shape[0],lat.shape[1]))
                apcp[:] = np.nan
                for f,fhhh in enumerate(fhhhs):
                    fname = conf['stormModel'].lower()+'.'+ymdh+'.'+fhhh+'.grb2'
                    grib2file = os.path.join(conf['COMmodels'][exp]+ymdh+'/00E/',fname)
                    print(grib2file)
                    grb = grib2io.open(grib2file,mode='r')
                    apcp[f,:,:] = grb.select(shortName='APCP')[0].data*0.0393701  # convert kg/m^2 to in

                apcp_24 = np.sum(apcp,axis=0)
    
                print('Obtaining apcp24_mean, apcp24_std, apcp24_max, apcp24_min')
                apcp24_mean[exp,c,F] = np.nanmean(apcp_24)
                apcp24_std[exp,c,F] = np.nanstd(apcp_24)
                apcp24_max[exp,c,F] = np.nanmax(apcp_24)
                apcp24_min[exp,c,F] = np.nanmin(apcp_24)

        end = time.perf_counter()
        print(f"Elapsed time: {end - start:.6f} seconds")

np.savez('precip24_2026.npz', apcp24_mean=apcp24_mean, apcp24_max=apcp24_max, apcp24_min=apcp24_min, apcp24_std=apcp24_std)
'''
'''
loaded_data_2023 = np.load('precip24_2023.npz')
apcp24_mean_2023 = loaded_data_2023['apcp24_mean']
loaded_data_2024 = np.load('precip24_2024.npz')
apcp24_mean_2024 = loaded_data_2024['apcp24_mean']
loaded_data_2025 = np.load('precip24_2025.npz')
apcp24_mean_2025 = loaded_data_2025['apcp24_mean']
loaded_data_2026 = np.load('precip24_2026.npz')
apcp24_mean_2026 = loaded_data_2026['apcp24_mean']
#apcp24_max = loaded_data['apcp24_max']
#apcp24_min = loaded_data['apcp24_min']
#apcp24_std = loaded_data['apcp24_std']

apcp24_mean = np.hstack([apcp24_mean_2023,apcp24_mean_2024,apcp24_mean_2025,apcp24_mean_2026])

np.savez('precip24.npz', apcp24_mean=apcp24_mean)
'''

loaded_data = np.load('precip24.npz')
apcp24_mean = loaded_data['apcp24_mean']/0.0393701 # Convert inches to mm

########
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]

fig,ax = plt.subplots(figsize = (13,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,apcp24_mean[exp,c,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=10)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('Precip24 (mm)',fontsize=18,labelpad=10)
plt.title('Mean Accumulated 24 Hours Precipitation ' + conf['ymdhs'][c],fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

##############
apcp24_diff = np.diff(apcp24_mean,axis=0)[0,:,:]

fig,ax = plt.subplots(figsize = (13,7))
plt.plot(ff,apcp24_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('Precip24 Diff (mm)',fontsize=18,labelpad=10)
plt.title('Mean Accumulated 24 Hours Precipitation Coupled - Uncoupled ',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

##############
fig,ax = plt.subplots(figsize = (13,7))
for c,ymdh in enumerate(conf['ymdhs']):
    plt.plot(ff,apcp24_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('Precip24 Diff (mm)',fontsize=18,labelpad=10)
plt.title('Mean Accumulated 24 Hours Precipitation Coupled - Uncoupled ',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.savefig('precip24_coupled_minus_uncoupled.png')



