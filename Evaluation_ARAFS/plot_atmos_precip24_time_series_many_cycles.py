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

np.savez('precip24_3.npz', apcp24_mean=apcp24_mean, apcp24_max=apcp24_max, apcp24_min=apcp24_min, apcp24_std=apcp24_std)
'''

loaded_data1 = np.load('precip24.npz')
apcp24_mean1 = loaded_data1['apcp24_mean']
ymdhs1 = ['2023010600','2023010700','2023010800','2023010900','2023011000','2023030700','2023030800','2023030900','2023031000','2023031100','2024021600','2024021700','2024021800','2024021900','2024022000','2025020100','2025020200','2025020300','2025020400','2025020500','2025022000','2025022100','2025022200','2025022300','2025022400','2026010300','2026010400','2026010500','2026010600','2026010700']

loaded_data2 = np.load('precip24_2.npz')
apcp24_mean2 = loaded_data2['apcp24_mean']
ymdhs2 = ['2023010500','2023030600','2024021500','2025013100','2025021900','2025120400', '2025120500','2025120600','2025120700','2025120800','2025120900','2025121000']

loaded_data3 = np.load('precip24_3.npz')
apcp24_mean3 = loaded_data3['apcp24_mean']
ymdhs3 = ['2025121100','2025121500','2025121600','2025121700','2025121800','2025121900','2025122000','2025122100','2025122200','2025122300','2025122400']

apcp24_mean0 = np.hstack([apcp24_mean1,apcp24_mean2,apcp24_mean3])
ymdhs0 = np.hstack([ymdhs1,ymdhs2,ymdhs3])

apcp24_mean = apcp24_mean0 * np.nan
ymdhs_mean = []
for c,ymdh in enumerate(conf['ymdhs']):
    print(ymdh)
    okc = np.where(ymdhs0 == ymdh)[0][0]
    apcp24_mean[:,c,:] = apcp24_mean0[:,okc,:]
    ymdhs_mean.append(ymdhs0[okc])

np.savez('precip24_53_cycles.npz',apcp24_mean=apcp24_mean,ymdhs_mean=ymdhs_mean)

loaded_data = np.load('precip24_53_cycles.npz')
apcp24_mean = loaded_data['apcp24_mean']
ymdhs_mean = loaded_data['ymdhs_mean']

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
    plt.plot(ff,apcp24_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10)
for c in [0,6,12,19,24,30,38,48]:
    plt.plot(ff,apcp24_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c][0:-4])
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



