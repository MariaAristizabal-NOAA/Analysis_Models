#!/usr/bin/env python3

"""This script is to plot the sensible and latent heat fluxes time series."""

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

'''
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

np.savez('shfl_lhfl_2.npz', shfl_mean=shfl_mean, shfl_max=shfl_max, shfl_min=shfl_min, shfl_std=shfl_std,lhfl_mean=lhfl_mean, lhfl_max=lhfl_max, lhfl_min=lhfl_min, lhfl_std=lhfl_std)
'''

'''
loaded_data1 = np.load('shfl_lhfl.npz')
shfl_mean1 = loaded_data1['shfl_mean']
lhfl_mean1 = loaded_data1['lhfl_mean']
ymdhs1 = ['2023010600','2023010700','2023010800','2023010900','2023011000','2023030700','2023030800','2023030900','2023031000','2023031100','2024021600','2024021700','2024021800','2024021900','2024022000','2025020100','2025020200','2025020300','2025020400','2025020500','2025022000','2025022100','2025022200','2025022300','2025022400','2026010300','2026010400','2026010500','2026010600','2026010700']

loaded_data2 = np.load('shfl_lhfl_2.npz')
shfl_mean2 = loaded_data2['shfl_mean']
lhfl_mean2 = loaded_data2['lhfl_mean']
ymdhs2 = ['2023010500','2023030600','2024021500','2025013100','2025021900','2025120400', '2025120500','2025120600','2025120700','2025120800','2025120900','2025121000','2025121100','2025121500','2025121600','2025121700','2025121800','2025121900','2025122000','2025122100','2025122200','2025122300','2025122400']

shfl_mean0 = np.hstack([shfl_mean1,shfl_mean2]) 
lhfl_mean0 = np.hstack([lhfl_mean1,lhfl_mean2]) 
ymdhs0 = np.hstack([ymdhs1,ymdhs2])

shfl_mean = shfl_mean0 * np.nan
lhfl_mean = lhfl_mean0 * np.nan
ymdhs_mean = []
for c,ymdh in enumerate(conf['ymdhs']):
    print(ymdh)
    okc = np.where(ymdhs0 == ymdh)[0][0]
    shfl_mean[:,c,:] = shfl_mean0[:,okc,:] 
    lhfl_mean[:,c,:] = lhfl_mean0[:,okc,:] 
    ymdhs_mean.append(ymdhs0[okc])

np.savez('shfl_lhfl_53_cycles.npz',shfl_mean=shfl_mean,lhfl_mean=lhfl_mean,ymdhs_mean=ymdhs_mean)
'''

loaded_data = np.load('shfl_lhfl_53_cycles.npz')
shfl_mean = loaded_data['shfl_mean']
lhfl_mean = loaded_data['lhfl_mean']
ymdhs_mean = loaded_data['ymdhs_mean']

########
#c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]
c = 0

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
    plt.plot(ff,shfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10)
for c in [0,6,12,19,24,30,38,48]:
    plt.plot(ff,shfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdhs'][c][0:-4],markeredgecolor='k',markersize=10)

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
    plt.plot(ff,lhfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10)
for c in [0,6,12,19,24,30,38,48]:
    plt.plot(ff,lhfl_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c][0:-4])
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


