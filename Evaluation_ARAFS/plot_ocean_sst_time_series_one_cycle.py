#!/usr/bin/env python3

"""This script is to plot out HAFS SST time series"""

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
import xarray as xr

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_ocean.yml:')
with open('plot_ocean.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['stormNumber'] = conf['stormID'][0:2]
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,3)
sst_mean = np.empty((len(conf['COMmodels']),len(ff)))
sst_mean[:] = np.nan
sst_max = np.empty((len(conf['COMmodels']),len(ff)))
sst_max[:] = np.nan
sst_min = np.empty((len(conf['COMmodels']),len(ff)))
sst_min[:] = np.nan
sst_std = np.empty((len(conf['COMmodels']),len(ff)))
sst_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    #===========================================================
    for f,fh in enumerate(ff):

        if len(str(fh))==3:
            fhour = str(fh)
        if len(str(fh))==2:
            fhour = '0'+str(fh)
        if len(str(fh))==1:
            fhour = '00'+str(fh)

        fname = conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+'f'+fhour+'.nc'
        ncfile = os.path.join(conf['COMmodels'][exp]+conf['ymdh']+'/00E/', fname)
        try: 
            nc = xr.open_dataset(ncfile)
            sst = np.asarray(nc['SST'][0,:,:])
            sst_mean[exp,f] = np.nanmean(sst)
            sst_max[exp,f] = np.nanmax(sst)
            sst_min[exp,f] = np.nanmin(sst)
            sst_std[exp,f] = np.nanstd(sst)
        except FileNotFoundError:
            print("The file does not exist.")

########
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    if np.isfinite(np.mean(sst_mean[exp,:])):
        plt.plot(ff,sst_mean[exp,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdh'],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Sea Surface Temperature',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,sst_mean[exp,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdh'],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,sst_min[exp,:],sst_max[exp,:],color=conf['ymdhs_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('SST ($^oC$)',fontsize=14,labelpad=10)
plt.title('Sea Surface Temperature',fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))


