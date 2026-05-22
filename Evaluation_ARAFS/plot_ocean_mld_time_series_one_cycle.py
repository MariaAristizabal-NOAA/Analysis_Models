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
mld_mean = np.empty((len(conf['COMmodels']),len(ff)))
mld_mean[:] = np.nan
mld_max = np.empty((len(conf['COMmodels']),len(ff)))
mld_max[:] = np.nan
mld_min = np.empty((len(conf['COMmodels']),len(ff)))
mld_min[:] = np.nan
mld_std = np.empty((len(conf['COMmodels']),len(ff)))
mld_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    #===========================================================
    # Read fluxes
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
            mld = np.asarray(nc['MLD_0125'])
            mld_mean[exp,f] = np.nanmean(mld)
            mld_max[exp,f] = np.nanmax(mld)
            mld_min[exp,f] = np.nanmin(mld)
            mld_std[exp,f] = np.nanstd(mld)
        except FileNotFoundError:
            print("The file does not exist.")

########
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    if np.isfinite(np.mean(mld_mean[exp,:])):
        plt.plot(ff,mld_mean[exp,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdh'],markeredgecolor='k',markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('MLD ($m$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Mixed Layer Depth',fontsize=14)
plt.xlim([0,121])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (9,6))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,mld_mean[exp,:],'o-',color=conf['ymdhs_colors'][c],label=conf['ymdh'],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,mld_min[exp,:],mld_max[exp,:],color=conf['ymdhs_colors'][exp],alpha=0.1)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('MLD ($m$)',fontsize=14,labelpad=10)
plt.title('Domain Mean Mixed Layer Depth',fontsize=14)
plt.xlim([0,121])
plt.xticks(np.arange(0,120,12))


