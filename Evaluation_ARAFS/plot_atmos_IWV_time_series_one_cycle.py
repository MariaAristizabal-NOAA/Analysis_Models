#!/usr/bin/env python3

#exp_names = ['HFSAv2a_baseline_latest','HKPP','HKP2']
#exp_labels = ['HFSAv2a_baseline','HKPP','HKP2']
#exp_colors = ['orange','green','dodgerblue']

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
def compute_iwv(grb, nlat, nlon, g=9.80665):
    """
    IWV (a.k.a. precipitable water) = (1/g) ∫ q dp
    Sign-safe: integrates with ascending pressure so dp > 0.
    Returns IWV
    """

    Ps_pa = np.array([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925, 950, 975, 1000])

    Qs = np.empty((len(Ps_pa),nlat,nlon))
    Qs[:] = np.nan

    for lv,Ps in enumerate(Ps_pa):
        print(lv,' ',Ps)
        Qs[lv,:,:] = grb.select(fullName="Specific Humidity",level=str(Ps)+' mb')[0].data * 100

    # trapezoidal layer means
    dp   = np.diff(Ps_pa)                                    # (L-1,), all > 0
    qbar = 0.5 * (Qs[:-1] + Qs[1:])
    dp3  = dp[:, None, None]

    IWV_kgm2 = np.sum(qbar * dp3, axis=0) / g           # kg m^-2
    IWV_mm   = IWV_kgm2                                     # 1 kg m^-2 == 1 mm

    # Clean tiny negative numerical noise
    IWV_mm = np.where(IWV_mm < 0, np.clip(IWV_mm, 0, None), IWV_mm)
    IWV = IWV_mm

    return IWV

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
IWV_mean = np.empty((len(conf['COMmodels']),len(ff)))
IWV_mean[:] = np.nan
IWV_max = np.empty((len(conf['COMmodels']),len(ff)))
IWV_max[:] = np.nan
IWV_min = np.empty((len(conf['COMmodels']),len(ff)))
IWV_min[:] = np.nan
IWV_std = np.empty((len(conf['COMmodels']),len(ff)))
IWV_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'])): 

    #================================================================
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
        lat = grb.select(shortName='NLAT')[0].data
        lon = grb.select(shortName='ELON')[0].data
     
        nlat = lat.shape[0]
        nlon = lat.shape[1]

        print('Calculating IWV')
        IWV = compute_iwv(grb,nlat, nlon)

        IWV_mean[exp,f] = np.nanmean(IWV)
        IWV_max[exp,f] = np.nanmax(IWV)
        IWV_min[exp,f] = np.nanmin(IWV)
        IWV_std[exp,f] = np.nanstd(IWV)

########
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,IWV_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=10)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('IWV (mm)',fontsize=18,labelpad=10)
plt.title('Mean IWV ' + conf['ymdh'],fontsize=18)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,IWV_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,IWV_min[exp,:],IWV_max[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(ff,IWV_min[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
    plt.plot(ff,IWV_max[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('IWV (mm)',fontsize=14,labelpad=10)
plt.title('Mean, Max and Min IWV ' + conf['ymdh'],fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

##############
fig,ax = plt.subplots(figsize = (11,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,IWV_mean[exp,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=7)
    ax.fill_between(ff,IWV_mean[exp,:]-IWV_std[exp,:],IWV_mean[exp,:]+IWV_std[exp,:],color=conf['exp_colors'][exp],alpha=0.1)
    plt.plot(ff,IWV_mean[exp,:]-IWV_std[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
    plt.plot(ff,IWV_mean[exp,:]+IWV_std[exp,:],'--',color=conf['exp_colors'][exp],markersize=7)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylabel('IWV (mm)',fontsize=14,labelpad=10)
plt.title('Mean and Std IWV ' + conf['ymdh'],fontsize=14)
plt.xlim([0,120])
plt.xticks(np.arange(0,120,12))

