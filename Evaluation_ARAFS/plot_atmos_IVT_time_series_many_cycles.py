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
def compute_ivt(grb, nlat, nlon, g=9.80665):

    Ps_pa = np.array([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 925, 950, 975, 1000])

    Qs = np.empty((len(Ps_pa),nlat,nlon))
    Qs[:] = np.nan
    Us = np.empty((len(Ps_pa),nlat,nlon))
    Us[:] = np.nan
    Vs = np.empty((len(Ps_pa),nlat,nlon))
    Vs[:] = np.nan

    for lv,Ps in enumerate(Ps_pa):
        print(lv,' ',Ps)
        Qs[lv,:,:] = grb.select(fullName="Specific Humidity",level=str(Ps)+' mb')[0].data * 100
        Us[lv,:,:] = grb.select(fullName="U-Component of Wind",level=str(Ps)+' mb')[0].data
        Vs[lv,:,:] = grb.select(fullName="V-Component of Wind",level=str(Ps)+' mb')[0].data

    # trapezoidal layer means
    dp   = np.diff(Ps_pa)                                    # (L-1,), all > 0
    qbar = 0.5 * (Qs[:-1] + Qs[1:])
    ubar = 0.5 * (Us[:-1] + Us[1:])
    vbar = 0.5 * (Vs[:-1] + Vs[1:])
    dp3  = dp[:, None, None]

    IVT_u = np.sum(qbar * ubar * dp3, axis=0) / g
    IVT_v = np.sum(qbar * vbar * dp3, axis=0) / g
    IVT  = np.hypot(IVT_u, IVT_v)

    return IVT_u, IVT_v, IVT

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)

xlim = conf['xlim']
ylim = conf['ylim']

#================================================================
# Cycle the different experiments

ff = np.arange(0,121,6)

'''
IVT_mean = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
IVT_mean[:] = np.nan
IVT_max = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
IVT_max[:] = np.nan
IVT_min = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
IVT_min[:] = np.nan
IVT_std = np.empty((len(conf['COMmodels']),len(conf['ymdhs']),len(ff)))
IVT_std[:] = np.nan

for exp in np.arange(len(conf['COMmodels'][0:2])): 
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
            #print(f'grib2file: {grib2file}')
            grb = grib2io.open(grib2file,mode='r')
        
            lat = grb.select(shortName='NLAT')[0].data
            lon = grb.select(shortName='ELON')[0].data
         
            nlat = lat.shape[0]
            nlon = lat.shape[1]
    
            print('Calculating IVT')
            _, _, IVT = compute_ivt(grb,nlat,nlon)
    
            IVT_mean[exp,c,f] = np.nanmean(IVT)
            IVT_max[exp,c,f] = np.nanmax(IVT)
            IVT_min[exp,c,f] = np.nanmin(IVT)
            IVT_std[exp,c,f] = np.nanstd(IVT)

            end = time.perf_counter()
            print(f"Elapsed time: {end - start:.6f} seconds")

np.savez('IVT.npz', IVT_mean=IVT_mean, IVT_max=IVT_max, IVT_min=IVT_min, IVT_std=IVT_std)
'''

loaded_data = np.load('IVT.npz')
IVT_mean = loaded_data['IVT_mean']
IVT_max = loaded_data['IVT_max']
IVT_min = loaded_data['IVT_min']
IVT_std = loaded_data['IVT_std']

########
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]
fig,ax = plt.subplots(figsize = (13,7))
for exp in np.arange(len(conf['COMmodels'])):
    plt.plot(ff,IVT_mean[exp,c,:],'o-',color=conf['exp_colors'][exp],label=conf['exp_labels'][exp],markeredgecolor='k',markersize=10)
ax.legend()
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('IVT ($Kg m^{-1} s^{-1}$)',fontsize=18,labelpad=10)
plt.title('Mean IVT ' + conf['ymdhs'][c],fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,120,12))

##############
IVT_diff = np.diff(IVT_mean,axis=0)[0,:,:]

fig,ax = plt.subplots(figsize = (13,7))
for c,ymdh in enumerate(conf['ymdhs']):
    plt.plot(ff,IVT_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('IVT Diff($Kg m^{-1} s^{-1}$)',fontsize=18,labelpad=10)
plt.title('Mean IVT Coupled - Uncoupled ',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.savefig('IVT_couple_minus_uncoupled.png')

##############
c = np.where([ymdh == conf['ymdh'] for ymdh in conf['ymdhs']])[0][0]

fig,ax = plt.subplots(figsize = (13,7))
plt.plot(ff,IVT_diff[c,:],'o-',color=conf['ymdhs_colors'][c],markeredgecolor='k',markersize=10,label=conf['ymdhs'][c])
ax.legend(bbox_to_anchor=(1.12, 1.1), loc='upper right')
ax.plot(ff,np.zeros(len(ff)),'-k',linewidth=2)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=18,labelpad=10)
plt.ylabel('IVT Diff($Kg m^{-1} s^{-1}$)',fontsize=18,labelpad=10)
plt.title('Mean IVT Coupled - Uncoupled ',fontsize=18)
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))


