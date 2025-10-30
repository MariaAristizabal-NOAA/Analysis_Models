#!/usr/bin/env python3

"""This script plots the storm track and the best track for variuos experiments"""
import sys
import yaml
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#================================================================
def get_storm_track_and_int(file_track,storm_num):

    ff = open(file_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    lead_time = []
    intt = []
    rmww = []
    for l in f:
        if l.split(',')[1].strip() == storm_num:
            lat = float(l.split(',')[6][0:4])/10
            if l.split(',')[6][4] == 'N':
                lat = lat
            else:
                lat = -lat
            lon = float(l.split(',')[7][0:5])/10
            if l.split(',')[7][4] == 'E':
                lon = lon
            else:
                lon = -lon
            latt.append(lat)
            lont.append(lon)
            lead_time.append(int(l.split(',')[5][1:4]))
            intt.append(float(l.split(',')[8]))
            rmww.append(float(l.split(',')[19]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    intt = np.asarray(intt)
    rmww = np.asarray(rmww)
    lead_time, ind = np.unique(lead_time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    rmw_track = rmww[ind]

    return lon_track, lat_track, lead_time, int_track, rmw_track
#================================================================
def get_best_track_and_int(file_best_track):

    ff = open(file_best_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    time = []
    intt = []
    for l in f:
        lat = float(l.split(',')[6][0:4])/10
        if l.split(',')[6][-1] == 'N':
            lat = lat
        else:
            lat = -lat
        lon = float(l.split(',')[7][0:5])/10
        if l.split(',')[7][-1] == 'E':
            lon = lon
        else:
            lon = -lon
        latt.append(lat)
        lont.append(lon)
        time.append(str(int(l.split(',')[2])))
        intt.append(float(l.split(',')[8]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    intt = np.asarray(intt)
    time_track, ind = np.unique(time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    name_storm = l.split(',')[-11]

    return lon_track, lat_track, time_track, int_track, name_storm

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos2.yml', 'rt') as f:
    conf = yaml.safe_load(f)
storm_id = conf['stormID']
storm_num = conf['stormID'][0:2]
inittime = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
cycle = conf['ymdh']

exp_labels = conf['exp_labels']
exp_colors = conf['exp_colors']
best_track_file = conf['best_track_file']

###################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

time_cycle = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

###################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

###################################################################
#%% Loop the experiments
ff = np.arange(0,127,3)
lon_forec_track = np.empty((len(conf['COMhafs']),len(ff)))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(conf['COMhafs']),len(ff)))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(conf['COMhafs']),len(ff)))
lead_time[:] = np.nan
int_track = np.empty((len(conf['COMhafs']),len(ff)))
int_track[:] = np.nan
for i,folder in enumerate(conf['COMhafs']):
    print(folder)
    #%% Get storm track from trak atcf files
    if exp_labels[i].split('_')[0] == 'HFSA':
        file_track = folder + '/' + storm_id.lower() +'.' + cycle + '.hfsa.trak.atcfunix'
        print(file_track)
    if exp_labels[i].split('_')[0] == 'HFSB':
        file_track = folder + '/' + storm_id.lower() +'.' + cycle + '.hfsb.trak.atcfunix'
        print(file_track)

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

###################################################################
#%% Figure track
lon_offset=0
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

fig = plt.figure(figsize=(8, 4))
ax = plt.axes(projection=myproj)
ax.axis('scaled')
okt = np.logical_and(time_best_track >= time_cycle[0],time_best_track <= time_cycle[-1])

for i in np.arange(len(conf['COMhafs'])): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.legend(loc='upper right',bbox_to_anchor=[1.2,1.0])
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.min(lon_best_track[okt])-5,np.max(lon_best_track[okt])+5])
plt.ylim([np.min(lat_best_track[okt])-5,np.max(lat_best_track[okt])+5])

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

#plt.savefig('track'+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

########################################################################
#%% Figure intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)
tt_cycle = np.asarray([t.replace(tzinfo=None) for t in time_cycle])
oktt = np.logical_and(tt_cycle[::2] >= time_best_track[0],tt_cycle[::2] <= time_best_track[-1])

fig,ax1 = plt.subplots(figsize=(7, 4))
for ii,i in enumerate(np.arange(len(conf['COMhafs']))): 
    plt.plot(lead_time[i,::2],int_track[i,::2],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)

plt.plot(lead_time[i,::2][oktt],int_best_track[okt],'o-k',label='Best')
#plt.legend(loc='upper right',bbox_to_anchor=[1.1,1.15])
#plt.legend(loc='lower right')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylim([10,190])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Cycle '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([10,190])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)
#plt.savefig('intensity'+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
