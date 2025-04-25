#%% User input

# forecasting cycle to be used
# Lee
cycle = '2023090706'
storm_id = '13l'
fhour = 'f006'
storm_name = 'lee'
basin = 'al'

exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1']
exp_labels = ['HFSA','HWRF','HFSAv1p1']
exp_colors = ['darkviolet','pink','steelblue']

# Transect lon and lat limits
lon_lim = [-59.8,-59.8]
lat_lim = [15.8,22.4]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_id + '/',scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_id + '/',scratch_folder + exp_names[2] + '/' + cycle + '/' + storm_id + '/']

best_track_file = abdeck_folder + 'btk/b' + basin + storm_id + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#################################################################################
#%% Loop the experiments

target_sst = np.empty((len(folder_exps),100))
target_sst[:] = np.nan
target_sst_diff = np.empty((len(folder_exps),100))
target_sst_diff[:] = np.nan
target_lat = np.empty((len(folder_exps),100))
target_lat[:] = np.nan
for i,folder in enumerate(folder_exps):
    print(folder)    

    if exp_labels[i] == 'HFSA':
        file0 = folder + storm_id + '.' + cycle + '.hfsa.hycom.3z.f000.nc'
        file = folder + storm_id + '.' + cycle + '.hfsa.hycom.3z.' + fhour + '.nc'
        #%% Read HAFS/HYCOM time
        MODEL0 = xr.open_dataset(file0)
        MODEL = xr.open_dataset(file)
        t = MODEL['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading HAFS/HYCOM grid)
        lon = np.asarray(MODEL['Longitude'][:])
        lat = np.asarray(MODEL['Latitude'][:])
        depth = np.asarray(MODEL['Z'][:])
    
        #%% Longitudinal transect
        lon = lon -360
        lat = lat
        depth = depth

        xlim = lon_lim
        ylim = lat_lim

        xmin = int(np.round(np.interp(xlim[0],lon,np.arange(len(lon)))))
        xmax = int(np.round(np.interp(xlim[1],lon,np.arange(len(lon)))))
        ymin = int(np.round(np.interp(ylim[0],lat,np.arange(len(lat)))))
        ymax = int(np.round(np.interp(ylim[1],lat,np.arange(len(lat)))))
        if xmin == xmax:
            xmax = xmax + 1

        latt = lat[ymin:ymax]
        target_lat[i,0:len(lat[ymin:ymax])] = latt

        temp0 = np.asarray(MODEL0['temperature'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp0[temp0 == 0] = np.nan
        temp = np.asarray(MODEL['temperature'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp[temp0 == 0] = np.nan
        diff = temp - temp0

        print(np.nanmin(diff))

        target_sst[i,0:len(latt)] = temp[0,:]
        target_sst_diff[i,0:len(latt)] = diff[0,:]

    if exp_labels[i] == 'HFSAv1p1':
        file003 = folder + storm_id + '.' + cycle + '.hfsa.mom6.f003.nc'
        file = folder + storm_id + '.' + cycle + '.hfsa.mom6.' + fhour + '.nc'
        #%% Read HAFS/HYCOM time
        MODEL0 = xr.open_dataset(file003)
        MODEL = xr.open_dataset(file)
        t = MODEL['time'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading HAFS/HYCOM grid
        lon = np.asarray(MODEL['xh'][:])
        lat = np.asarray(MODEL['yh'][:])
        depth = np.asarray(MODEL['z_l'][:])

        #%% Longitudinal transect
        lon = lon
        lat = lat
        depth = depth

        xlim = lon_lim
        ylim = lat_lim

        xmin = int(np.round(np.interp(xlim[0],lon,np.arange(len(lon)))))
        xmax = int(np.round(np.interp(xlim[1],lon,np.arange(len(lon)))))
        ymin = int(np.round(np.interp(ylim[0],lat,np.arange(len(lat)))))
        ymax = int(np.round(np.interp(ylim[1],lat,np.arange(len(lat)))))
        if xmin == xmax:
            xmax = xmax + 1

        latt = lat[ymin:ymax]
        target_lat[i,0:len(lat[ymin:ymax])] = latt

        temp0 = np.asarray(MODEL0['temp'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp0[temp0 == 0] = np.nan
        temp = np.asarray(MODEL['temp'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp[temp0 == 0] = np.nan
        diff = temp - temp0

        print(np.nanmin(diff))

        target_sst[i,0:len(latt)] = temp[0,:]
        target_sst_diff[i,0:len(latt)] = diff[0,:]

    if exp_labels[i] == 'HWRF':
        fnu = int(int(fhour[1:])/6)
        if fnu <= 10:
            fnum = '0' + str(fnu)
        else:
            fnum = str(fnu)  
        file0 = folder + storm_name + storm_id + '.' + cycle + '.pom.0000.nc'
        file = folder + storm_name + storm_id + '.' + cycle + '.pom.00' + fnum + '.nc'
        grid = folder + storm_name + storm_id + '.' + cycle + '.pom.grid.nc'
        #%% Read HAFS/HYCOM time
        MODEL0 = xr.open_dataset(file0)
        MODEL = xr.open_dataset(file)
        GRID = xr.open_dataset(grid)
        t = MODEL['time'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading HAFS/HYCOM grid
        lonp = np.asarray(GRID['east_e'][:])
        latp = np.asarray(GRID['north_e'][:])
        zlev_pom = np.asarray(GRID['zz'][:])
        hpom = np.asarray(GRID['h'][:])
        zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
        zmatrix = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])     
  
        #%% Longitudinal transect
        lon = lonp[0,:]
        lat = latp[:,0]
        depth = zmatrix

        xlim = lon_lim
        ylim = lat_lim

        xmin = int(np.round(np.interp(xlim[0],lon,np.arange(len(lon)))))
        xmax = int(np.round(np.interp(xlim[1],lon,np.arange(len(lon)))))
        ymin = int(np.round(np.interp(ylim[0],lat,np.arange(len(lat)))))
        ymax = int(np.round(np.interp(ylim[1],lat,np.arange(len(lat)))))
        if xmin == xmax:
            xmax = xmax + 1

        latt = lat[ymin:ymax]
        target_lat[i,0:len(lat[ymin:ymax])] = latt

        temp0 = np.asarray(MODEL0['t'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp0[temp0 == 0] = np.nan
        temp = np.asarray(MODEL['t'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp[temp0 == 0] = np.nan
        diff = temp - temp0

        print(np.nanmin(diff))

        target_sst[i,0:len(latt)] = temp[0,:]
        target_sst_diff[i,0:len(latt)] = diff[0,:]

fig,ax = plt.subplots(figsize=(10,5))
for i,label in enumerate(exp_labels):
    ax.plot(target_lat[i,:],target_sst[i,:],'o-',color = exp_colors[i],markeredgecolor='k',label=label)
ax.legend()
ax.set_ylabel('SST',fontsize=14)
ax.set_xlabel('Latitude',fontsize=14)


