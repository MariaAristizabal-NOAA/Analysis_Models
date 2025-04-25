#%% User input

url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2020/NDBC_buoys_Pamlico_Sound.nc'
#url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2020/41025h2020.nc'

exp_names = ['GFSv17']
exp_labels = ['GFSv17']
exp_colors = ['darkviolet']

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

scratch_folder = '/scratch1/NCEPDEV/climate/Jiande.Wang/working/scratch/HR5-download-HPSS/OCEAN/daily-mean/'

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

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
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_var_from_model_following_trajectory

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
'''
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder1 + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')
'''

################################################################################
#%% Time window
'''
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

################################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]
'''
#################################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#################################################################################
#%% Read GFS output
GFS_file = scratch_folder + 'HR5.20200111.daily.nc'
GFS = xr.open_dataset(GFS_file)

time = np.asarray(GFS.time)
lon = np.asarray(GFS.xh)
lat = np.asarray(GFS.yh)
sst_gfs = np.asarray(GFS.SST)

#################################################################################
#%% Read NDBC data
url = url_NDBC

gdata = xr.open_dataset(url)#,decode_times=False)

station_buoy = np.asarray(gdata.station)
lat_buoy = np.asarray(gdata.latitude)
lon_buoy = np.asarray(gdata.longitude)
time_buoy = np.asarray(gdata.time)
timestamp_buoy = mdates.date2num(time_buoy)

#sst_buoy = np.asarray(gdata.sea_surface_temperature[:,0,0])
sst_buoy = np.asarray(gdata.wtmp)

#times = np.asarray(gdata.time)
#timestamp = mdates.date2num(time)
#times = np.asarray(mdates.num2date(timestamps))

oktimeg = np.logical_and(time_buoy >= time[0],time_buoy <= time[-1])

# Fields within time window
stationB = station_buoy[oktimeg]
timeB = time_buoy[oktimeg]
timestampB = timestamp_buoy[oktimeg]
sstB = sst_buoy[oktimeg]
lonB = lon_buoy[oktimeg]
latB = lat_buoy[oktimeg]

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
plt.plot(lonB[0],latB[0],'*',color='orange',label='NDBC Buoys ')
plt.plot(lonB,latB,'*',color='orange')
plt.axis('scaled')
plt.xlim([-80,-70])
plt.ylim([30,40])
plt.legend()

#################################################################################

file_model = GFS_file
var_name = 'SST'
time_name='time'
lon_name='xh'
lat_name='yh'

if file_model.split('.')[-1] == 'nc':
    model = xr.open_dataset(file_model)
    time_model = np.asarray(model[time_name][:])
    timestamp_model = mdates.date2num(time_model)
    lon_model = np.asarray(model[lon_name][:])
    lat_model = np.asarray(model[lat_name][:])
        
    if model[var_name].ndim == 4:
        var = np.asarray(model[var_name][:,depth_level,:,:])
    if model[var_name].ndim == 3:
        var = np.asarray(model[var_name][:,:,:])

    for station in np.unique(stationB):
        print(station)
        okstat = stationB == station
        lon_obs = lonB[okstat]
        lat_obs = latB[okstat]
        time_obs = timeB[okstat] 
        timestamp_obs = timestampB[okstat]
        sst_obs = sstB[okstat]
                
        # Interpolating lat_obs and lon_obs into model grid
        sublon = np.interp(timestamp_model,timestamp_obs,lon_obs)
        sublat = np.interp(timestamp_model,timestamp_obs,lat_obs)
        
        oklon = np.round(np.interp(sublon,lon_model,np.arange(len(lon_model)))).astype('int')
        oklat = np.round(np.interp(sublat,lat_model,np.arange(len(lat_model)))).astype('int')
        
        if file_model.split('.')[-1] == 'nc':
            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name][:,depth_level,:,:])
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name][:,:,:])
        
        target_var = []
        for t in np.arange(len(oklon)):
            target_var.append(var[t,oklat[t],oklon[t]])

        fig,ax = plt.subplots(figsize=(8, 4))
        plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
        plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
        plt.plot(lon_obs,lat_obs,'*',color='orange',label='NDBC Buoy '+station)
        plt.axis('scaled')
        plt.xlim([-80,-70])
        plt.ylim([30,40])
        plt.legend()

        fig,ax = plt.subplots(figsize=(10, 4))
        plt.plot(time_obs,sst_obs,'.-',color='blue',label='NDBC Buoy '+station)
        plt.plot(time_model,target_var,'o-',color='orange',markeredgecolor='k',label='GFSv17',markersize=7)
        #plt.legend(loc='lower left')
        plt.legend()
        plt.title('Sea Surface Temperature '+GFS_file.split('/')[-1],fontsize=18)
        plt.ylabel('$^oC$',fontsize=14)
        date_form = DateFormatter("%m-%d")
        ax.xaxis.set_major_formatter(date_form)
       #plt.ylim([970,1020])


