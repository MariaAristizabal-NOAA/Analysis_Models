#%% User input

# forecasting cycle to be used
cycles = ['20240901','20240902','20240903','20240904','20240905','20240905','20240906','20240907','20240908','20240909','20240910','20240911','20240912','20240913','20240914','20240915','20240916','20240917','20240918','20240919','20240920','20240921','20240922','20240923','20240923','20240924','20240925','20240926','20240927','20240928','20240929','20240930','20241001','20241001','20241002','20241003','20241004','20241005','20241006','20241007','20241008','20241009','20241010','20241011','20241012','20241013','20241014','20241015']

lon_lim = [-100,-60]
lat_lim = [0,40]

scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
folder_rtofs1 = scratch_folder + 'RTOFS/'
folder_rtofs2 = scratch_folder + 'RTOFS_v2.5.test01/'

forec_hours = ['n024','f006','f012','f018']
prefix = 'rtofs_glo_3dz_' 

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# folder utils for Hycom
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################

import sys
import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

#sys.path.append(folder_myutils)
#from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
#                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
#                            Haversine

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Time window
'''
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]
tini = datetime.strptime(date_ini,'%Y/%m/%d')
tend = tini + timedelta(hours=192)
date_end = tend.strftime('%Y/%m/%d')
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
#%% Read rtofs

file1 = glob.glob(folder_rtofs1 + 'rtofs.' + cycles[0] + '/' + prefix + forec_hours[0] + '*.nc')[0]
RTOFS1 = xr.open_dataset(file1)
lat = np.asarray(RTOFS1['Latitude'][:])

time = []
timestamp = []
mean_sst1 = []
std_sst1 = []
mean_sst2 = []
std_sst2 = []
mean_bias = []
std_bias = []
bias_SST = np.empty((len(cycles),len(forec_hours),lat.shape[0],lat.shape[1]))
bias_SST[:] = np.nan
for c,cycle in enumerate(cycles):
    print(cycle)
    for f,fh in enumerate(forec_hours):
        print(fh)
        file1 = glob.glob(folder_rtofs1 + 'rtofs.' + cycle + '/' + prefix + fh + '*.nc')[0]
        file2 = glob.glob(folder_rtofs2 + 'rtofs.' + cycle + '/' + prefix + fh + '*.nc')[0]
        RTOFS1 = xr.open_dataset(file1)
        RTOFS2 = xr.open_dataset(file2)
        t = RTOFS1['MT'][:]
        timest= mdates.date2num(t)[0]
        time.append(mdates.num2date(timest))
        timestamp.append(timest)

        lon = np.asarray(RTOFS1['Longitude'][:])
        lat = np.asarray(RTOFS1['Latitude'][:])
         
        SST1 = np.asarray(RTOFS1['temperature'][0,0,:,:])
        mean_sst1.append(np.nanmean(SST1))
        std_sst1.append(np.nanstd(SST1))
    
        SST2 = np.asarray(RTOFS2['temperature'][0,0,:,:])
        mean_sst2.append(np.nanmean(SST2))
        std_sst2.append(np.nanstd(SST2))

        bias_sst = SST2 - SST1
        bias_SST[c,f,:,:] = bias_sst
        mean_bias.append(np.nanmean(bias_sst))
        std_bias.append(np.nanstd(bias_sst))

        '''
        cflevels = np.arange(-2,2.1,0.2)
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lon,lat,bias_sst,levels=cflevels,cmap='RdBu_r',extend='both')
        plt.colorbar(extendrect=True)
        plt.axis('scaled')
        plt.ylim(lat_lim)
        plt.xlim(lon_lim)
        plt.title('SST Diff RTOFSv2.5 - RTOFSv2.4  '+ str(time[f])[0:13])

        plt.figure()
        plt.plot(SST1,SST2,'.',color='green')
        plt.plot(np.arange(-5,36),np.arange(-5,36),color='grey',linewidth=3)
        plt.xlabel('RTOFSv2.4')
        plt.ylabel('RTOFSv2.5')
        plt.title(str(time[f])[0:13])
        '''

Bias_SST = np.reshape(bias_SST,(48*4,1710,742))
time_mean_Bias_SST = np.nanmean(Bias_SST,axis=0)
Bias_SST_std = np.nanstd(Bias_SST,axis=0)

cflevels = np.arange(-1,1.1,0.2)
plt.figure()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lon,lat,time_mean_Bias_SST,levels=cflevels,cmap='RdBu_r',extend='both')
plt.colorbar(extendrect=True)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('SST Bias RTOFSv2.5 - RTOFSv2.4')

cflevels = np.arange(0,1.1,0.1)
plt.figure()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lon,lat,Bias_SST_std,levels=cflevels,cmap='RdBu_r',extend='max')
plt.colorbar(extendrect=True)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('STD: SST Bias RTOFSv2.5 - RTOFSv2.4')

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(time,mean_sst1,'.-',color='magenta',label='RTOFSv2.4')
plt.plot(time,mean_sst2,'.-',color='salmon',label='RTOFSv2.5')
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.grid(True)
plt.title('Domain Mean SST')
plt.ylabel('$^oC$',fontsize=14)

#ax.fill_between(time,np.asarray(mean_sst1)+np.asarray(std_sst1),np.asarray(mean_sst1)-np.asarray(std_sst1),color='blue',alpha=0.3)

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(time,np.zeros(len(time)),'-k')
plt.plot(time,mean_bias,'.-',color='cyan')
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.ylim([-0.2,0.2])
plt.grid(True)
plt.title('SST Bias RTOFSv2.5 - RTOFSv2.4')
plt.ylabel('$^oC$',fontsize=14)

