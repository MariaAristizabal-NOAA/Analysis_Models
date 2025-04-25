#%% User input

# forecasting cycle to be used
cycle = '2022090212'
storm_num = '05'
basin = 'al'

exp_names = ['hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_rt']
exp_colors = ['c']

#exp_names = ['hafsv0p2a_2021rt_natl','HAFSv0p3_baseline','HAFSv0p3_HFAB_HF3A']
#exp_labels = ['HAFSv0.2a','H3BL','HF3A']
#exp_colors = ['limegreen','indianred','orange']

# Transect lon and lat limits
lon_lim = [-44,-44]
lat_lim = [34,42]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_HYCOM,\
                            figure_transect_temp, glider_data_vector_to_array,\
                            grid_glider_data


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

################################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

#################################################################################
#%% Loop the experiments

target_temp_10m = np.empty((len(folder_exps),22))
target_temp_10m[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Read HAFS/HYCOM time
    time_hycom = []
    timestamp_hycom = []
    for n,file in enumerate(files_hafs_hycom):
        print(file)
        HYCOM = xr.open_dataset(file)
        t = HYCOM['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time_hycom.append(mdates.num2date(timestamp))
        timestamp_hycom.append(timestamp)

    time_hycom = np.asarray(time_hycom)
    timestamp_hycom = np.asarray(timestamp_hycom)

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])
    
    #%% Longitudinal transect
    lon = lon_hafs_hycom
    lat = lat_hafs_hycom
    depth = depth_hafs_hycom

    xlim = lon_lim
    ylim = lat_lim

    xmin = int(np.round(np.interp(xlim[0],lon,np.arange(len(lon)))))
    xmax = int(np.round(np.interp(xlim[1],lon,np.arange(len(lon)))))
    ymin = int(np.round(np.interp(ylim[0],lat,np.arange(len(lat)))))
    ymax = int(np.round(np.interp(ylim[1],lat,np.arange(len(lat)))))
    if xmin == xmax:
        xmax = xmax + 1

    latt = lat[ymin:ymax]

    MODEL0 = xr.open_dataset(files_hafs_hycom[0])
    temp0 = temp = np.asarray(MODEL0['temperature'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
    temp0[temp0 == 0] = np.nan

    for n,file in enumerate(files_hafs_hycom):
        print(file)
        hycom = xr.open_dataset(file)
        ff = str(n*6)
        if len(ff)==1:
            fhour = '00' + ff
        if len(ff)==2:
            fhour = '0' + ff
        if len(ff)==3:
            fhour = ff

        MODEL = xr.open_dataset(file)
        temp = np.asarray(MODEL['temperature'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp[temp == 0] = np.nan

        print(np.nanmin(temp-temp0))

        #############################################################
        # Temp

        kw = dict(levels=np.arange(15,31.1,0.5))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(latt,-depth,temp,cmap='Spectral_r',**kw,extend='both')
        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(latt,-depth,temp,[26],colors='k')
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_title('HAFS: HYCOM Forecast for '+ storm_num + ' Init: ' + cycle + ' F' + fhour + '\n' + 'Temperature')
        fname = storm_num + '.' + cycle + '.hycom.temp_trans.f' + fhour + '.png'
        fig.savefig(fname,bbox_inches='tight')
        plt.close()

        #####################################################################
        # Figure temp-temp_f000

        kw = dict(levels=np.arange(-4,4.1,0.2))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(latt,-depth,temp-temp0,cmap='seismic',**kw,extend='both')
        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(latt,-depth,temp-temp0,[0],colors='k',alpha=0.3)
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_title('HAFS: HYCOM Forecast for '+ storm_num + ' Init: ' + cycle + ' F' + fhour + '\n' + 'Temperature difference')
        ax.set_xlabel('dt min = ' + str(np.round(np.nanmin(temp-temp0),2)),fontsize=14)
        fname = storm_num + '.' + cycle + '.hycom.temp_diff_trans.f' + fhour + '.png'
        fig.savefig(fname,bbox_inches='tight')
        plt.close()

    ############################################################################
