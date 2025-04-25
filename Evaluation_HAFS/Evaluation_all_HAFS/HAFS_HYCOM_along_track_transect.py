#%% User input
# forecasting cycle to be used

# Danielle
cycle = '2022090306'
storm_num = '05'
basin = 'al'

#cycle = '2021082706'
#storm_num = '09'
#basin = 'al'

exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start']
exp_labels = ['HAFSv0.3a_fixed_wind_stresses']
exp_colors = ['c']

lon_lim = [-100,-40.0]
lat_lim = [0.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

##########################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob
import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, GOFS_coor_to_geo_coord 

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#########################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#########################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

########################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

########################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

temp_along_storm_track_hycom = np.empty((len(folder_exps),26,300))
temp_along_storm_track_hycom[:] = np.nan
temp_along_storm_track_hycom0 = np.empty((len(folder_exps),26,300))
temp_along_storm_track_hycom0[:] = np.nan
lon_forec_track_int_hycom = np.empty((len(folder_exps),300)) 
lon_forec_track_int_hycom[:] = np.nan
lat_forec_track_int_hycom = np.empty((len(folder_exps),300)) 
lat_forec_track_int_hycom[:] = np.nan
#%% Loop the experiments
for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

    #%% Read HAFS/HYCOM time
    time_hycom = []
    timestamp_hycom = []
    for n,file in enumerate(files_hafs_hycom):
        #print(file)
        HYCOM = xr.open_dataset(file)
        t = HYCOM['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time_hycom.append(mdates.num2date(timestamp))
        timestamp_hycom.append(timestamp)

    time_hycom = np.asarray(time_hycom)
    timestamp_hycom = np.asarray(timestamp_hycom)

    #%% HAFS-HYCOM temp along storm path
    lon_forec_track_interp = np.interp(lat_hafs_hycom,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_hafs_hycom)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_int_hycom[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_int_hycom[i,0:len(lat_forec_track_int)] = lat_forec_track_int
    lon_forec_track_int_hycom_g = lon_forec_track_int_hycom[i,:]
    lat_forec_track_int_hycom_g = lat_forec_track_int_hycom[i,:]
    oklon = np.round(np.interp(lon_forec_track_int,lon_hafs_hycom,np.arange(len(lon_hafs_hycom)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_hafs_hycom,np.arange(len(lat_hafs_hycom)))).astype(int)

    file = files_hafs_hycom[0]
    HYCOM = xr.open_dataset(file)

    for x in np.arange(len(lon_forec_track_int)):
        temp_along_storm_track_hycom0[i,:,x] = np.asarray(HYCOM['temperature'][0,:,oklat[x],oklon[x]])

    for n,file in enumerate(files_hafs_hycom):
        print(file)
        HYCOM = xr.open_dataset(file)
        fhour = file.split('/')[-1].split('.')[-2]

        for x in np.arange(len(lon_forec_track_int)):
            temp_along_storm_track_hycom[i,:,x] = np.asarray(HYCOM['temperature'][0,:,oklat[x],oklon[x]])

        #############################################################
        # Temp
        okl = np.isfinite(lon_forec_track_int_hycom_g)
        lat_forec_track = lat_forec_track_int_hycom_g[okl]
        temp_along_storm_track = temp_along_storm_track_hycom[i,:,okl]
        
        kw = dict(levels=np.arange(15,31.1,0.5))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(lat_forec_track,-depth_hafs_hycom,temp_along_storm_track.T,cmap='Spectral_r',**kw,extend='both')
        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(lat_forec_track,-depth_hafs_hycom,temp_along_storm_track.T,[26],colors='k')
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_title('HAFS: HYCOM Forecast for '+ storm_num + ' Init: ' + cycle + ' ' + fhour + '\n' + 'Temperature')
        fname = storm_num + '.' + cycle + '.hycom.temp_trans_along_track.' + fhour + '.png'
        fig.savefig(fname,bbox_inches='tight')
        plt.close()

        #############################################################
        # Temp differen
        okl = np.isfinite(lon_forec_track_int_hycom_g)
        lat_forec_track = lat_forec_track_int_hycom_g[okl]
        temp_diff = temp_along_storm_track_hycom[i,:,okl] - temp_along_storm_track_hycom0[i,:,okl]

        kw = dict(levels=np.arange(-4,4.1,0.2))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(lat_forec_track,-depth_hafs_hycom,temp_diff.T,cmap='seismic',**kw,extend='both')
        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(lat_forec_track,-depth_hafs_hycom,temp_diff.T,[0],colors='k')
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_xlabel('dT min = ' + str(np.round(np.nanmin(temp_diff),2)),fontsize=14) 
        ax.set_title('HAFS: HYCOM Forecast for '+ storm_num + ' Init: ' + cycle + ' ' + fhour + '\n' + 'Temperature Difference')
        fname = storm_num + '.' + cycle + '.hycom.temp_diff_trans_along_track.' + fhour + '.png'
        fig.savefig(fname,bbox_inches='tight')
        plt.close()

