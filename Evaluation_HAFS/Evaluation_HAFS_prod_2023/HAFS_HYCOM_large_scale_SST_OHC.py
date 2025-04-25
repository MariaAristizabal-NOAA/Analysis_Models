#%% User input

# Lee
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'

exp_names = ['HFSA_oper']
exp_labels = ['HFSA']
exp_colors = ['darkviolet']

xlim = [-70,-45.0]
ylim = [12.5,27.5]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_id + '/']

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = scratch_folder + 'adeck/a' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

#folder_fig = home_folder + 'Analysis/Evaluation_HAFS/Evaluation_HAFSv0p2a_phase3/Laura_2020/Figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import os.path
import glob
import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine 

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
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= xlim[0],bath_lat <= ylim[-1])
oklonbath = np.logical_and(bath_lon >= xlim[0],bath_lon <= ylim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read GFS track
#lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)

#################################################################################
#%% Loop the experiments
dtemp =  0.2
ref_depth = 10

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    file_track = glob.glob(os.path.join(folder,'*hfsa.trak.atcfunix'))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    #if len(files_hafs_fv3) == 0:
    #    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))
    
    #%% file to plot mlt and ohc
    f=0

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[f],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

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

    xlimh, ylimh  = geo_coord_to_HYCOM_coord(xlim,ylim)

    if np.min(lon_hafs_hycom) < 0:
        oklon = np.where(np.logical_and(lon_hafs_hycom>xlim[0],lon_hafs_hycom<xlim[1]))[0]
    else:
        oklon = np.where(np.logical_and(lon_hafs_hycom>xlimh[0],lon_hafs_hycom<xlimh[1]))[0]
    oklat = np.where(np.logical_and(lat_hafs_hycom>ylimh[0],lat_hafs_hycom<ylimh[1]))[0]

    hycom = xr.open_dataset(files_hafs_hycom[f])
    target_temp = np.asarray(hycom['temperature'][0,:,oklat,:][:,:,oklon])
    target_salt = np.asarray(hycom['salinity'][0,:,oklat,:][:,:,oklon])
    sst = target_temp[0,:,:]
    #target_temp = np.asarray(hycom['temperature'][0,:,oklat,:][:,:,oklon])
    #target_salt = np.asarray(hycom['salinity'][0,:,oklat,:][:,:,oklon])

    mlt_hycom = np.empty((len(oklat),len(oklon)))
    mlt_hycom[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            _, mlt_hycom[y,x] = MLD_temp_crit(dtemp,ref_depth,depth_hafs_hycom,target_temp[:,y,x])

    ohc_hycom = np.empty((len(oklat),len(oklon)))
    ohc_hycom[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],depth_hafs_hycom)
            ohc_hycom[y,x] = OHC_from_profile(depth_hafs_hycom,target_temp[:,y,x],dens_prof)
      
    if np.min(lon_hafs_hycom) < 0:
        lon_HAFS = lon_hafs_hycom[oklon]
    else:    
        lon_HAFS = lon_hafs_hycom[oklon] -360
    lat_HAFS = lat_hafs_hycom[oklat]

    kw = dict(levels=np.linspace(19, 32, 27))
    plt.figure(figsize=(10, 5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon_HAFS,lat_HAFS,sst,np.linspace(19, 32, 27),colors='grey',alpha=0.5)
    plt.contourf(lon_HAFS,lat_HAFS,sst,cmap='jet',**kw,extend='both')
    cbar = plt.colorbar(pad=0.008,extendrect=True)
    cbar.ax.set_ylabel('$^oC$',fontsize = 14)
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.axis('scaled')
    plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
    plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
    plt.title('SST' + exp_labels[i] + ' ' + str(time_fv3[f])[0:13])
    plt.legend(loc='lower left')

    kw = dict(levels=np.arange(0,151,10))
    plt.figure(figsize=(10, 5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver',extend='both')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon_HAFS,lat_HAFS,ohc_hycom,[50,100],colors='k')
    plt.contourf(lon_HAFS,lat_HAFS,ohc_hycom,cmap='jet',**kw)
    cbar = plt.colorbar(pad=0.005,extendrect=True)
    cbar.ax.set_ylabel('$Kj/cm^2$',fontsize = 14)
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.axis('scaled')
    plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
    plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
    plt.title('OHC '+ exp_labels[i] + ' ' + str(time_fv3[f])[0:13])
    plt.legend(loc='lower left')

#################################################################################
