#%% User input
# forecasting cycle to be used

# Fiona
#cycle = '2022091806'
#storm_num = '07'
#basin = 'al'
#file_ohc = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/OSPO_OHC/' + cycle[0:4] + '/ospo_OHC_north_atlant_Sep_2022.nc'

# Ian
#cycle = '2022092300'
cycle = '2022092806'
storm_num = '09'
basin = 'al'
file_ohc = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/OSPO_OHC/' + cycle[0:4] + '/ospo_OHC_north_atlant_Sep_2022.nc'

#cycle = '2021082706'
#storm_num = '09'
#basin = 'al'

#exp_names = ['hafsv0p3a_2022rt_natl']
#exp_labels = ['HAFSv0.3a_rt']
exp_names = ['hafsv0p3a_2022rt_cpl']
exp_labels = ['HAFSv0.3a']
exp_colors = ['c']

lon_lim = [-100,-40.0]
lat_lim = [0.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch1_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

bath_file = scratch1_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

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

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read OHC file
OHC = xr.open_dataset(file_ohc)

time_ospo = np.asarray(OHC['time'][:])
okt = np.where(mdates.date2num(time_ospo) == mdates.date2num(datetime(tini.year,tini.month,tini.day)))[0][0]

lat_ospo =  np.asarray(OHC['latitude'][:])
lon_ospo =  np.asarray(OHC['longitude'][:])
ohc_ospo =  np.asarray(OHC['ohc'][okt,:,:])

#################################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

ohc_along_storm_track_hycom = np.empty((len(folder_exps),300))
ohc_along_storm_track_hycom[:] = np.nan
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
        print(file)
        HYCOM = xr.open_dataset(file)
        t = HYCOM['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time_hycom.append(mdates.num2date(timestamp))
        timestamp_hycom.append(timestamp)

    time_hycom = np.asarray(time_hycom)
    timestamp_hycom = np.asarray(timestamp_hycom)

    file = files_hafs_hycom[0]
    MODEL = xr.open_dataset(file)
 
    lon_limh, lat_limh  = geo_coord_to_HYCOM_coord(lon_lim,lat_lim)

    #############  OHC large scale
    if np.min(lon_hafs_hycom) < 0:
        oklon = np.where(np.logical_and(lon_hafs_hycom>lon_lim[0],lon_hafs_hycom<lon_lim[1]))[0]
    else:
        oklon = np.where(np.logical_and(lon_hafs_hycom>lon_limh[0],lon_hafs_hycom<lon_limh[1]))[0]
    oklat = np.where(np.logical_and(lat_hafs_hycom>lat_limh[0],lat_hafs_hycom<lat_limh[1]))[0]

    target_temp = np.asarray(MODEL['temperature'][0,:,oklat,:][:,:,oklon])
    target_salt = np.asarray(MODEL['salinity'][0,:,oklat,:][:,:,oklon])
    lonn_hafs_hycom = lon_hafs_hycom[oklon]
    latt_hafs_hycom = lat_hafs_hycom[oklat]

    ohc_hafs_hycom = np.empty((len(oklat),len(oklon)))
    ohc_hafs_hycom[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],depth_hafs_hycom)
            ohc_hafs_hycom[y,x] = OHC_from_profile(depth_hafs_hycom,target_temp[:,y,x],dens_prof)

    #############  OHC along forecasted track
    #%% HAFS-HYCOM OHC along storm path
    lon_forec_track_interp = np.interp(lat_hafs_hycom,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_hafs_hycom)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_int_hycom[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_int_hycom[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    #lon_forec_track_int_hycom_g, lat_forec_track_int_hycom_g  = geo_coord_to_HYCOM_coord(lon_forec_track_int_hycom,lat_forec_track_int_hycom)
    lon_forec_track_int_hycom_g = lon_forec_track_int_hycom[i,:]
    lat_forec_track_int_hycom_g = lat_forec_track_int_hycom[i,:]

    oklon = np.round(np.interp(lon_forec_track_int,lon_hafs_hycom,np.arange(len(lon_hafs_hycom)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_hafs_hycom,np.arange(len(lat_hafs_hycom)))).astype(int)

    target_temp = np.asarray(MODEL['temperature'][0,:,oklat,:][:,:,oklon])
    target_salt = np.asarray(MODEL['salinity'][0,:,oklat,:][:,:,oklon])

    for x in np.arange(len(oklon)):
        dens_prof = sw.dens(target_salt[:,x,x],target_temp[:,x,x],depth_hafs_hycom)
        ohc_along_storm_track_hycom[i,x]  = OHC_from_profile(depth_hafs_hycom,target_temp[:,x,x],dens_prof)
    ##############################################

    #%% Figure OHC all domain HAFS
    lev = np.arange(-9000,9100,100)
    kw = dict(levels=np.arange(0,161,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lonn_hafs_hycom,latt_hafs_hycom,ohc_hafs_hycom,lev,colors='grey',alpha=0.5)
    plt.contour(lonn_hafs_hycom,latt_hafs_hycom,ohc_hafs_hycom,[100],colors='k')
    plt.contourf(lonn_hafs_hycom,latt_hafs_hycom,ohc_hafs_hycom,cmap='Spectral_r',extend='max',**kw)
    cbar = plt.colorbar(extendrect=True)
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.axis('scaled')
    plt.ylim([2,40])
    plt.xlim([-98,-40])
    plt.title('OHC HAFS-HYCOM '+ str(time_hycom[0])[0:13])
    plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
    plt.legend()
    file_name = folder_fig + 'OHC_hafs' + file.split('/')[-1].split('.')[-2]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% OSPO ohc along storm path

ohc_along_storm_track_ospo = np.empty((len(exp_names),300))
ohc_along_storm_track_ospo[:] = np.nan
lat_forec_track_int_ospo = np.empty((len(exp_names),300))
lat_forec_track_int_ospo[:] = np.nan
lon_forec_track_int_ospo = np.empty((len(exp_names),300))
lon_forec_track_int_ospo[:] = np.nan
for i in np.arange(len(exp_names)):
    lon_forec_track_interp = np.interp(lat_ospo,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_ospo)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_int_ospo[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_int_ospo[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    oklon = np.round(np.interp(lon_forec_track_int,lon_ospo,np.arange(len(lon_ospo)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_ospo,np.arange(len(lat_ospo)))).astype(int)

    for x in np.arange(len(lon_forec_track_int)):
        ohc_along_storm_track_ospo[i,x] = ohc_ospo[oklat[x],oklon[x]]

#################################################################################
#%% Figure OHC all domain OSPO

kw = dict(levels=np.arange(0,161,20))
plt.figure(figsize=(8,5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contour(lon_ospo,lat_ospo,ohc_ospo,lev,colors='grey',alpha=0.5)
plt.contour(lon_ospo,lat_ospo,ohc_ospo,[100],colors='k',alpha=0.5)
plt.contourf(lon_ospo,lat_ospo,ohc_ospo,cmap='Spectral_r',extend='max',**kw)
cbar = plt.colorbar(extendrect=True)
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('Ocean Heat Content OSPO '+ str(time_ospo[okt])[0:13])
plt.legend()
file_name = folder_fig + 'OHC_ospo' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% time series sst along storm path

fig,ax = plt.subplots(figsize=(10, 5))
for i in np.arange(len(exp_names)):
    plt.plot(lat_forec_track_int_hycom[i,:],ohc_along_storm_track_hycom[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    plt.plot(lat_forec_track_int_ospo[i,:],ohc_along_storm_track_ospo[i,:],'o-',color='k',label='OSPO OHC',markeredgecolor='k',markersize=7)
#'o-',color='c',label='hafsv0.3a',markeredgecolor='k',markersize=7)
plt.legend()
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.ylabel('($kJ/cm^2$)',fontsize=14)
plt.title('Ocean Heat Content Along Track',fontsize=14)


