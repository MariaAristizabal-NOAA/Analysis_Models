#%% User input
# forecasting cycle to be used

cycle = '2022092718'
storm_num = '09'
basin = 'al'
file_oisst = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/OISST/' + cycle[0:-4] + '/oisst-avhrr-v02r01.20220918_preliminary.nc'

#cycle = '2021082706'
#storm_num = '09'
#basin = 'al'

exp_names = ['hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_rt']
exp_colors = ['c']

lon_lim = [-100,-40.0]
lat_lim = [0.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# RTOFS grid file name
#rtofs_grid = scratch_folder + 'RTOFS/' + 'GRID_DEPTH/regional.grid'

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
#%% Read OISST file
OISST = xr.open_dataset(file_oisst)
time_oisst = np.asarray(OISST['time'][:])[0]
latoisst =  np.asarray(OISST['lat'][:])
lonoisst =  np.asarray(OISST['lon'][:])
sstoisst =  np.asarray(OISST['sst'][0,0,:,:])

lonoisstt, lat_oisst = GOFS_coor_to_geo_coord(lonoisst,latoisst)
oklon = np.argsort(lonoisstt)
lon_oisst = lonoisstt[oklon]
sst_oisst = sstoisst[:,oklon]

#################################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

sst_along_storm_track_hycom = np.empty((len(folder_exps),300))
sst_along_storm_track_hycom[:] = np.nan
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

    #%% HAFS-HYCOM sst along storm path
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

    file = files_hafs_hycom[0]
    MODEL = xr.open_dataset(file)

    for x in np.arange(len(lon_forec_track_int)):
        sst_along_storm_track_hycom[i,x] = np.asarray(MODEL['temperature'][0,0,oklat[x],oklon[x]])

#################################################################################
#%% Figure SST all domain HAFS
n = 0
file = files_hafs_hycom[n]
HYCOM0 = xr.open_dataset(files_hafs_hycom[n])
SST_hafs_hycom = np.asarray(HYCOM0['temperature'][0,0,:,:])

if cycle == '2021070100':
    kw = dict(levels=np.arange(22,31,0.5))
    lev = [22,23,24,25,26,27,28,29,30,31]

if cycle == '2021082800' or '2021093006':
    kw = dict(levels=np.arange(22,32.5,0.5))
    lev = [22,23,24,25,26,27,28,29,30,31,32]

plt.figure(figsize=(8,5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contour(lon_hafs_hycom,lat_hafs_hycom,SST_hafs_hycom,lev,colors='grey',alpha=0.5)
plt.contourf(lon_hafs_hycom,lat_hafs_hycom,SST_hafs_hycom,cmap='Spectral_r',**kw) 
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
cbar = plt.colorbar()
cbar.ax.set_ylabel('$^oC$',fontsize = 14)
plt.axis('scaled')
plt.ylim([2,40])
plt.xlim([-98,-40])
plt.title('SST HAFS-HYCOM '+ str(time_hycom[n])[0:13])
plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
plt.legend()
file_name = folder_fig + 'SST_hafs' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% Figure SST all domain OISST
    
y = str(time_hycom[n].year)
m = [str(time_hycom[n].month) if len(str(time_hycom[n].month))>1 else '0'+str(time_hycom[n].month)][0]
d = [str(time_hycom[n].day) if len(str(time_hycom[n].day))>1 else '0'+str(time_hycom[n].day)][0]

OISST = xr.open_dataset(file_oisst)
time_oisst = np.asarray(OISST['time'][:])[0]
latoisst =  np.asarray(OISST['lat'][:])
lonoisst =  np.asarray(OISST['lon'][:])
sstoisst =  np.asarray(OISST['sst'][0,0,:,:])

lonoisstt, lat_oisst = GOFS_coor_to_geo_coord(lonoisst,latoisst)
oklon_oisstt = np.argsort(lonoisstt)
lon_oisst = lonoisstt[oklon_oisstt]
sst_oisst = sstoisst[:,oklon_oisstt]
oklon_oisst = np.where(np.logical_and(lon_oisst>lon_lim[0],lon_oisst<lon_lim[1]))[0]
oklat_oisst = np.where(np.logical_and(lat_oisst>lat_lim[0],lat_oisst<lat_lim[1]))[0]
lon_oisst = lon_oisst[oklon_oisst]
lat_oisst = lat_oisst[oklat_oisst]
sst_oisst = sst_oisst[oklat_oisst,:][:,oklon_oisst]

plt.figure(figsize=(8,5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contour(lon_oisst,lat_oisst,sst_oisst,lev,colors='grey',alpha=0.5)
plt.contourf(lon_oisst,lat_oisst,sst_oisst,cmap='Spectral_r',**kw) #,vmin=-3.0,vmax=3.0)
#plt.contourf(lon_oisst,lat_oisst,sst_oisst,cmap=cmocean.cm.thermal,**kw) #,vmin=-3.0,vmax=3.0)
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
cbar = plt.colorbar()
cbar.ax.set_ylabel('$^oC$',fontsize = 14)
plt.axis('scaled')
plt.ylim([2,40])
plt.xlim([-98,-40])
plt.title('SST OISST '+ y + '-' + m + '-' + d )
plt.legend()
file_name = folder_fig + 'SST_oisst' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% OISST sst along storm path

sst_along_storm_track_oisst = np.empty((len(exp_names),300))
sst_along_storm_track_oisst[:] = np.nan
lat_forec_track_int_oisst = np.empty((len(exp_names),300))
lat_forec_track_int_oisst[:] = np.nan
lon_forec_track_int_oisst = np.empty((len(exp_names),300))
lon_forec_track_int_oisst[:] = np.nan
for i in np.arange(len(exp_names)):
    lon_forec_track_interp = np.interp(lat_oisst,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_oisst)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_int_oisst[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_int_oisst[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    oklon = np.round(np.interp(lon_forec_track_int,lon_oisst,np.arange(len(lon_oisst)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_oisst,np.arange(len(lat_oisst)))).astype(int)

    for x in np.arange(len(lon_forec_track_int)):
        sst_along_storm_track_oisst[i,x] = sst_oisst[oklat[x],oklon[x]]

################################################################################
#%% time series sst along storm path

fig,ax = plt.subplots(figsize=(10, 5))
for i in np.arange(len(exp_names)):
    plt.plot(lat_forec_track_int_hycom[i,:],sst_along_storm_track_hycom[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    plt.plot(lat_forec_track_int_oisst[i,:],sst_along_storm_track_oisst[i,:],'o-',color='blue',label='OISST',markeredgecolor='k',markersize=7)
#'o-',color='c',label='hafsv0.3a',markeredgecolor='k',markersize=7)
plt.legend()
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.ylabel('($^oC$)',fontsize=14)
plt.title('SST Along Track',fontsize=14)

if cycle == '2021082800':
   plt.xlim([22,30])
   plt.ylim([28.5,31.7])

if cycle == '2021093006':
   plt.xlim([20,40])
   plt.ylim([24,29.5])


