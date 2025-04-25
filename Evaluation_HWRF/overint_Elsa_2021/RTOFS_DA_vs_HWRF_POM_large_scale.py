#%% User input

# forecasting cycle to be used
cycle = '2021070200'
storm_id = '05l'

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

folder_hwrf = scratch_folder + 'HWRF2021_oper/' + storm_id + '/' + cycle
folder_rtofs_da = scratch_folder + 'RTOFS_DA/rtofs_da.' + cycle[0:-2] + '/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = scratch_folder + 'bdeck/bal052021.dat'
GFS_track_file = scratch_folder + 'adeck/aal052021.dat'

# RTOFS grid file name
#rtofs_grid = scratch_folder + 'RTOFS/' + 'GRID_DEPTH/regional.grid'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import sys
import os
import os.path
import glob
import sys
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
#%% Read GFS track
lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)
 
#################################################################################
#%% Get list files
#files_rtofs_da = sorted(glob.glob(os.path.join(folder_rtofs_da,'*z.f*archv.a.*')))
files_rtofs_da = sorted(glob.glob(os.path.join(folder_rtofs_da,'*_US_east.nc')))
files_hwrf_pom = sorted(glob.glob(os.path.join(folder_hwrf,'*pom*00*.nc')))
file_hwrf_pom_grid = sorted(glob.glob(os.path.join(folder_hwrf,'*pom*grid*.nc')))[0]
files_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*hafsprs*.nc')))

################################################################################
#%% Reading RTOFS grid
#lines_grid = [line.rstrip() for line in open(rtofs_grid+'.b')]
#lon_rtofs = np.array(readgrids(rtofs_grid,'plon:',[0]))
#lat_rtofs = np.array(readgrids(rtofs_grid,'plat:',[0]))

# Extracting the longitudinal and latitudinal size array
#idm=int([line.split() for line in lines_grid if 'longitudinal' in line][0][0])
#jdm=int([line.split() for line in lines_grid if 'latitudinal' in line][0][0])

rtofs_grid = xr.open_dataset(files_rtofs_da[0],decode_times=False)
lon_rtofs = np.asarray(rtofs_grid['Longitude'][:])
lat_rtofs = np.asarray(rtofs_grid['Latitude'][:])
depth_rtofs = np.asarray(rtofs_grid['Depth'][:])

################################################################################
#%% Reading HWRF/POM grid
print('Retrieving coordinates from POM')
pom_grid = xr.open_dataset(file_hwrf_pom_grid,decode_times=False)
lon_pom = np.asarray(pom_grid['east_e'][:])
lat_pom = np.asarray(pom_grid['north_e'][:])
zlev_pom = np.asarray(pom_grid['zz'][:])
hpom = np.asarray(pom_grid['h'][:])
#zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1))
#zmatrix_pom = zmatrix.reshape(hpom.shape[0],hpom.shape[1],zlev_pom.shape[0])
zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
zmatrix_pom = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])

#################################################################################
#%% Get storm track from trak atcf files

file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*')))[0]

lon_forec_track_hwrf, lat_forec_track_hwrf, lead_time_hwrf, int_track_hwrf = \
get_storm_track_and_int(file_track_hwrf)

#################################################################################
#%% Read HWRF time
time_hwrf = []
for n,file in enumerate(files_hwrf):
    print(file)
    hwrf = xr.open_dataset(file)
    t = hwrf.variables['time'][:]
    timestamp = mdates.date2num(t)[0]
    time_hwrf.append(mdates.num2date(timestamp))

time_hwrf = np.asarray(time_hwrf)

#################################################################################
#%% Read POM time
time_pom = []
timestamp_pom = []
for n,file in enumerate(files_hwrf_pom):
    print(file)
    pom = xr.open_dataset(file)
    t = pom['time'][:]
    timestamp = mdates.date2num(t)[0]
    time_pom.append(mdates.num2date(timestamp))
    timestamp_pom.append(timestamp)

time_pom = np.asarray(time_pom)
timestamp_pom = np.asarray(timestamp_pom)    

#################################################################################
#%% Read RTOFS time
time_rtofs = []
timestamp_rtofs = []
for n,file in enumerate(files_rtofs_da):
    print(file)
    RTOFS = xr.open_dataset(file)
    t = RTOFS['MT'][:]
    timestamp = mdates.date2num(t)[0]
    time_rtofs.append(mdates.num2date(timestamp))
    timestamp_rtofs.append(timestamp)

time_rtofs = np.asarray(time_rtofs)
timestamp_rtofs = np.asarray(timestamp_rtofs)
okt = np.argsort(timestamp_rtofs)

timestamp_rtofs = timestamp_rtofs[okt]
time_rtofs = time_rtofs[okt]
files_rtofs_da = np.asarray(files_rtofs_da)[okt]
   
################################################################################
#%% Figure SST all domain RTOFS-DA
okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

oktime_rtofs = np.where(np.logical_and(timestamp_rtofs >= timestamp_pom[0],\
                                       timestamp_rtofs <= timestamp_pom[-1]))[0]

n = 0 
rtofs0 = xr.open_dataset(np.asarray(files_rtofs_da)[oktime_rtofs][0])
SST0_rtofs = np.asarray(rtofs0['temperature'][0,0,:,:])
#kw = dict(levels=np.arange(24,32,0.5))
kw = dict(levels=np.arange(19,34,1))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
#plt.contourf(lon_rtofs,lat_rtofs,SST0_rtofs,cmap=cmocean.cm.thermal,**kw,vmin=24,vmax=32)
plt.contourf(lon_rtofs,lat_rtofs,SST0_rtofs,cmap='RdYlBu_r',**kw,vmin=19,vmax=34)
plt.plot(lon_best_track[okt], lat_best_track[okt],'.-',color='k',markeredgecolor='k',label='Best Track',markersize=7)
plt.plot(lon_best_track[okt][::2][n], lat_best_track[okt][::2][n],'o-',color='red',markeredgecolor='k',markersize=7)
plt.colorbar()
#plt.contourf(lon_rtofs,lat_rtofs,SST0_rtofs,[26,28,29,30,31],colors='k',alpha=0.3)
plt.contour(lon_rtofs,lat_rtofs,SST0_rtofs,[26],colors='r',alpha=1)
plt.axis('scaled')
#plt.ylim([0,45])
plt.ylim([0,40])
#plt.xlim([-100,np.max(lon_rtofs)])
plt.xlim([-100,-40])
plt.title('SST RTOFSv2 '+ str(time_rtofs[oktime_rtofs][::2][n])[0:13])
plt.text(-95,5,file.split('/')[-1].split('_')[3],fontsize=14)
plt.legend()
file_name = folder_fig + 'SST0_RTOFS_' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
'''
#%% Figure SST diff all domain RTOFS-DA
for n,file in enumerate(np.asarray(files_rtofs_da)[oktime_rtofs][::2]):
    print(n,file)
    HYCOM0 = xr.open_dataset(np.asarray(files_rtofs_da)[oktime_rtofs][0])
    HYCOM = xr.open_dataset(file)
    diff_SST_hafs_hycom = np.asarray(HYCOM['temperature'][0,0,:,:])-np.asarray(HYCOM0['temperature'][0,0,:,:])
    kw = dict(levels=np.arange(-3.0,3.1,0.1))
    plt.figure()
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    #plt.contour(lon_hafs_hycom-360,lat_hafs_hycom,diff_SST_hafs_hycom) #,[24,26,28,29,30,31],colors='grey',alpha=0.5)
    plt.contourf(lon_rtofs,lat_rtofs,diff_SST_hafs_hycom,cmap='seismic',**kw,vmin=-3.0,vmax=3.1)
    plt.plot(lon_best_track[okt], lat_best_track[okt],'.-',color='k',markeredgecolor='k',label='Best Track',markersize=7)
    plt.plot(lon_best_track[okt][::2][n], lat_best_track[okt][::2][n],'o-',color='red',markeredgecolor='k',markersize=7)
    plt.colorbar()
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,np.max(lon_rtofs)])
    plt.title('Diff SST RTOFS '+ str(time_rtofs[oktime_rtofs][::2][n])[0:13])
    plt.text(-95,5,file.split('/')[-1].split('_')[3],fontsize=14)
    plt.legend()
    file_name = folder_fig + 'diff_SST_RTOFS_' + file.split('/')[-1].split('.')[-2]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
'''
################################################################################
#%% Figure SST all domain POM
n=0
file = files_hwrf_pom[0]
pom0 = xr.open_dataset(files_hwrf_pom[0])
SST0_pom = np.asarray(pom0['t'][0,0,:,:])
#kw = dict(levels=np.arange(24,32,0.5))
kw = dict(levels=np.arange(19,34,1))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lon_pom,lat_pom,SST0_pom,cmap='RdYlBu_r',**kw,vmin=19,vmax=34)
plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'.-',color='green',markeredgecolor='k',label='HWRF',markersize=7)
plt.plot(lon_forec_track_hwrf[::4][n],lat_forec_track_hwrf[::4][n],'o-',color='red',markeredgecolor='k',markersize=7)
plt.colorbar()
#plt.contourf(lon_pom,lat_pom,SST0_pom,[26,28,29,30,31],colors='k',alpha=0.3)
plt.contourf(lon_pom,lat_pom,SST0_pom,[26,26.1],colors='r',alpha=1)
plt.axis('scaled')
#plt.ylim([0,45])
plt.ylim([0,40])
#plt.xlim([-100,np.max(lon_rtofs)])
plt.xlim([-100,-40])
plt.title('SST HWRF-POM '+ str(time_pom[::2][n])[0:13])
plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
plt.legend()
file_name = folder_fig + 'SST_hwrf_pom' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
'''
#%% Figure SST diff all domain POM
for n,file in enumerate(files_hwrf_pom[::2]):
    print(n,file)
    pom0 = xr.open_dataset(files_hwrf_pom[0])
    pom = xr.open_dataset(file)
    diff_SST_pom = np.asarray(pom['t'][0,0,:,:])-np.asarray(pom0['t'][0,0,:,:])
    kw = dict(levels=np.arange(-3.0,3.1,0.1))
    plt.figure()
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    #plt.contour(lon_hafs_hycom-360,lat_hafs_hycom,diff_SST_hafs_hycom) #,[24,26,28,29,30,31],colors='grey',alpha=0.5)
    plt.contourf(lon_pom,lat_pom,diff_SST_pom,cmap='seismic',**kw,vmin=-3.0,vmax=3.0)
    plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'.-',color='green',markeredgecolor='k',label='HWRF',markersize=7)
    plt.plot(lon_forec_track_hwrf[::4][n],lat_forec_track_hwrf[::4][n],'o-',color='red',markeredgecolor='k',markersize=7)
    plt.colorbar()
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,np.max(lon_rtofs)])
    plt.title('Diff SST HWRF-POM '+ str(time_pom[::2][n])[0:13])
    plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
    plt.legend()
    file_name = folder_fig + 'diff_SST_hwrf_pom' + file.split('/')[-1].split('.')[-2]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
'''
################################################################################

#%% Figure OHC all domain RTOFS
n=0
file = np.asarray(files_rtofs_da)[oktime_rtofs][0]
print(file)
xlim = [-100,-40]
ylim = [0,40]

oklon = np.where(np.logical_and(lon_rtofs[0,:]>xlim[0],lon_rtofs[0,:]<xlim[1]))[0]
oklat = np.where(np.logical_and(lat_rtofs[:,0]>ylim[0],lat_rtofs[:,0]<ylim[1]))[0]

lon_RTOFS = lon_rtofs[oklat,:][:,oklon]
lat_RTOFS = lat_rtofs[oklat,:][:,oklon]

RTOFS = xr.open_dataset(file)
target_temp = np.asarray(RTOFS['temperature'][0,:,oklat,:][:,:,oklon])
target_salt = np.asarray(RTOFS['salinity'][0,:,oklat,:][:,:,oklon])

OHC_rtofs = np.empty((len(oklat),len(oklon)))
OHC_rtofs[:] = np.nan
for x in np.arange(len(oklon)):
    #print(x)
    for y in np.arange(len(oklat)):
        dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],depth_rtofs)
        OHC_rtofs[y,x] = OHC_from_profile(depth_rtofs,target_temp[:,y,x],dens_prof)

#kw = dict(levels=np.arange(0,101,10))
kw = dict(levels=np.arange(0,160,10))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
#plt.contour(lon_RTOFS,lat_RTOFS,OHC_rtofs,[50,100],colors='k')
plt.contour(lon_RTOFS,lat_RTOFS,OHC_rtofs,[50,100],colors='k')
plt.contourf(lon_RTOFS,lat_RTOFS,OHC_rtofs,cmap='jet',**kw)
c=plt.colorbar()
c.set_label('$KJ/cm^2$',rotation=90, labelpad=10, fontsize=14)
#plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
plt.xlim([-100,-40])
plt.ylim([0,40])
plt.title('OHC RTOFSv2 '+ str(time_rtofs[oktime_rtofs[n]])[0:13])
plt.legend()

################################################################################
#%% Figure OHC all domain HWRF-POM
n=0
file = files_hwrf_pom[0]
print(file)
xlim = [-100,-40]
ylim = [0,40]

oklon = np.where(np.logical_and(lon_pom[0,:]>xlim[0],lon_pom[0,:]<xlim[1]))[0]
oklat = np.where(np.logical_and(lat_pom[:,0]>ylim[0],lat_pom[:,0]<ylim[1]))[0]

lon = lon_pom[oklat,:][:,oklon]
lat = lat_pom[oklat,:][:,oklon]

MODEL = xr.open_dataset(file)
target_temp = np.asarray(MODEL['t'][0,:,oklat,:][:,:,oklon])
target_salt = np.asarray(MODEL['s'][0,:,oklat,:][:,:,oklon])
target_depth = -zmatrix_pom[:,oklat,:][:,:,oklon]

OHC_pom = np.empty((len(oklat),len(oklon)))
OHC_pom[:] = np.nan
for x in np.arange(len(oklon)):
    #print(x)
    for y in np.arange(len(oklat)):
        dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],target_depth[:,y,x])
        OHC_pom[y,x] = OHC_from_profile(target_depth[:,y,x],target_temp[:,y,x],dens_prof)

#kw = dict(levels=np.arange(0,101,10))
kw = dict(levels=np.arange(0,160,10))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lon,lat,OHC_pom,cmap='jet',**kw)
c=plt.colorbar()
#plt.contour(lon,lat,OHC_pom,[60],colors='k')
plt.contour(lon,lat,OHC_pom,[50,100],colors='k')
c.set_label('$KJ/cm^2$',rotation=90, labelpad=10, fontsize=14)
plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
plt.xlim([-100,-40])
plt.ylim([0,40])
plt.title('OHC HWRF-POM '+ str(time_pom[n])[0:13])
plt.legend()

################################################################################
#%% Figure SSH all domain HWRF-POM
n=0
file = files_hwrf_pom[0]
print(file)
xlim = [-100,-40]
ylim = [0,40]

oklon = np.where(np.logical_and(lon_pom[0,:]>xlim[0],lon_pom[0,:]<xlim[1]))[0]
oklat = np.where(np.logical_and(lat_pom[:,0]>ylim[0],lat_pom[:,0]<ylim[1]))[0]

lon = lon_pom[oklat,:][:,oklon]
lat = lat_pom[oklat,:][:,oklon]

MODEL = xr.open_dataset(file)
target_ssh = np.asarray(MODEL['elb'][0,oklat,:][:,oklon])

kw = dict(levels=np.arange(-1.25,1.26,0.25))
#kw = dict(levels=np.arange(-100,101,5))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lon,lat,target_ssh,cmap='RdYlBu_r',**kw)
#plt.contourf(lon,lat,target_ssh*100,cmap='RdYlBu_r',**kw,vmax=30,vmin=-30)
c=plt.colorbar()
plt.contour(lon,lat,target_ssh,[0.75],colors='k')
c.set_label('$m$',rotation=90, labelpad=10, fontsize=14)
plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
plt.xlim([-100,-40])
plt.ylim([0,40])
plt.title('SSH HWRF-POM '+ str(time_pom[n])[0:13])
plt.legend()

