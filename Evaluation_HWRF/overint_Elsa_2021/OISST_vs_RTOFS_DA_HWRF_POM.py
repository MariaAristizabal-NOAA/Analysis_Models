#%% User input
# forecasting cycle to be used
#cycle = '2021070200'
#cycle = '2021070212'
#cycle = '2021070218'
cycle = '2021070312'
storm_id = '05l'

lon_lim = [-98.5,-50.0]
lat_lim = [10.0,35.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

folder_hwrf = scratch_folder + 'HWRF2021_oper/' + storm_id + '/' + cycle
folder_rtofs_da = scratch_folder + 'RTOFS_DA/rtofs_da.' + cycle[0:-2] + '/'
folder_oisst = scratch_folder + 'Data/OISST/' + cycle[0:-4] + '/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = scratch_folder + 'bdeck/bal052021.dat'
GFS_track_file = scratch_folder + 'adeck/aal052021.dat'

# RTOFS grid file name
#rtofs_grid = scratch_folder + 'RTOFS/' + 'GRID_DEPTH/regional.grid'

# folder utils for Hycom
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
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
files_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*hwrfprs*.nc')))
files_oisst = sorted(glob.glob(os.path.join(folder_oisst,'*oisst*.nc')))[0]

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
zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
zmatrix_pom = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])

################################################################################
#%% Reading HWRF grid
print('Retrieving coordinates from POM')
hwrf = xr.open_dataset(files_hwrf[0],decode_times=False)
lon_hwrf = np.asarray(hwrf['longitude'][:])
lat_hwr = np.asarray(hwrf['latitude'][:])

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

#################################################################################
#%% Read OISST file
'''
OISST = xr.open_dataset(file_oisst)
time_oisst = np.asarray(OISST['time'][:])[0]
latoisst =  np.asarray(OISST['lat'][:])
lonoisst =  np.asarray(OISST['lon'][:])
sstoisst =  np.asarray(OISST['sst'][0,0,:,:])

lonoisstt, lat_oisst = GOFS_coor_to_geo_coord(lonoisst,latoisst)    
oklon = np.argsort(lonoisstt)
lon_oisst = lonoisstt[oklon]
sst_oisst = sstoisst[:,oklon]

plt.figure()
plt.contourf(lon_oisst,lat_oisst,sst_oisst)
plt.colorbar()
'''

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
#plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
#plt.legend(loc='lower right',bbox_to_anchor=(1.2,0))
plt.legend(loc='upper right')
plt.title('Track Forecast Laura cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([lon_lim[0],lon_lim[1]])
plt.ylim([lat_lim[0],lat_lim[1]])
#################################################################################
#%% Figure intensity
okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)
tt_hwrf = np.asarray([t.replace(tzinfo=None) for t in time_hwrf])
oktt = np.logical_and(tt_hwrf[::2] >= time_best_track[0],tt_hwrf[::2] <= time_best_track[-1])

fig,ax1 = plt.subplots(figsize=(10, 4))
plt.plot(lead_time_hwrf[::2][oktt],int_best_track[okt],'o-k',label='Best')
plt.plot(lead_time_hwrf[::2],int_track_hwrf[::2],'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
plt.plot(lead_time_GFS, int_GFS_track,'o-',color='blue',label='GFS')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
#ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.legend(loc='upper right',fontsize=14)
plt.ylim([10,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Cycle '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([10,145])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

#################################################################################
#%% Get MLT around 2 degrees of storm eye HWRF-POM

dtemp =  0.2
ref_depth = 10

MLT_pom_mean = []
MLT_pom_min = []
MLT_pom_max = []
SST_pom_mean = []
SST_pom_min = []
SST_pom_max = []
SST_oisst_mean = []
SST_oisst_min = []
SST_oisst_max = []
for n,file in enumerate(files_hwrf_pom):
    print(file)
    lon_forec_track = lon_forec_track_hwrf
    lat_forec_track = lat_forec_track_hwrf
    lon = lon_pom
    lat = lat_pom
    time = time_pom

    xlim = [lon_forec_track[::2][n]-2,lon_forec_track[::2][n]+2]
    ylim = [lat_forec_track[::2][n]-2,lat_forec_track[::2][n]+2]

    #xlimh, ylimh  = geo_coord_to_HYCOM_coord(xlim,ylim)

    oklon = np.where(np.logical_and(lon[0,:]>xlim[0],lon[0,:]<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat[:,0]>ylim[0],lat[:,0]<ylim[1]))[0]

    MODEL = xr.open_dataset(file)
    target_temp = np.asarray(MODEL['t'][0,:,oklat,:][:,:,oklon])
    target_sst = np.asarray(MODEL['t'][0,0,oklat,:][:,oklon])
    target_sst[target_sst == 0] = np.nan
    target_salt = np.asarray(MODEL['s'][0,:,oklat,:][:,:,oklon])
    target_depth = -zmatrix_pom[:,oklat,:][:,:,oklon]
    target_dens = np.asarray(pom['rho'][0,:,oklat,:][:,:,oklon]) * 1000 + 1000
    target_dens[target_dens==1000.0] = np.nan
    target_t = np.asarray(MODEL['time'][:])
    timestamp = mdates.date2num(target_t)[0]
    target_time = mdates.num2date(timestamp)

    MLT = np.empty((len(oklat),len(oklon)))
    MLT[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            if target_depth[:,y,x][-1] <= 1.0:
                MLT[y,x] = np.nan
            else:
               _, MLT[y,x] = MLD_temp_crit(dtemp,ref_depth,target_depth[:,y,x],target_temp[:,y,x])   
    MLT_pom_mean.append(np.nanmean(MLT))
    MLT_pom_min.append(np.nanmin(MLT))
    MLT_pom_max.append(np.nanmax(MLT))

    SST_pom_mean.append(np.nanmean(target_sst))
    SST_pom_min.append(np.nanmin(target_sst))
    SST_pom_max.append(np.nanmax(target_sst))

    y = str(target_time.year)
    m = [str(target_time.month) if len(str(target_time.month))>1 else '0'+str(target_time.month)][0]
    d = [str(target_time.day) if len(str(target_time.day))>1 else '0'+str(target_time.day)][0]

    file_oisst = sorted(glob.glob(os.path.join(folder_oisst,'*oisst*'+y+m+d+'*.nc')))[0]
    OISST = xr.open_dataset(file_oisst)
    time_oisst = np.asarray(OISST['time'][:])[0]
    latoisst =  np.asarray(OISST['lat'][:])
    lonoisst =  np.asarray(OISST['lon'][:])
    sstoisst =  np.asarray(OISST['sst'][0,0,:,:])

    lonoisstt, lat_oisst = GOFS_coor_to_geo_coord(lonoisst,latoisst)
    oklon_oisstt = np.argsort(lonoisstt)
    lon_oisst = lonoisstt[oklon_oisstt]
    sst_oisst = sstoisst[:,oklon_oisstt]
    oklon_oisst = np.where(np.logical_and(lon_oisst>xlim[0],lon_oisst<xlim[1]))[0]
    oklat_oisst = np.where(np.logical_and(lat_oisst>ylim[0],lat_oisst<ylim[1]))[0]
    target_lon_oisst = lon_oisst[oklon_oisst]
    target_lat_oisst = lat_oisst[oklat_oisst]
    target_sst_oisst = sst_oisst[oklat_oisst,:][:,oklon_oisst]
    SST_oisst_mean.append(np.nanmean(target_sst_oisst))
    SST_oisst_min.append(np.nanmin(target_sst_oisst))
    SST_oisst_max.append(np.nanmax(target_sst_oisst))
     
    if np.logical_or(n==6,n==14):
        lonn = lon[0,oklon]
        latt = lat[oklat,0]

        kw = dict(levels=np.arange(18,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn,latt,MLT,cmap=cmocean.cm.thermal,**kw,vmin=26)
        plt.colorbar()
        plt.contour(lonn,latt,MLT,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track[::2], lat_forec_track[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(latt),np.max(latt)])
        plt.xlim([np.min(lonn),np.max(lonn)])
        plt.title('Mixed Layer Temp. HWRF-POM '+ str(time[n])[0:13])        
        plt.legend()

        kw = dict(levels=np.arange(18,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn,latt,target_sst,cmap=cmocean.cm.thermal,**kw,vmin=26)
        plt.colorbar()
        plt.contour(lonn,latt,target_sst,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track[::2], lat_forec_track[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(latt),np.max(latt)])
        plt.xlim([np.min(lonn),np.max(lonn)])
        plt.title('SST HWRF-POM '+ str(time[n])[0:13])
        plt.legend()

        kw = dict(levels=np.arange(18,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(target_lon_oisst,target_lat_oisst,target_sst_oisst,cmap=cmocean.cm.thermal,**kw,vmin=26)
        plt.colorbar()
        plt.contour(target_lon_oisst,target_lat_oisst,target_sst_oisst,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track[::2], lat_forec_track[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(latt),np.max(latt)])
        plt.xlim([np.min(lonn),np.max(lonn)])
        plt.title('SST OISST '+ str(time_oisst)[0:13])
        plt.legend()

#################################################################################
#%% Get MLT around 2 degrees of storm eye RTOFS
dtemp =  0.2
ref_depth = 10

oktime_rtofs = np.where(np.logical_and(timestamp_rtofs >= timestamp_pom[0],\
                                       timestamp_rtofs <= timestamp_pom[-1]))[0] 

MLT_RTOFS_mean = []
MLT_RTOFS_min = []
MLT_RTOFS_max = [] 
SST_RTOFS_mean = []
SST_RTOFS_min = []
SST_RTOFS_max = []
for n,file in enumerate(np.asarray(files_rtofs_da)[oktime_rtofs]):
    print(file)
    lon_forec_track = lon_forec_track_hwrf
    lat_forec_track = lat_forec_track_hwrf
    lon = lon_rtofs
    lat = lat_rtofs
    time = time_pom

    xlim = [lon_forec_track[::2][n]-2,lon_forec_track[::2][n]+2]
    ylim = [lat_forec_track[::2][n]-2,lat_forec_track[::2][n]+2]

    #xlimh, ylimh  = geo_coord_to_HYCOM_coord(xlim,ylim)

    oklon = np.where(np.logical_and(lon[0,:]>xlim[0],lon[0,:]<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat[:,0]>ylim[0],lat[:,0]<ylim[1]))[0]

    RTOFS = xr.open_dataset(file)
    target_temp = np.asarray(RTOFS['temperature'][0,:,oklat,:][:,:,oklon])
    target_salt = np.asarray(RTOFS['salinity'][0,:,oklat,:][:,:,oklon])
    SST_rtofs = np.asarray(RTOFS['temperature'][0,0,oklat,:][:,oklon])
    MLT_rtofs = np.empty((len(oklat),len(oklon)))
    MLT_rtofs[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            _, MLT_rtofs[y,x] = MLD_temp_crit(dtemp,ref_depth,depth_rtofs,target_temp[:,y,x])
   
    MLT_RTOFS_mean.append(np.nanmean(MLT_rtofs))
    MLT_RTOFS_min.append(np.nanmin(MLT_rtofs))
    MLT_RTOFS_max.append(np.nanmax(MLT_rtofs))

    SST_RTOFS_mean.append(np.nanmean(SST_rtofs))
    SST_RTOFS_min.append(np.nanmin(SST_rtofs))
    SST_RTOFS_max.append(np.nanmax(SST_rtofs))
    
    if np.logical_or(n==6,n==14):
        lon_RTOFS = lon[oklat,:][:,oklon]
        lat_RTOFS = lat[oklat,:][:,oklon]

        kw = dict(levels=np.arange(18,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lon_RTOFS,lat_RTOFS,MLT_rtofs,cmap=cmocean.cm.thermal,**kw,vmin=26)
        plt.colorbar()
        plt.contour(lon_RTOFS,lat_RTOFS,MLT_rtofs,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(lat_RTOFS),np.max(lat_RTOFS)])
        plt.xlim([np.min(lon_RTOFS),np.max(lon_RTOFS)])
        plt.title('Mixed Layer Temp. RTOFSv2 '+ str(time_rtofs[oktime_rtofs[n]])[0:13])
        plt.legend()

        kw = dict(levels=np.arange(18,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lon_RTOFS,lat_RTOFS,SST_rtofs,cmap=cmocean.cm.thermal,**kw,vmin=26)
        plt.colorbar()
        plt.contour(lon_RTOFS,lat_RTOFS,SST_rtofs,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(lat_RTOFS),np.max(lat_RTOFS)])
        plt.xlim([np.min(lon_RTOFS),np.max(lon_RTOFS)])
        plt.title('SST RTOFSv2 '+ str(time_rtofs[oktime_rtofs[n]])[0:13])
        plt.legend()

#################################################################################
#%% Figure mean MLT storm-scale
fig,ax = plt.subplots(figsize = (8,5))

plt.plot(lead_time_hwrf[::2][0:-1],SST_pom_mean,'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],SST_pom_min,SST_pom_max,color='darkviolet',alpha=0.1)


plt.plot(lead_time_hwrf[::2][0:-1],SST_RTOFS_mean,'o-',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7)
plt.plot(lead_time_hwrf[::2][0:-1],SST_RTOFS_max,'--',color='gold',markeredgecolor='k',markersize=7,alpha=0.3)
plt.plot(lead_time_hwrf[::2][0:-1],SST_RTOFS_min,'--',color='gold',markeredgecolor='k',markersize=7,alpha=0.3)
ax.fill_between(lead_time_hwrf[::2][0:-1],SST_RTOFS_min,SST_RTOFS_max,color='gold',alpha=0.1)

plt.plot(lead_time_hwrf[::2][0:-1],SST_oisst_mean,'o-',color='forestgreen',label='OISST',markeredgecolor='k',markersize=7)
plt.plot(lead_time_hwrf[::2][0:-1],SST_oisst_max,'--',color='forestgreen',markeredgecolor='k',markersize=7)
plt.plot(lead_time_hwrf[::2][0:-1],SST_oisst_min,'--',color='forestgreen',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],SST_oisst_min,SST_oisst_max,color='forestgreen',alpha=0.1)

plt.legend()

plt.title('SST Cycle '+ cycle,fontsize=16)
plt.ylabel('$^oC$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

#################################################################################
#%% Figure MLT vs wtmp vs tmp surface from HRWF storm-scale
'''
fig,ax = plt.subplots(figsize = (8,5))

plt.plot(lead_time_hwrf[::2][0:-1],MLT_pom_mean,'o-',color='darkviolet',label='Mixed Layer Temp. POM',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],MLT_pom_min,MLT_pom_max,color='darkviolet',alpha=0.1)

plt.plot(lead_time_hwrf[::2][0:-1],SST_pom_mean,'o-',color='plum',label='SST POM',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],SST_pom_min,MLT_pom_max,color='plum',alpha=0.1)

plt.plot(lead_time_hwrf[::2],wtmp_mean[::2]-273.15,'o-',color='darkcyan',label='WTMP HWRF',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2],wtmp_min[::2]-273.15,wtmp_max[::2]-273.15,color='darkcyan',alpha=0.1)

plt.plot(lead_time_hwrf[::2],tmp_mean[::2]-273.15,'o-',color='yellowgreen',label='TMP HWRF',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2],tmp_min[::2]-273.15,tmp_max[::2]-273.15,color='yellowgreen',alpha=0.1)

plt.legend()

plt.title('Cycle '+ cycle,fontsize=16)
plt.ylabel('$^oC$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))
'''

###############################################################################
