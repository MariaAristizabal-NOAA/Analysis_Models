#%% User input
# forecasting cycle to be used
#cycle = '2021070200'
#cycle = '2021070212'
cycle = '2021070218'
storm_id = '05l'

lon_lim = [-98.5,-50.0]
lat_lim = [10.0,35.0]

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
files_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*hwrfprs*.nc')))

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
#%% Figure RTOFS temp transect along storm path
lon_forec_track_interp = np.interp(lat_rtofs[:,0],lat_forec_track_hwrf,lon_forec_track_hwrf,left=np.nan,right=np.nan)
lat_forec_track_interp = np.copy(lat_rtofs[:,0])
lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
lon_forec_track_int_rtofs = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
lat_forec_track_int_rtofs = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
oklon = np.round(np.interp(lon_forec_track_int_rtofs,lon_rtofs[0,:],np.arange(len(lon_rtofs[0,:])))).astype(int)
oklat = np.round(np.interp(lat_forec_track_int_rtofs,lat_rtofs[:,0],np.arange(len(lat_rtofs[:,0])))).astype(int)
okdepth = np.where(depth_rtofs <= 200)[0]
    
n = 3
file = files_rtofs_da[n]
MODEL = xr.open_dataset(file)

trans_temp_rtofs = np.empty((len(depth_rtofs[okdepth]),len(lon_forec_track_int_rtofs)))
trans_temp_rtofs[:] = np.nan
trans_salt_rtofs = np.empty((len(depth_rtofs[okdepth]),len(lon_forec_track_int_rtofs)))
trans_salt_rtofs[:] = np.nan
for x in np.arange(len(lon_forec_track_int_rtofs)):
    trans_temp_rtofs[:,x] = np.asarray(MODEL['temperature'][0,okdepth,oklat[x],oklon[x]])
    trans_salt_rtofs[:,x] = np.asarray(MODEL['salinity'][0,okdepth,oklat[x],oklon[x]])

kw = dict(levels = np.arange(14,31,1))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(lat_rtofs[oklat,0],-depth_rtofs[okdepth],trans_temp_rtofs,cmap=cmocean.cm.thermal,**kw)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.contour(lat_rtofs[oklat,0],-depth_rtofs[okdepth],trans_temp_rtofs,[26,27,28,29,30],color='k',alpha=0.3)
plt.contour(lat_rtofs[oklat,0],-depth_rtofs[okdepth],trans_temp_rtofs,[26],color='k')
cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.title('RTOFSv2 Temp. along Forecasted Storm Track \n on '+str(time_rtofs[n])[0:13],fontsize=16)

kw = dict(levels = np.arange(33.5,37.6,0.2))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(lat_rtofs[oklat,0],-depth_rtofs[okdepth],trans_salt_rtofs,cmap=cmocean.cm.haline,**kw)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.contour(lat_rtofs[oklat,0],-depth_rtofs[okdepth],trans_salt_rtofs,[34.5,35],colors='k')
cbar.ax.set_ylabel(' ',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.title('RTOFSv2 Salinity along Forecasted Storm Track \n on '+str(time_rtofs[n])[0:13],fontsize=16)

#################################################################################
#%% Figure POM temp transect along storm path
lon_forec_track_interp = np.interp(lat_pom[:,0],lat_forec_track_hwrf,lon_forec_track_hwrf,left=np.nan,right=np.nan)
lat_forec_track_interp = np.copy(lat_pom[:,0])
lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

lon_forec_track_int_pom = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
lat_forec_track_int_pom = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]

oklon = np.round(np.interp(lon_forec_track_int_pom,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
oklat = np.round(np.interp(lat_forec_track_int_pom,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)
#okdepth = np.where(depth_pom <= 200)[0]

n = 0
file = files_hwrf_pom[n]
MODEL = xr.open_dataset(file)

trans_temp_pom = np.empty((len(zmatrix),len(lon_forec_track_int_pom)))
trans_temp_pom[:] = np.nan
trans_salt_pom = np.empty((len(zmatrix),len(lon_forec_track_int_pom)))
trans_salt_pom[:] = np.nan
trans_depth_pom = np.empty((len(zmatrix),len(lon_forec_track_int_pom)))
trans_depth_pom[:] = np.nan
for x in np.arange(len(lon_forec_track_int_pom)):
    trans_temp_pom[:,x] = np.asarray(MODEL['t'][0,:,oklat[x],oklon[x]])
    trans_salt_pom[:,x] = np.asarray(MODEL['s'][0,:,oklat[x],oklon[x]])
    trans_depth_pom[:,x] = zmatrix_pom[:,oklat[x],:][:,oklon[x]]

trans_temp_pom[trans_temp_pom==0] = np.nan
trans_salt_pom[trans_salt_pom==0] = np.nan

lat_matrix = np.tile(lat_pom[oklat,0],(len(zmatrix),1))
kw = dict(levels = np.arange(14,31,1))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(lat_matrix,trans_depth_pom,trans_temp_pom,cmap=cmocean.cm.thermal,**kw)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.contour(lat_matrix,trans_depth_pom,trans_temp_pom,[26,27,28,29,30],color='k',alpha=0.3)
plt.contour(lat_matrix,trans_depth_pom,trans_temp_pom,[26],color='k')
cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.title('POM Temp. along Forecasted Storm Track \n on '+str(time_pom[n])[0:13],fontsize=16)
plt.ylim([-200,0])

lat_matrix = np.tile(lat_pom[oklat,0],(len(zmatrix),1))
kw = dict(levels = np.arange(33.5,37.6,0.2))
fig,ax = plt.subplots(figsize=(10, 5))
plt.contourf(lat_matrix,trans_depth_pom,trans_salt_pom,cmap=cmocean.cm.haline,**kw)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
plt.contour(lat_matrix,trans_depth_pom,trans_salt_pom,[34.5,35],colors='k')
cbar.ax.set_ylabel(' ',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.title('POM Salit. along Forecasted Storm Track \n on '+str(time_pom[n])[0:13],fontsize=16)
plt.ylim([-200,0])


#################################################################################
#%% Time series temp and salinity at 10 m 
okdr = np.where(depth_rtofs <= 5.0)[0]
temp_rtofs = np.nanmean(trans_temp_rtofs[okdr,:],axis=0)    
salt_rtofs = np.nanmean(trans_salt_rtofs[okdr,:],axis=0)    

temp_pom = np.empty(trans_depth_pom.shape[1])
temp_pom[:] = np.nan
salt_pom = np.empty(trans_depth_pom.shape[1])
salt_pom[:] = np.nan
for x in range(trans_depth_pom.shape[1]): 
    okdp = np.where(-trans_depth_pom[:,x] < 5.0)[0]
    temp_pom[x] = np.mean(trans_temp_pom[okdp,x],axis=0)
    salt_pom[x] = np.mean(trans_salt_pom[okdp,x],axis=0)


fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(lat_forec_track_int_rtofs,temp_rtofs,'o-',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7)
plt.plot(lat_forec_track_int_pom,temp_pom,'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
plt.legend()
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.ylabel('($^oC$)',fontsize=14)
plt.title('SST Along Track',fontsize=14)

fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(lat_forec_track_int_rtofs,salt_rtofs,'o-',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7)
plt.plot(lat_forec_track_int_pom,salt_pom,'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
plt.legend()
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.ylabel('($^oC$)',fontsize=14)
plt.title('SSS Along Track',fontsize=14)

'''
#################################################################################
#%% Get MLT around 2 degrees of storm eye HWRF-POM
dtemp =  0.2
ref_depth = 10

MLT_pom_mean = []
MLT_pom_min = []
MLT_pom_max = []
OHC_pom_mean = []
OHC_pom_min = []
OHC_pom_max = []
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
    target_salt = np.asarray(MODEL['s'][0,:,oklat,:][:,:,oklon])
    target_depth = -zmatrix_pom[:,oklat,:][:,:,oklon]
    target_dens = np.asarray(pom['rho'][0,:,oklat,:][:,:,oklon]) * 1000 + 1000
    target_dens[target_dens==1000.0] = np.nan
    
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

    OHC = np.empty((len(oklat),len(oklon)))
    OHC[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],target_depth[:,y,x])
            OHC[y,x] = OHC_from_profile(target_depth[:,y,x],target_temp[:,y,x],dens_prof)   
 
    OHC_pom_mean.append(np.nanmean(OHC))
    OHC_pom_min.append(np.nanmin(OHC))
    OHC_pom_max.append(np.nanmax(OHC))

    if np.logical_or(n==6,n==14):
        lonn = lon[0,oklon]
        latt = lat[oklat,0]

        kw = dict(levels=np.arange(26,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn,latt,MLT,cmap=cmocean.cm.thermal,**kw)
        plt.colorbar()
        plt.contour(lonn,latt,MLT,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track[::2], lat_forec_track[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(latt),np.max(latt)])
        plt.xlim([np.min(lonn),np.max(lonn)])
        plt.title('Mixed Layer Temp. HWRF-POM '+ str(time[n])[0:13])        
        plt.legend()

        kw = dict(levels=np.arange(0,110,10))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn,latt,OHC,cmap='magma',**kw)
        plt.colorbar()
        plt.contour(lonn,latt,OHC,[60],colors='k')
        plt.plot(lon_forec_track[::2], lat_forec_track[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWrf-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(latt),np.max(latt)])
        plt.xlim([np.min(lonn),np.max(lonn)])
        plt.title('OHC HWRF-POM '+ str(time[n])[0:13])        
        plt.legend(loc='upper right')

#################################################################################
#%% Get MLT around 2 degrees of storm eye RTOFS
dtemp =  0.2
ref_depth = 10

oktime_rtofs = np.where(np.logical_and(timestamp_rtofs >= timestamp_pom[0],\
                                       timestamp_rtofs <= timestamp_pom[-1]))[0] 

MLT_RTOFS_mean = []
MLT_RTOFS_min = []
MLT_RTOFS_max = [] 
OHC_RTOFS_mean = []
OHC_RTOFS_min = []
OHC_RTOFS_max = [] 
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
    OHC_rtofs = np.empty((len(oklat),len(oklon)))
    OHC_rtofs[:] = np.nan
    MLT_rtofs = np.empty((len(oklat),len(oklon)))
    MLT_rtofs[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            _, MLT_rtofs[y,x] = MLD_temp_crit(dtemp,ref_depth,depth_rtofs,target_temp[:,y,x])
   
    MLT_RTOFS_mean.append(np.nanmean(MLT_rtofs))
    MLT_RTOFS_min.append(np.nanmin(MLT_rtofs))
    MLT_RTOFS_max.append(np.nanmax(MLT_rtofs))
    
    OHC_rtofs = np.empty((len(oklat),len(oklon)))
    OHC_rtofs[:] = np.nan
    for x in np.arange(len(oklon)):
        #print(x)
        for y in np.arange(len(oklat)):
            dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],depth_rtofs)
            OHC_rtofs[y,x] = OHC_from_profile(depth_rtofs,target_temp[:,y,x],dens_prof)
   
    OHC_RTOFS_mean.append(np.nanmean(OHC_rtofs))
    OHC_RTOFS_min.append(np.nanmin(OHC_rtofs))
    OHC_RTOFS_max.append(np.nanmax(OHC_rtofs))
    
    if np.logical_or(n==6,n==14):
        lon_RTOFS = lon[oklat,:][:,oklon]
        lat_RTOFS = lat[oklat,:][:,oklon]

        kw = dict(levels=np.arange(26,31.1,0.2))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lon_RTOFS,lat_RTOFS,MLT_rtofs,cmap=cmocean.cm.thermal,**kw)
        plt.colorbar()
        plt.contour(lon_RTOFS,lat_RTOFS,MLT_rtofs,[26,27,28,29,30],colors='k',alpha=0.3)
        plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(lat_RTOFS),np.max(lat_RTOFS)])
        plt.xlim([np.min(lon_RTOFS),np.max(lon_RTOFS)])
        plt.title('Mixed Layer Temp. RTOFSv2 '+ str(time_rtofs[oktime_rtofs[n]])[0:13])
        plt.legend()
 
        kw = dict(levels=np.arange(0,110,10))
        plt.figure()
        plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lon_RTOFS,lat_RTOFS,OHC_rtofs,cmap='magma',**kw)
        plt.colorbar()
        plt.contour(lon_RTOFS,lat_RTOFS,OHC_rtofs,[60],colors='k')
        plt.plot(lon_forec_track_hwrf[::2], lat_forec_track_hwrf[::2],'o-',color='darkviolet',markeredgecolor='k',label='HWRF-POM',markersize=7)
        plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
        plt.plot(lon_GFS_track, lat_GFS_track,'o-',color='blue',label='GFS Track')
        plt.ylim([np.min(lat_RTOFS),np.max(lat_RTOFS)])
        plt.xlim([np.min(lon_RTOFS),np.max(lon_RTOFS)])
        plt.title('OHC RTOFSv2 '+ str(time_rtofs[oktime_rtofs[n]])[0:13])
        plt.legend(loc='upper right')
        #plt.savefig('OHC',bbox_inches = 'tight',pad_inches = 0.1)  

#################################################################################
#%% Figure mean MLT storm-scale

fig,ax = plt.subplots(figsize = (8,5))

plt.plot(lead_time_hwrf[::2][0:-1],MLT_pom_mean,'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],MLT_pom_min,MLT_pom_max,color='darkviolet',alpha=0.1)


plt.plot(lead_time_hwrf[::2][0:-1],MLT_RTOFS_mean,'o-',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7)
plt.legend()
plt.plot(lead_time_hwrf[::2][0:-1],MLT_RTOFS_max,'--',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7,alpha=0.3)
plt.plot(lead_time_hwrf[::2][0:-1],MLT_RTOFS_min,'--',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7,alpha=0.3)
ax.fill_between(lead_time_hwrf[::2][0:-1],MLT_RTOFS_min,MLT_RTOFS_max,color='gold',alpha=0.1)

plt.title('Mixed Layer Temperature Cycle '+ cycle,fontsize=16)
plt.ylabel('$^oC$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

#################################################################################
#%% Figure mean ohc storm-scale

fig,ax = plt.subplots(figsize = (8,5))

plt.plot(lead_time_hwrf[::2][0:-1],OHC_pom_mean,'o-',color='darkviolet',label='HWRF-POM',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hwrf[::2][0:-1],OHC_pom_min,OHC_pom_max,color='darkviolet',alpha=0.1)


plt.plot(lead_time_hwrf[::2][0:-1],OHC_RTOFS_mean,'o-',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7)
plt.legend()
plt.plot(lead_time_hwrf[::2][0:-1],OHC_RTOFS_max,'--',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7,alpha=0.3)
plt.plot(lead_time_hwrf[::2][0:-1],OHC_RTOFS_min,'--',color='gold',label='RTOFSv2',markeredgecolor='k',markersize=7,alpha=0.3)
ax.fill_between(lead_time_hwrf[::2][0:-1],OHC_RTOFS_min,OHC_RTOFS_max,color='gold',alpha=0.1)

plt.title('OHC Cycle '+ cycle,fontsize=16)
plt.ylabel('$kJ/cm^2$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

################################################################################
'''
