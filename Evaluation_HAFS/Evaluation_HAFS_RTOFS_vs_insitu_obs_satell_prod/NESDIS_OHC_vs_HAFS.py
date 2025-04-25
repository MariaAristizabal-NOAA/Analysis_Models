#%% User input
# forecasting cycle to be used

# Milton
cycle = '2024100606'
storm_num = '14'
storm_id = '14l'
basin = 'al'
file_ohc = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NESDIS_OHC/2024/' + '/ohc_na14QG3_2024_280.nc'

exp_names = ['HFSA_oper','HFSB_oper']
exp_labels = ['HFSA','HFSB']
exp_colors = ['purple','limegreen']
hafs = ['hfsa','hfsb']

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom
folder_myutils= '/home/Maria.Aristizabal/Utils/'

#================================================================
# Calculate ocean heat content
def ohc(temp,salt,zl):
    cp = 3985 #Heat capacity in J/(kg K)
    zl_array = np.reshape(np.tile(zl,(temp.shape[1]*temp.shape[2],1)).T,(temp.shape[0],temp.shape[1],temp.shape[2]))
    no26 = temp < 26
    temp[no26] = np.nan
    salt[no26] = np.nan
    density = dens(salt,temp,zl_array)
    rho0 = np.nanmean(density,axis=0)
    zl_array_fac = (zl_array[0:-1,:,:] + zl_array[1:,:,:])/2
    zero_array = np.zeros((1,temp.shape[1],temp.shape[2]))
    bott_array = np.ones((1,temp.shape[1],temp.shape[2]))* zl_array_fac[-1,0,0] + (zl[-1] - zl[-2])
    zl_array_face = np.vstack((zero_array,zl_array_fac,bott_array))
    dz_array = np.diff(zl_array_face,axis=0)
    ohc = np.abs(cp * rho0 * np.nansum((temp-26)*dz_array,axis=0)) * 10**(-7) # in kJ/cm^2

    return ohc

################################################################################
import sys
import os
import glob
import xarray as xr
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, GOFS_coor_to_geo_coord 

from eos80 import dens

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

################################################################################
'''
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')
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
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read OHC file
OHC = xr.open_dataset(file_ohc)

time_nesdis = np.asarray(OHC['time'][:])
#okt = np.where(mdates.date2num(time_nesdis) == mdates.date2num(datetime(tini.year,tini.month,tini.day)))[0][0]

lat_nesdis =  np.asarray(OHC['latitude'][:])
lon_nesdis =  np.asarray(OHC['longitude'][:])
#ohc_nesdis =  np.asarray(OHC['ohc'][okt,:,:])
ohc_nesdis =  np.asarray(OHC['ohc'][0,:,:])
ohc_nesdis[ohc_nesdis == 0] = np.nan

#################################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

ohc_forec_track_hafs = np.empty((len(folder_exps),300))
ohc_forec_track_hafs[:] = np.nan
ohc_forec_track_hafs_mean = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_mean[:] = np.nan
ohc_forec_track_hafs_max = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_max[:] = np.nan
ohc_forec_track_hafs_min = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_min[:] = np.nan
lon_forec_track_hafs = np.empty((len(folder_exps),300)) 
lon_forec_track_hafs[:] = np.nan
lat_forec_track_hafs = np.empty((len(folder_exps),300)) 
lat_forec_track_hafs[:] = np.nan
#%% Loop the experiments
for i,folder in enumerate(folder_exps):
    print(folder)

    oceanf = glob.glob(os.path.join(folder,'*f000.nc'))[0].split('/')[-1].split('.')
    ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

    #%% Get list files
    if ocean == 'mom6':
        file_hafs_ocean = folder + '/' + storm_id + '.' + cycle + '.' + hafs[i] + '.mom6.f000.nc'
        nc = xr.open_dataset(file_hafs_ocean)
        temp = np.asarray(nc['temp'][0,:,:])
        salt = np.asarray(nc['so'][0,:,:])
        lon = np.asarray(nc['xh'])
        lat = np.asarray(nc['yh'])
        zl = np.asarray(nc['z_l'])
        time = np.asarray(nc['time'])
        ohc_hafs = ohc(temp,salt,zl)

    if ocean == 'hycom':
        file_hafs_ocean = folder + '/' + storm_id + '.' + cycle + '.' + hafs[i] + '.hycom.3z.f000.nc'
        nc = xr.open_dataset(file_hafs_ocean)
        #ohc_hafs = np.asarray(nc['ocean_heat_content'][0,:,:])
        #ohc_hafs[ohc_hafs<-1000] = np.nan
        temp = np.asarray(nc['temperature'][0,:,:,:])
        salt = np.asarray(nc['salinity'][0,:,:,:])
        lon = np.asarray(nc['Longitude'])
        lat = np.asarray(nc['Latitude'])
        depth = np.asarray(nc['Z'][:])
        time = np.asarray(nc['MT'][:])
        ohc_hafs = ohc(temp,salt,depth)
        #timestamp = mdates.date2num(t)[0]
        #time_hycom.append(mdates.num2date(timestamp))
        #timestamp_hycom.append(timestamp)

    lonmin_raw = np.min(lon)
    lonmax_raw = np.max(lon)
    print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

    #================================================================
    # Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
    lon[lon>180] = lon[lon>180] - 360
    lon[lon<-180] = lon[lon<-180] + 360
    sort_lon = np.argsort(lon)
    lon = lon[sort_lon]

    # define grid boundaries
    lonmin_new = np.min(lon)
    lonmax_new = np.max(lon)
    latmin = np.min(lat)
    latmax = np.max(lat)
    print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

    # Sort OHC according with sorted lon 
    ohc_hafs = ohc_hafs[:,sort_lon]

    #%% Get storm track from trak atcf files
    file_track = folder + '/' + storm_id + '.' + cycle + '.' + hafs[i] + '.trak.atcfunix'

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

    #############  OHC large scale
    if np.min(lon) < 0:
        oklon = np.where(np.logical_and(lon>lon_lim[0],lon<lon_lim[1]))[0]
    else:
        oklon = np.where(np.logical_and(lon>lon_lim[0],lon<lon_lim[1]))[0]
    oklat = np.where(np.logical_and(lat>lat_lim[0],lat<lat_lim[1]))[0]

    #############  OHC along forecasted track
    '''
    oklon = np.round(np.interp(lon_forec_track[i,:],lon,np.arange(len(lon)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track[i,:],lat,np.arange(len(lat)))).astype(int)
    lon_forec_track_hafs[i,:] = lon[oklon]
    lat_forec_track_hafs[i,:] = lat[oklat]

    for x in np.arange(lon_forec_track.shape[1]):
        ohc_forec_track_hafs[i,x] = ohc_hafs[oklat,:][:,oklon][x,x]
    '''

    lon_forec_track_interp = np.interp(lat,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.interp(lon,lon_forec_track[i,:],lat_forec_track[i,:],left=np.nan,right=np.nan)
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    if len(lat_forec_track_int) > len(lon_forec_track_int):
        lon_forec_track_interp = np.copy(lon)
        lon_forec_track_interp[np.isnan(lat_forec_track_interp)] = np.nan
    else:
        lat_forec_track_interp = np.copy(lat)
        lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_hafs[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_hafs[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    oklon = np.round(np.interp(lon_forec_track_int,lon,np.arange(len(lon)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat,np.arange(len(lat)))).astype(int)

    for x in np.arange(len(lon_forec_track_int)):
        ohc_forec_track_hafs[i,x] = ohc_hafs[oklat[x],oklon[x]]
        
        lon[oklon[x]]
        lat[oklat[x]]
        dlon_lower = np.mean(np.diff(lon_nesdis))
        dlat_lower = np.mean(np.diff(lat_nesdis))
        dlon_higher = np.mean(np.diff(lon[oklon[x]-10:oklon[x]+10]))
        dlat_higher = np.mean(np.diff(lat[oklat[x]-10:oklat[x]+10]))
        resolx_diff = int((dlon_lower/dlon_higher-1)/2)
        resoly_diff = int((dlat_lower/dlat_higher-1)/2)

        ohc_forec_track_hafs_mean[i,x] = np.nanmean(ohc_hafs[oklat[x]-resoly_diff:oklat[x]+resoly_diff+1,:][:,oklon[x]-resolx_diff:oklon[x]+resolx_diff+1])
        ohc_forec_track_hafs_max[i,x] = np.nanmax(ohc_hafs[oklat[x]-resoly_diff:oklat[x]+resoly_diff+1,:][:,oklon[x]-resolx_diff:oklon[x]+resolx_diff+1])
        ohc_forec_track_hafs_min[i,x] = np.nanmin(ohc_hafs[oklat[x]-resoly_diff:oklat[x]+resoly_diff+1,:][:,oklon[x]-resolx_diff:oklon[x]+resolx_diff+1])
        
    #################################
    # Spatial interpolation of ohc_hafs to the NESDIS spatial resolution
    lono,lato = np.meshgrid(lon_nesdis,lat_nesdis)
    lonh,lath = np.meshgrid(lon,lat)
    
    interpolator = LinearNDInterpolator(list(zip(np.ravel(lonh),np.ravel(lath))),np.ravel(ohc_hafs))
    ohc_from_hafs_to_nesdis = interpolator((lono,lato))
    
    ##############################################
    #%% Figure OHC all domain HAFS
    lev = np.arange(0,161,20)
    kw = dict(levels=np.arange(0,161,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon,lat,ohc_hafs,lev,colors='grey',alpha=0.5)
    plt.contour(lon,lat,ohc_hafs,[100],colors='k')
    plt.contourf(lon,lat,ohc_hafs,cmap='Spectral_r',extend='max',**kw)
    cbar = plt.colorbar(extendrect=True)
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('OHC HAFS '+ cycle)
    plt.legend()
    #plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
    #plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    ##############################################
    #%% Figure OHC all domain HAFS on NESDIS resolution
    kw = dict(levels=np.arange(0,161,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon_nesdis,lat_nesdis,ohc_from_hafs_to_nesdis,lev,colors='grey',alpha=0.5)
    plt.contour(lon_nesdis,lat_nesdis,ohc_from_hafs_to_nesdis,[100],colors='k',alpha=0.5)
    plt.contourf(lon_nesdis,lat_nesdis,ohc_from_hafs_to_nesdis,cmap='Spectral_r',extend='max',**kw)
    cbar = plt.colorbar(extendrect=True)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('Ocean Heat Content HAFS to NESDIS resolution '+ str(time_nesdis)[2:15])
    plt.legend()

    ########################################################
    # difference between OHC NESDIS - OHC HAFS
    kw = dict(levels=np.arange(-50,51,5))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon_nesdis,lat_nesdis,ohc_nesdis-ohc_from_hafs_to_nesdis,cmap='Spectral_r',extend='both',**kw)
    cbar = plt.colorbar(extendrect=True)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('OHC NESDIS - OHC HAFS '+ str(time_nesdis)[2:15])
    plt.legend()

################################################################################
#%% NESDIS ohc along storm path
ohc_forec_track_nesdis = np.empty((len(exp_names),100))
ohc_forec_track_nesdis[:] = np.nan
lon_forec_track_nesdis = np.empty((len(exp_names),100))
lon_forec_track_nesdis[:] =np.nan
lat_forec_track_nesdis = np.empty((len(exp_names),100))
lat_forec_track_nesdis[:] =np.nan
for i in np.arange(len(exp_names)):
    '''
    oklon = np.round(np.interp(lon_forec_track[i,:],lon_nesdis,np.arange(len(lon_nesdis)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track[i,:],lat_nesdis,np.arange(len(lat_nesdis)))).astype(int)
    lon_forec_track_nesdis[i,:] = lon_nesdis[oklon]
    lat_forec_track_nesdis[i,:] = lat_nesdis[oklat]

    for x in np.arange(lon_forec_track.shape[1]):
        ohc_forec_track_nesdis[i,x] = ohc_nesdis[oklat[x],oklon[x]]
    '''

    ###########
    lon_forec_track_interp = np.interp(lat_nesdis,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.interp(lon_nesdis,lon_forec_track[i,:],lat_forec_track[i,:],left=np.nan,right=np.nan)
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    if len(lat_forec_track_int) > len(lon_forec_track_int):
        lon_forec_track_interp = np.copy(lon_nesdis)
        lon_forec_track_interp[np.isnan(lat_forec_track_interp)] = np.nan
    else:
        lat_forec_track_interp = np.copy(lat_nesdis)
        lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_nesdis[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_nesdis[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    oklon = np.round(np.interp(lon_forec_track_int,lon_nesdis,np.arange(len(lon_nesdis)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_nesdis,np.arange(len(lat_nesdis)))).astype(int)

    for x in np.arange(len(lon_forec_track_int)):
        ohc_forec_track_nesdis[i,x] = ohc_nesdis[oklat[x],oklon[x]]

#################################################################################
#%% Figure OHC all domain NESDIS

lev = np.arange(0,161,20)
kw = dict(levels=np.arange(0,161,20))
plt.figure(figsize=(8,5))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contour(lon_nesdis,lat_nesdis,ohc_nesdis,lev,colors='grey',alpha=0.5)
plt.contour(lon_nesdis,lat_nesdis,ohc_nesdis,[100],colors='k',alpha=0.5)
plt.contourf(lon_nesdis,lat_nesdis,ohc_nesdis,cmap='Spectral_r',extend='max',**kw)
cbar = plt.colorbar(extendrect=True)
cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('Ocean Heat Content NESDIS '+ str(time_nesdis)[2:15])
plt.legend()
#plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% time series sst along storm path
for i in np.arange(len(exp_names)):
    dist_forec_track_hafs = np.sqrt((lon_forec_track_hafs[i,:]-lon_forec_track_hafs[i,0])**2 + (lat_forec_track_hafs[i,:]-lat_forec_track_hafs[i,0])**2)*111
    dist_forec_track_nesdis = np.sqrt((lon_forec_track_nesdis[i,:]-lon_forec_track_nesdis[i,0])**2 + (lat_forec_track_nesdis[i,:]-lat_forec_track_nesdis[i,0])**2)*111

    fig,ax = plt.subplots(figsize=(10, 5))
    plt.plot(dist_forec_track_hafs,ohc_forec_track_hafs_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=5)
    ax.fill_between(dist_forec_track_hafs,ohc_forec_track_hafs_max[i,:],ohc_forec_track_hafs_min[i,:],color=exp_colors[i],alpha=0.3)
    plt.plot(dist_forec_track_nesdis,ohc_forec_track_nesdis[i,:],'o-',color='orange',label='NESDIS OHC',markeredgecolor='k',markersize=5)
    ax.fill_between(dist_forec_track_nesdis,ohc_forec_track_nesdis[i,:]+13.5,ohc_forec_track_nesdis[i,:]-13.5,color='orange',alpha=0.3)
    plt.legend()
    plt.xlabel('Along Track Distance ($km$)',fontsize=14)
    plt.ylabel('($kJ/cm^2$)',fontsize=14)
    plt.title('Ocean Heat Content Along Track',fontsize=14)

    fig,ax = plt.subplots(figsize=(10, 5))
    plt.plot(lon_forec_track_hafs[i,:],ohc_forec_track_hafs[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    plt.plot(lon_forec_track_nesdis[i,:],ohc_forec_track_nesdis[i,:],'o-',color='k',label='NESDIS OHC',markeredgecolor='r',markersize=7)
    plt.legend()
    plt.xlabel('Longitude',fontsize=14)
    plt.ylabel('($kJ/cm^2$)',fontsize=14)
    plt.title('Ocean Heat Content Along Track',fontsize=14)

    plt.figure()
    plt.plot(lon_forec_track_hafs[i,:],lat_forec_track_hafs[i,:],'*-',color=exp_colors[i],label=exp_labels[i])
    plt.plot(lon_forec_track_nesdis[i,:],lat_forec_track_nesdis[i,:],'.-',label='forec track NESDIS')
    plt.plot(lon_forec_track[i,:],lat_forec_track[i,:],'.-',color='k',label=exp_labels[i])
    plt.xlabel('Longitude',fontsize=14)
    plt.ylabel('Latitude',fontsize=14)
    plt.legend()
