#%% User input
# forecasting cycle to be used

# Milton
cycle = '2024100606'
storm_num = '14'
storm_id = '14l'
basin = 'al'
#model = 'hfsa'
model = 'hfsb'
file_ohc = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/OSPO_OHC/2024/' + '/ohc_na14QG3_2024_280.nc'

#exp_names = ['HFSA_oper']
#exp_labels = ['HFSA']
#exp_colors = ['purple']
exp_names = ['HFSB_oper']
exp_labels = ['HFSB']
exp_colors = ['limegreen']

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_id]

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
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
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

sys.path.append(folder_uom)
from Upper_ocean_metrics import  OHC_from_profile

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
#%% Read OHC file
OHC = xr.open_dataset(file_ohc)

time_ospo = np.asarray(OHC['time'][:])
#okt = np.where(mdates.date2num(time_ospo) == mdates.date2num(datetime(tini.year,tini.month,tini.day)))[0][0]

lat_ospo =  np.asarray(OHC['latitude'][:])
lon_ospo =  np.asarray(OHC['longitude'][:])
#ohc_ospo =  np.asarray(OHC['ohc'][okt,:,:])
ohc_ospo =  np.asarray(OHC['ohc'][0,:,:])

#################################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

ohc_along_storm_track = np.empty((len(folder_exps),300))
ohc_along_storm_track[:] = np.nan
lon_forec_track_int_ocean = np.empty((len(folder_exps),300)) 
lon_forec_track_int_ocean[:] = np.nan
lat_forec_track_int_ocean = np.empty((len(folder_exps),300)) 
lat_forec_track_int_ocean[:] = np.nan
#%% Loop the experiments
for i,folder in enumerate(folder_exps):
    print(folder)

    oceanf = glob.glob(os.path.join(folder,'*f000.nc'))[0].split('/')[-1].split('.')
    ocean = [f for f in oceanf if f == 'hycom' or f == 'mom6'][0]

    #%% Get list files
    if ocean == 'mom6':
        file_hafs_ocean = folder + '/' + storm_id + '.' + cycle + '.' + model + '.mom6.f000.nc'
        nc = xr.open_dataset(file_hafs_ocean)
        temp = np.asarray(nc['temp'][0,:,:])
        salt = np.asarray(nc['so'][0,:,:])
        lon = np.asarray(nc['xh'])
        lat = np.asarray(nc['yh'])
        zl = np.asarray(nc['z_l'])
        time = np.asarray(nc['time'])
        ohc_hafs = ohc(temp,salt,zl)

    if ocean == 'hycom':
        file_hafs_ocean = folder + '/' + storm_id + '.' + cycle + '.' + model + '.hycom.3z.f000.nc'
        nc = xr.open_dataset(file_hafs_ocean)
        ohc_hafs = np.asarray(nc['ocean_heat_content'][0,:,:])
        ohc_hafs[ohc_hafs<-1000] = np.nan
        #var = np.asarray(nc['temperature'][0,0,:,:])
        lon = np.asarray(nc['Longitude'])
        lat = np.asarray(nc['Latitude'])
        depth = np.asarray(nc['Z'][:])
        time = np.asarray(nc['MT'][:])
        #timestamp = mdates.date2num(t)[0]
        #time_hycom.append(mdates.num2date(timestamp))
        #timestamp_hycom.append(timestamp)

    #%% Get storm track from trak atcf files
    file_track = folder + '/' + storm_id + '.' + cycle + '.' + model + '.trak.atcfunix'
    #14l.2024100606.hfsa.trak.atcfunix

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
 
    lon_limh, lat_limh  = geo_coord_to_HYCOM_coord(lon_lim,lat_lim)

    #############  OHC large scale
    if np.min(lon) < 0:
        oklon = np.where(np.logical_and(lon>lon_lim[0],lon<lon_lim[1]))[0]
    else:
        oklon = np.where(np.logical_and(lon>lon_limh[0],lon<lon_limh[1]))[0]
    oklat = np.where(np.logical_and(lat>lat_limh[0],lat<lat_limh[1]))[0]

    '''
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
    '''

    #############  OHC along forecasted track
    #%% HAFS-HYCOM OHC along storm path
    lon_forec_track_interp = np.interp(lat,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    lon_forec_track_int_ocean[i,0:len(lon_forec_track_int)] = lon_forec_track_int
    lat_forec_track_int_ocean[i,0:len(lat_forec_track_int)] = lat_forec_track_int

    #lon_forec_track_int_hycom_g, lat_forec_track_int_hycom_g  = geo_coord_to_HYCOM_coord(lon_forec_track_int_hycom,lat_forec_track_int_hycom)
    lon_forec_track_int_ocean_g = lon_forec_track_int_ocean[i,:]
    lat_forec_track_int_ocean_g = lat_forec_track_int_ocean[i,:]

    oklon = np.round(np.interp(lon_forec_track_int,lon,np.arange(len(lon)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat,np.arange(len(lat)))).astype(int)

    #target_temp = np.asarray(MODEL['temperature'][0,:,oklat,:][:,:,oklon])
    #target_salt = np.asarray(MODEL['salinity'][0,:,oklat,:][:,:,oklon])
    for x in np.arange(len(oklon)):
        ohc_along_storm_track = ohc_hafs[oklat,:][:,oklon][x,x]
        #dens_prof = sw.dens(target_salt[:,x,x],target_temp[:,x,x],depth_hafs_hycom)
        #ohc_along_storm_track_hycom[i,x]  = OHC_from_profile(depth_hafs_hycom,target_temp[:,x,x],dens_prof)
    ##############################################

    #%% Figure OHC all domain HAFS
    lev = np.arange(-9000,9100,100)
    kw = dict(levels=np.arange(0,161,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon,lat,ohc_hafs,lev,colors='grey',alpha=0.5)
    plt.contour(lon,lat,ohc_hafs,[100],colors='k')
    plt.contourf(lon,lat,ohc_hafs,cmap='Spectral_r',extend='max',**kw)
    cbar = plt.colorbar(extendrect=True)
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('OHC HAFS '+ cycle)
    #plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
    plt.legend()
    #plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

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
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('Ocean Heat Content OSPO '+ str(time_ospo)[2:15])
plt.legend()
#file_name = folder_fig + 'OHC_ospo' + file.split('/')[-1].split('.')[-2]
#plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

################################################################################
#%% time series sst along storm path

fig,ax = plt.subplots(figsize=(10, 5))
for i in np.arange(len(exp_names)):
    plt.plot(lat_forec_track_int_ocean[i,:],ohc_along_storm_track_ocean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    plt.plot(lat_forec_track_int_ospo[i,:],ohc_along_storm_track_ospo[i,:],'o-',color='k',label='OSPO OHC',markeredgecolor='k',markersize=7)
#'o-',color='c',label='hafsv0.3a',markeredgecolor='k',markersize=7)
plt.legend()
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.ylabel('($kJ/cm^2$)',fontsize=14)
plt.title('Ocean Heat Content Along Track',fontsize=14)


