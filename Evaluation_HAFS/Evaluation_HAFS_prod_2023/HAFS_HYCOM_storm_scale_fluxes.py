#%% User input
# forecasting cycle to be used

# Danielle
cycle = '2022090212'
storm_num = '05'
basin = 'al'

exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start','hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_fixed_wind_stresses','HAFSv0.3a_rt']
exp_colors = ['limegreen','c']

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',
               scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/'] 

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = scratch_folder + 'adeck/a' + basin + storm_num + cycle[0:4] + '.dat'

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

okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

#################################################################################
#%% Read GFS track
#lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)

#################################################################################
#%% Loop the experiments

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan
shtfl_mean = np.empty((len(folder_exps),43))
shtfl_mean[:] = np.nan
shtfl_min = np.empty((len(folder_exps),43)) 
shtfl_min[:] = np.nan
shtfl_max = np.empty((len(folder_exps),43))
shtfl_max[:] = np.nan
lhtfl_mean = np.empty((len(folder_exps),43))
lhtfl_mean[:] = np.nan
lhtfl_min = np.empty((len(folder_exps),43))
lhtfl_min[:] = np.nan
lhtfl_max = np.empty((len(folder_exps),43))
lhtfl_max[:] = np.nan
uflx_mean = np.empty((len(folder_exps),43))
uflx_mean[:] = np.nan
uflx_min = np.empty((len(folder_exps),43))
uflx_min[:] = np.nan
uflx_max = np.empty((len(folder_exps),43))
uflx_max[:] = np.nan
vflx_mean = np.empty((len(folder_exps),43))
vflx_mean[:] = np.nan
vflx_min = np.empty((len(folder_exps),43))
vflx_min[:] = np.nan
vflx_max = np.empty((len(folder_exps),43))
vflx_max[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    #files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    #if files:
    #    file_track = files[0]
    #else:
    #    file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    file_track = folder + storm_num + basin[-1] + '.' + cycle + '.hafs.trak.atcfunix'
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    #if len(files_hafs_fv3) == 0:
    #    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs.grid01*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Reading HAFS/FV3 grid
    #hafs_fv3_grid = xr.open_dataset(files_hafs_fv3[0],decode_times=False)
    #lon_hafs_fv3 = np.asarray(hafs_fv3_grid['longitude'][:])
    #lat_hafs_fv3 = np.asarray(hafs_fv3_grid['latitude'][:])
     
    fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
    lon_hafs_fv3 = np.asarray(fv3.lon_0) - 360
    lat_hafs_fv3 = np.asarray(fv3.lat_0)
    #t0 = FV3.variables['TMP_P0_L1_GLL0'].attrs['initial_time']
    #dt = FV3.variables['TMP_P0_L1_GLL0'].attrs['forecast_time'][0]
    #time_fv3.append(datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt)))

    #%% Read HAFS/Fv3 time
    '''
    time_fv3 = []
    for n,file in enumerate(files_hafs_fv3):
        print(file)
        FV3 = xr.open_dataset(file)
        t = FV3.variables['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_fv3.append(mdates.num2date(timestamp))

    time_fv3 = np.asarray(time_fv3)
    '''

    for n,file in enumerate(files_hafs_fv3):
       print(file)

       if np.logical_and(np.isfinite(lon_forec_track[i,:][n]),np.isfinite(lat_forec_track[i,:][n])):
           xlim = [lon_forec_track[i,:][n]-2,lon_forec_track[i,:][n]+2]
           ylim = [lat_forec_track[i,:][n]-2,lat_forec_track[i,:][n]+2]

           oklon = np.where(np.logical_and(lon_hafs_fv3>xlim[0],lon_hafs_fv3<xlim[1]))[0]
           oklat = np.where(np.logical_and(lat_hafs_fv3>ylim[0],lat_hafs_fv3<ylim[1]))[0]

           #fv3 = xr.open_dataset(file)
           #target_shtfl = np.asarray(fv3['SHTFL_surface'][0,oklat,:][:,oklon])
           #target_lhtfl = np.asarray(fv3['LHTFL_surface'][0,oklat,:][:,oklon])
           #target_uflx = np.asarray(fv3['UFLX_surface'][0,oklat,:][:,oklon])
           #target_vflx = np.asarray(fv3['VFLX_surface'][0,oklat,:][:,oklon])
           #target_ugrd = np.asarray(fv3['UGRD_10maboveground'][0,oklat,:][:,oklon])
           #target_vgrd = np.asarray(fv3['VGRD_10maboveground'][0,oklat,:][:,oklon])

           fv3 = xr.open_dataset(file,engine="pynio")
           #target_shtfl = np.asarray(fv3['SHTFL_P8_L1_GLL0_avg'])[oklat,:][:,oklon]
           #target_lhtfl = np.asarray(fv3['LHTFL_P8_L1_GLL0_avg'])[oklat,:][:,oklon]
           target_lhtfl = np.asarray(fv3['LHTFL_P0_L1_GLL0'])[oklat,:][:,oklon]
           target_shtfl = np.asarray(fv3['SHTFL_P0_L1_GLL0'])[oklat,:][:,oklon]
           target_uflx = np.asarray(fv3['UFLX_P0_L1_GLL0'])[oklat,:][:,oklon]
           target_vflx = np.asarray(fv3['VFLX_P0_L1_GLL0'])[oklat,:][:,oklon] 
           target_ugrd = np.asarray(fv3['UGRD_P0_L103_GLL0'])[oklat,:][:,oklon]
           target_vgrd = np.asarray(fv3['VGRD_P0_L103_GLL0'])[oklat,:][:,oklon]

           shtfl_mean[i,n] = np.nanmean(target_shtfl)
           shtfl_min[i,n] = np.nanmin(target_shtfl)
           shtfl_max[i,n] = np.nanmax(target_shtfl)

           lhtfl_mean[i,n] = np.nanmean(target_lhtfl)
           lhtfl_min[i,n] = np.nanmin(target_lhtfl)
           lhtfl_max[i,n] = np.nanmax(target_lhtfl)

           uflx_mean[i,n] = np.nanmean(target_uflx)
           uflx_min[i,n] = np.nanmin(target_uflx)
           uflx_max[i,n] = np.nanmax(target_uflx)

           vflx_mean[i,n] = np.nanmean(target_vflx)
           vflx_min[i,n] = np.nanmin(target_vflx)
           vflx_max[i,n] = np.nanmax(target_vflx)        
       else:
           shtfl_mean[i,n] = np.nan
           shtfl_min[i,n] = np.nan
           shtfl_max[i,n] = np.nan

           lhtfl_mean[i,n] = np.nan
           lhtfl_min[i,n] = np.nan
           lhtfl_max[i,n] = np.nan

           uflx_mean[i,n] = np.nan
           uflx_min[i,n] = np.nan
           uflx_max[i,n] = np.nan

           vflx_mean[i,n] = np.nan
           vflx_min[i,n] = np.nan
           vflx_max[i,n] = np.nan

       
       if n==20:
           ff = file.split('/')[-1].split('.')[-2]
           #############################################
           lon_HAFS = lon_hafs_fv3[oklon]
           lat_HAFS = lat_hafs_fv3[oklat]

           kw = dict(levels=np.arange(-200,310,50))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           ctr = plt.contourf(lon_HAFS,lat_HAFS,target_shtfl,cmap='Spectral_r',**kw,extend='both')
           cbar = plt.colorbar(ctr,extendrect=True)
           cbar.ax.set_ylabel('$W/m^2$',fontsize = 14)
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
           plt.plot(lon_best_track[okt],lat_best_track[okt],'o-',color='k',label='Best Track')
           plt.plot(lon_forec_track[i,:][n], lat_forec_track[i,:][n],'o-',color='red',markeredgecolor='k',markersize=7)
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('Sensible Heat Flux '+ exp_labels[i] + ' ' + ff)
           plt.legend(loc='lower left')
           file_name = storm_num + '.' + cycle + '.storm.sshfl.' + ff + '.' + exp_names[i]
           plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

           #############################################
           kw = dict(levels=np.arange(0,1010,100))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           ctr = plt.contourf(lon_HAFS,lat_HAFS,target_lhtfl,cmap='Spectral_r',**kw,extend='both')
           cbar = plt.colorbar(ctr,extendrect=True)
           cbar.ax.set_ylabel('$W/m^2$',fontsize = 14)
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
           plt.plot(lon_best_track[okt],lat_best_track[okt],'o-',color='k',label='Best Track')
           plt.plot(lon_forec_track[i,:][n], lat_forec_track[i,:][n],'o-',color='red',markeredgecolor='k',markersize=7)
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('Latent Heat Flux '+ exp_labels[i] + ' ' + ff)
           plt.legend(loc='lower left')
           file_name = storm_num + '.' + cycle + '.storm.lthfl.' + ff + '.' + exp_names[i]
           plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
       
           #########################################
           lon_HAFS = lon_hafs_fv3[oklon]
           lat_HAFS = lat_hafs_fv3[oklat]

           kw = dict(levels=np.arange(-300,310,50))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           ctr = plt.contourf(lon_HAFS,lat_HAFS,target_uflx,cmap='seismic',**kw,extend='both') 
           cbar = plt.colorbar(ctr,extendrect=True)
           cbar.ax.set_ylabel('$N/m^2$',fontsize = 14)
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
           plt.plot(lon_best_track[okt],lat_best_track[okt],'o-',color='k',label='Best Track')
           plt.plot(lon_forec_track[i,:][n], lat_forec_track[i,:][n],'o-',color='red',markeredgecolor='k',markersize=7)
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('U Momentum Flux '+ exp_labels[i] + ' ' + ff)
           plt.legend(loc='lower left')
           file_name = storm_num + '.' + cycle + '.storm.umfl.' + ff + '.' + exp_names[i]
           plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
           
           #########################################
           lon_HAFS = lon_hafs_fv3[oklon]
           lat_HAFS = lat_hafs_fv3[oklat]
           Vgrd = np.sqrt(target_ugrd**2 + target_vgrd**2)

           kw = dict(levels=np.arange(0,51,5))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           ctr = plt.contourf(lon_HAFS,lat_HAFS,Vgrd,cmap='Spectral_r',**kw,extend='max')
           cbar = plt.colorbar(ctr,extendrect=True)
           cbar.ax.set_ylabel('m/s',fontsize = 14)
           Q = plt.quiver(lon_HAFS[::10],lat_HAFS[::10],target_ugrd[::10,::10],target_vgrd[::10,::10]) #,**kw)
           plt.quiverkey(Q,X=1,Y=-0.1, U= 30.0,label='30 m/s')
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[i,:][::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
           plt.plot(lon_best_track[okt],lat_best_track[okt],'o-',color='k',label='Best Track')
           plt.plot(lon_forec_track[i,:][n], lat_forec_track[i,:][n],'o-',color='red',markeredgecolor='k',markersize=7)
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('Speed at 10 m Above Ground '+ exp_labels[i] + ' ' + ff)
           plt.legend(loc='lower left',bbox_to_anchor=[-0.1,-0.1])
           file_name = storm_num + '.' + cycle + '.storm.sp10m.' + ff + '.' + exp_names[i]
           plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% Figure mean SHTFL along storm track

fig,ax = plt.subplots(figsize = (8,5))

for i,folder in enumerate(folder_exps):
    plt.plot(lead_time[i,:],shtfl_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time[i,:],shtfl_min[i,:],shtfl_max[i,:],color=exp_colors[i],alpha=0.1)

plt.legend(loc='lower left')
plt.title('Sensible Heat Flux Cycle '+ cycle,fontsize=16)
plt.ylabel('$W/m^2$')
plt.xlabel('Forecast Lead Time (Hr)')
plt.ylim([-500,500])

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

file_name = storm_num + '.' + cycle + '.storm.snhfl.time_series.' + exp_names[i]
plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% Figure mean SLHTFL along storm track

fig,ax = plt.subplots(figsize = (8,5))

for i,folder in enumerate(folder_exps):
    plt.plot(lead_time[i,:],lhtfl_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time[i,:],lhtfl_min[i,:],lhtfl_max[i,:],color=exp_colors[i],alpha=0.1)

plt.legend(loc='lower left')
plt.title('Latent Heat Flux Cycle '+ cycle,fontsize=16)
plt.ylabel('$W/m^2$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

file_name = storm_num + '.' + cycle + '.storm.lthfl.time_series.' + exp_names[i]
plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% Figure mean UFLx along storm track

fig,ax = plt.subplots(figsize = (8,5))

for i,folder in enumerate(folder_exps):
    plt.plot(lead_time[i,:],uflx_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time[i,:],uflx_min[i,:],uflx_max[i,:],color=exp_colors[i],alpha=0.1)

plt.legend(loc='lower left')
plt.title('U Momentum Flux Cycle '+ cycle,fontsize=16)
plt.ylabel('$W/m^2$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

file_name = storm_num + '.' + cycle + '.storm.umfl.time_series.' + exp_names[i]
plt.savefig(file_name+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


