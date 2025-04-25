# forecasting cycle to be used

# Danielle
cycle = '2022090212'
storm_num = '05'
basin = 'al'

#cycle = '2021082706'
#storm_num = '09'
#basin = 'al'

#exp_names = ['hafsv0p2a_2021rt_natl','HAFSv0p3_baseline','HAFSv0p3_HFAB_HF3A']
#exp_labels = ['HAFSv0.2a','H3BL','HF3A']
#exp_colors = ['limegreen','indianred','orange']

exp_names = ['hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_rt']
exp_colors = ['c']

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

#folder_exps = [scratch_folder + exp_names[0] + '/' + cycle,
#               scratch_folder + exp_names[1] + '/' + cycle,\
#               scratch_folder + exp_names[2] + '/' + cycle]

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
import cmocean
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

mlt_mean = np.empty((len(folder_exps),22))
mlt_mean[:] = np.nan
mlt_min = np.empty((len(folder_exps),22)) 
mlt_min[:] = np.nan
mlt_max = np.empty((len(folder_exps),22))
mlt_max[:] = np.nan
ohc_mean = np.empty((len(folder_exps),22))
ohc_mean[:] = np.nan
ohc_min = np.empty((len(folder_exps),22))
ohc_min[:] = np.nan
ohc_max = np.empty((len(folder_exps),22))
ohc_max[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    #if len(files_hafs_fv3) == 0:
    #    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs.grid01*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
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

    for n,file in enumerate(files_hafs_hycom):
       print(file)
       if np.logical_and(np.isfinite(lon_forec_track[i,:][::2][n]),np.isfinite(lat_forec_track[i,:][::2][n])):
           xlim = [lon_forec_track[i,:][::2][n]-4,lon_forec_track[i,:][::2][n]+4]
           ylim = [lat_forec_track[i,:][::2][n]-4,lat_forec_track[i,:][::2][n]+4]

           xlimh, ylimh  = geo_coord_to_HYCOM_coord(xlim,ylim)

           if np.min(lon_hafs_hycom) < 0:
               oklon = np.where(np.logical_and(lon_hafs_hycom>xlim[0],lon_hafs_hycom<xlim[1]))[0]
           else:
               oklon = np.where(np.logical_and(lon_hafs_hycom>xlimh[0],lon_hafs_hycom<xlimh[1]))[0]
           oklat = np.where(np.logical_and(lat_hafs_hycom>ylimh[0],lat_hafs_hycom<ylimh[1]))[0]

           hycom = xr.open_dataset(file)
           target_temp = np.asarray(hycom['temperature'][0,:,oklat,:][:,:,oklon])
           target_salt = np.asarray(hycom['salinity'][0,:,oklat,:][:,:,oklon])

           mlt_hycom = np.empty((len(oklat),len(oklon)))
           mlt_hycom[:] = np.nan
           for x in np.arange(len(oklon)):
               #print(x)
               for y in np.arange(len(oklat)):
                   _, mlt_hycom[y,x] = MLD_temp_crit(dtemp,ref_depth,depth_hafs_hycom,target_temp[:,y,x])

           mlt_mean[i,n] = np.nanmean(mlt_hycom)
           mlt_min[i,n] = np.nanmin(mlt_hycom)
           mlt_max[i,n] = np.nanmax(mlt_hycom)

           ohc_hycom = np.empty((len(oklat),len(oklon)))
           ohc_hycom[:] = np.nan
           for x in np.arange(len(oklon)):
               #print(x)
               for y in np.arange(len(oklat)):
                   dens_prof = sw.dens(target_salt[:,y,x],target_temp[:,y,x],depth_hafs_hycom)
                   ohc_hycom[y,x] = OHC_from_profile(depth_hafs_hycom,target_temp[:,y,x],dens_prof)

           ohc_mean[i,n] = np.nanmean(ohc_hycom)
           ohc_min[i,n] = np.nanmin(ohc_hycom)
           ohc_max[i,n] = np.nanmax(ohc_hycom)
       else:
           mlt_mean[i,n] = np.nan
           mlt_min[i,n] = np.nan
           mlt_mean[i,n] = np.nan
           ohc_min[i,n] = np.nan
           ohc_max[i,n] = np.nan
           ohc_max[i,n] = np.nan

       if n==8:
           lon_HAFS = lon_hafs_hycom[oklon]
           lat_HAFS = lat_hafs_hycom[oklat]

           okti = np.where(mdates.date2num(time_best_track) == mdates.date2num(time_hycom[0]))[0][0]
           oktf = np.where(mdates.date2num(time_best_track) == mdates.date2num(time_hycom[-1]))[0][0]
            
           kw = dict(levels=np.arange(22,32,0.5))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           plt.contour(lon_HAFS,lat_HAFS,mlt_hycom,[26],colors='k')
           plt.contourf(lon_HAFS,lat_HAFS,mlt_hycom,cmap='Spectral_r',**kw,extend='both')
           cbar = plt.colorbar(extendrect=True)
           cbar.ax.set_ylabel('$^oC$',fontsize = 14)
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[0,:][::2],'o-',color=exp_colors[0],markeredgecolor='k',label=exp_labels[0],markersize=7)
           plt.plot(lon_forec_track[i,:][::2][n], lat_forec_track[0,:][::2][n],'X',color='r',markeredgecolor='k',markersize=7)
           plt.plot(lon_best_track[okti:oktf], lat_best_track[okti:oktf],'o-',color='k',label='Best Track')
           plt.axis('scaled')
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('Mixed Layer Temp.' + exp_labels[i] + ' ' + str(time_fv3[n])[0:13])
           plt.legend(loc='lower left')

           kw = dict(levels=np.arange(0,170,20))
           plt.figure()
           plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
           plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
           plt.contour(lon_HAFS,lat_HAFS,ohc_hycom,[100],colors='k')
           plt.contour(lon_HAFS,lat_HAFS,ohc_hycom,[60],colors='grey')
           plt.contourf(lon_HAFS,lat_HAFS,ohc_hycom,cmap='Spectral_r',**kw,extend='max')
           cbar = plt.colorbar(extendrect=True)
           cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
           plt.plot(lon_forec_track[i,:][::2], lat_forec_track[0,:][::2],'o-',color=exp_colors[0],markeredgecolor='k',label=exp_labels[0],markersize=7)
           plt.plot(lon_forec_track[i,:][::2][n], lat_forec_track[0,:][::2][n],'X',color='red',markeredgecolor='k',markersize=7)
           plt.plot(lon_best_track[okti:oktf], lat_best_track[okti:oktf],'o-',color='k',label='Best Track')
           plt.axis('scaled')
           plt.ylim([np.min(lat_HAFS),np.max(lat_HAFS)])
           plt.xlim([np.min(lon_HAFS),np.max(lon_HAFS)])
           plt.title('OHC '+ exp_labels[i] + ' ' + str(time_fv3[n])[0:13])
           plt.legend(loc='lower left')

#################################################################################
#%% Figure mean MLT along storm track

fig,ax = plt.subplots(figsize = (8,5))

for i,folder in enumerate(folder_exps):
    plt.plot(lead_time[i,:][::2],mlt_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time[i,:][::2],mlt_min[i,:],mlt_max[i,:],color=exp_colors[i],alpha=0.1)

plt.legend(loc='lower left')
plt.title('Mixed Layer Temperature Cycle '+ cycle,fontsize=16)
plt.ylabel('$^oC$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

#################################################################################
#%% Figure mean OHC along storm track

fig,ax = plt.subplots(figsize = (8,5))

for i,folder in enumerate(folder_exps):
    plt.plot(lead_time[i,:][::2],ohc_mean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time[i,:][::2],ohc_min[i,:],ohc_max[i,:],color=exp_colors[i],alpha=0.1)

plt.legend()
plt.title('OHC Cycle '+ cycle,fontsize=16)
plt.ylabel('$kJ/cm^2$')
plt.xlabel('Forecast Lead Time (Hr)')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,127,12))

#################################################################################

