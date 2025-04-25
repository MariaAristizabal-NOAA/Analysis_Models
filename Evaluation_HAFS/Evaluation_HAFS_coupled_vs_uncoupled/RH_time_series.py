
#%% User input
storm_num = '06'
cycle = '2022090618'

#storm_num = '07'
#cycle = '2022091900'

#storm_num = '09'
#cycle = '2022092618'

#storm_num = '09'
#cycle = '2021082800'

scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/AOML/aoml-hafs1/Lew.Gramer/ocean/hafsv0p3s_2022rt_flux_fix/'

exp_names = ['hafsv0p3s_2022rt_flux_fix/coupled','hafsv0p3s_2022rt_flux_fix/uncoupled']
exp_labels = ['Coupled HAFS-S','Uncoupled HAFS-S']
exp_colors = ['blue','red']

exp_folders = [scratch_folder1 + exp_names[0], scratch_folder1 + exp_names[1]]
name_trak_file = storm_num + 'l.' + cycle + '.hafs.grid02.trak.atcfunix.all'
files_adeck = [scratch_folder2 + 'coupled/adeck/' + name_trak_file, scratch_folder2 + '/uncoupled/adeck/' + name_trak_file]

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

dir_figures = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Figures_coupled_uncoupled/'

#%%
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import netCDF4
import glob
import os
from datetime import datetime,timedelta
import matplotlib.dates as mdates
import sys
import metpy.calc
from metpy.units import units
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine

#################################################################################

def get_R_and_tan_theta_around_8degrees_from_eye(lat, lon, lon_track, lat_track):

    xlim = [lon_track-8,lon_track+8]
    ylim = [lat_track-8,lat_track+8]

    oklon = np.where(np.logical_and(lon>xlim[0],lon<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat>ylim[0],lat<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon[oklon],lat[oklat])

    eye_lon = np.tile(lon_track,meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_track,meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    delta_x = np.diff(lon[oklon])[0]
    delta_y = np.diff(lat[oklat])[0]
    area = delta_x * delta_y * np.cos(lat[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area)

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    tan_theta = np.empty(lat_lon_matrix.shape[1])
    tan_theta[:] = np.nan
    cos_theta = np.empty(lat_lon_matrix.shape[1])
    cos_theta[:] = np.nan
    sin_theta = np.empty(lat_lon_matrix.shape[1])
    sin_theta[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i],tan_theta[i],cos_theta[i], sin_theta[i] = Haversine(lat_track,lon_track,lat_lon_matrix[1,i],lat_lon_matrix[0,i])
   
    return oklat, oklon, R, tan_theta, cos_theta, sin_theta

#################################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= 16,bath_lat <= 20)
oklonbath = np.logical_and(bath_lon >= -67.5,bath_lon <= -65)

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#########################################################################################
max_rh_0_15km= np.empty((len(exp_names),43))
max_rh_0_15km[:] = np.nan
min_rh_0_15km= np.empty((len(exp_names),43))
min_rh_0_15km[:] = np.nan
mean_rh_0_15km= np.empty((len(exp_names),43))
mean_rh_0_15km[:] = np.nan

delta_r = 15 # Kilometers
ndr = int(560/delta_r) + 1
nvert_lev = 45

for f,ff in enumerate(exp_folders):
    print(ff)

    #%% Get storm track from trak atcf files
    file_adeck = files_adeck[f]
    lon_forec_track, lat_forec_track, lead_time, int_track, rmw_track = get_storm_track_and_int(file_adeck,storm_num)
        
    #%% Get file for all forecast times
    files_hafs_RH = sorted(glob.glob(os.path.join(ff,cycle,'*RH.nc')))

    #%% Loop throug forecast times
    for h,fl in enumerate(files_hafs_RH):
        print(h)
        #%% Get list name of variables to read
        hafs_rh = xr.open_dataset(fl)
        sname = str(hafs_rh.variables.items()).split('short_name')
        vars_list_RH = []
        for n,sn in enumerate(sname[1:]): 
    	    if sn.split()[1].split('_')[0] == 'RH':
                if sn.split()[1].split('_')[1][-2:] == 'mb':
            	    vars_list_RH.append(sn.split()[1])

        #################################################################################
        #%% Read profile of humidity from hafs files
        hafs = xr.open_dataset(fl)
        time = np.asarray(hafs.variables['time'][:])
        lat = np.asarray(hafs.variables['latitude'][:])
        lon = np.asarray(hafs.variables['longitude'][:])
        
        #ind = np.where(lead_time == int(hh))[0][0]
        lon_track = lon_forec_track[h]
        lat_track = lat_forec_track[h]

        oklat, oklon, R, tan_theta, cos_theta, sin_theta = get_R_and_tan_theta_around_8degrees_from_eye(lat, lon, lon_track, lat_track)

    	# Obtaining points between 0 and 15 km of radius from eye
        okR0_15 = np.logical_and(R <= 15,R >= 0)

    	# Interpolating lon_track and lat_track into hafs grid
        oklon_track = int(np.round(np.interp(lon_track,lon,np.arange(len(lon)))))
        oklat_track = int(np.round(np.interp(lat_track,lat,np.arange(len(lat)))))

        #################################################################################
        #%% 
        meanrh_0_15km = np.empty((nvert_lev))
        meanrh_0_15km[:] = np.nan
        maxrh_0_15km = np.empty((nvert_lev))
        maxrh_0_15km[:] = np.nan
        minrh_0_15km = np.empty((nvert_lev))
        minrh_0_15km[:] = np.nan
        height = np.empty((nvert_lev))
        height[:] = np.nan
        for z,var_name in enumerate(vars_list_RH):
            hafs = xr.open_dataset(fl)
            pres = int(var_name.split('_')[1].split('mb')[0])
            height[z] = metpy.calc.pressure_to_height_std(pres*units.millibar).magnitude #in Km
            meanrh_0_15km[z] = np.nanmean(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR0_15])
            maxrh_0_15km[z] = np.nanmax(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR0_15])
            minrh_0_15km[z] = np.nanmin(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR0_15])
            
        okh = np.where(height <= 0.12)[0]
        if okh.size > 0:
            mean_rh_0_15km[f,h] = np.mean(meanrh_0_15km[okh])
            max_rh_0_15km[f,h] = np.max(maxrh_0_15km[okh])
            min_rh_0_15km[f,h] = np.min(minrh_0_15km[okh])
        else:
            mean_rh_0_15km[f,h] = np.nan
            max_rh_0_15km[f,h] = np.nan
            min_rh_0_15km[f,h] = np.nan

#####################################################################################

fig,ax = plt.subplots(figsize=(10, 5))
for f,ff in enumerate(exp_folders):
    plt.plot(lead_time,mean_rh_0_15km[f,:],'X-',color=exp_colors[f],label=exp_names[f] ,markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time,max_rh_0_15km[f,:],min_rh_0_15km[f,:],color=exp_colors[f],alpha=0.1)

plt.legend()
plt.title('Near Surface Relative Humidity '+storm_num +'L '+cycle,fontsize=16)
plt.ylabel('(%)',fontsize = 14)
plt.xlabel('Forecast Lead Time',fontsize=14)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.savefig(dir_figures + 'RH_near_surface_'+cycle,bbox_inches = 'tight',pad_inches = 0.1)

