
#%% User input

#model = 'hfsb'
#storm_num = '07'
#storm_id = '07l'
#cycle = '2022091900'

#model = 'hfsb'
#storm_num = '12'
#storm_id = '12l'
#cycle = '2021090206'

model = 'hfsb'
storm_num = '09'
storm_id = '09l'
cycle = '2022092706'

scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/AOML/aoml-hafs1/Lew.Gramer/ocean/hafsv1_fnl_hfsb/'

exp_names = ['coupled','uncoupled']
exp_labels = ['Coupled HFSB','Uncoupled HFSB']
exp_colors = ['blue','red']

exp_folders = [scratch_folder2 + exp_names[0] + '/', scratch_folder2 + exp_names[1] + '/']

track_file = storm_id + '.' + cycle + '.' + model + '.storm.trak.atcfunix'

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

###############################################################################
Int_track = np.empty((len(exp_names),43))
Int_track[:] = np.nan
max_warm_core_anomaly = np.empty((len(exp_names),43))
max_warm_core_anomaly[:] = np.nan
min_warm_core_anomaly = np.empty((len(exp_names),43))
min_warm_core_anomaly[:] = np.nan
mean_warm_core_anomaly = np.empty((len(exp_names),43))
mean_warm_core_anomaly[:] = np.nan
std_warm_core_anomaly = np.empty((len(exp_names),43))
std_warm_core_anomaly[:] = np.nan
height_max_warm_core_anomaly = np.empty((len(exp_names),43))
height_max_warm_core_anomaly[:] = np.nan

delta_r = 15 # Kilometers
ndr = int(560/delta_r) + 1
nvert_lev = 45

for f,ff in enumerate(exp_folders):
    print(ff)

    #%% Get storm track from trak atcf files
    file_adeck = ff + 'adeck/' + track_file
    lon_forec_track, lat_forec_track, lead_time, int_track, rmw_track = get_storm_track_and_int(file_adeck,storm_num)
        
    #%% Get file for all forecast times
    files_hafs = sorted(glob.glob(os.path.join(ff,'grb2',cycle,'*.hfsb.parent.atm.f*.grb2')))

    #%% Loop throug forecast times
    for h,fl in enumerate(files_hafs):
        print(h)

        #%% Read profile of temperature from hafs files
        hafs = xr.open_dataset(fl,engine="pynio")
 
        t0 = hafs['TMP_P0_L100_GLL0'].attrs['initial_time']
        dt = hafs['TMP_P0_L100_GLL0'].attrs['forecast_time'][0]
        time = datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt))
        timestamp = mdates.date2num(time)

        lon = np.asarray(hafs.lon_0) - 360
        lat = np.asarray(hafs.lat_0)
        press = np.asarray(hafs.lv_ISBL0)/100
        height = metpy.calc.pressure_to_height_std(press*units.millibar).magnitude #in Km
        
        lon_track = lon_forec_track[h]
        lat_track = lat_forec_track[h]

        oklat, oklon, R, tan_theta, cos_theta, sin_theta = get_R_and_tan_theta_around_8degrees_from_eye(lat, lon, lon_track, lat_track)

    	# Obtaining points between 200 and 300 km of radius from eye
        okR200_300 = np.logical_and(R <= 300,R >= 200)
        okR0_15 = np.logical_and(R <= 15,R >= 0)

    	# Interpolating lon_track and lat_track into hafs grid
        oklon_track = int(np.round(np.interp(lon_track,lon,np.arange(len(lon)))))
        oklat_track = int(np.round(np.interp(lat_track,lat,np.arange(len(lat)))))

        #################################################################################
        #%% Perturbation temperature
        warm_core_anomaly = np.empty((nvert_lev))
        warm_core_anomaly[:] = np.nan
        for z in np.arange(len(press)):
            pres = press[z]
            temp_mean200_300 = np.nanmean(np.ravel(np.asarray(hafs['TMP_P0_L100_GLL0'][z,oklat,oklon]))[okR200_300])
            pot_temp_mean200_300 = metpy.calc.potential_temperature(pres * units.mbar, temp_mean200_300 * units.kelvin).magnitude - 273.15 
            temp_mean0_15 = np.nanmean(np.ravel(np.asarray(hafs['TMP_P0_L100_GLL0'][z,oklat,oklon]))[okR0_15])
            pot_temp_mean0_15 = metpy.calc.potential_temperature(pres * units.mbar, temp_mean0_15 * units.kelvin).magnitude - 273.15 
            warm_core_anomaly[z] = pot_temp_mean0_15 - pot_temp_mean200_300
        
        mca = np.max(warm_core_anomaly) 
        max_warm_core_anomaly[f,h] = mca
        ok_mca = np.where(warm_core_anomaly == mca)[0][0]
        height_max_warm_core_anomaly[f,h] = height[ok_mca]  
        min_warm_core_anomaly[f,h] = np.min(warm_core_anomaly)
        mean_warm_core_anomaly[f,h] = np.mean(warm_core_anomaly)
        std_warm_core_anomaly[f,h] = np.std(warm_core_anomaly)

    Int_track[f,:] = int_track
        
#####################################################################################

fig,ax = plt.subplots(figsize=(10, 5))
for f,ff in enumerate(exp_folders):
    plt.plot(lead_time,max_warm_core_anomaly[f,:],'X-',color=exp_colors[f],label=exp_names[f] ,markeredgecolor='k',markersize=7)
plt.legend()
plt.title('Max. Warm Core Anomaly '+storm_num +'L '+cycle,fontsize=16)
plt.ylabel('($^oC$)',fontsize = 14)
plt.xlabel('Forecast Lead Time',fontsize=14)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.savefig(dir_figures + 'max_warm_core_anomaly_'+cycle,bbox_inches = 'tight',pad_inches = 0.1)

fig,ax = plt.subplots(figsize=(10, 5))
for f,ff in enumerate(exp_folders):
    plt.plot(lead_time,mean_warm_core_anomaly[f,:],'X-',color=exp_colors[f],label=exp_names[f] ,markeredgecolor='k',markersize=7)
    ax.fill_between(lead_time,mean_warm_core_anomaly[f,:]-std_warm_core_anomaly[f,:],mean_warm_core_anomaly[f,:]+std_warm_core_anomaly[f,:],color=exp_colors[f],alpha=0.1)
plt.legend()
plt.title('Warm Core Anomaly '+storm_num +'L '+cycle,fontsize=16)
plt.ylabel('($^oC$)',fontsize = 14)
plt.xlabel('Forecast Lead Time',fontsize=14)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.savefig(dir_figures + 'std_warm_core_anomaly_'+cycle,bbox_inches = 'tight',pad_inches = 0.1)

###################################################################
# Correlation

fig,ax = plt.subplots(figsize=(7, 6))
for f,ff in enumerate(exp_folders):
    r = np.corrcoef(Int_track[f,:],mean_warm_core_anomaly[f,:])
    plt.plot(Int_track[f,:],mean_warm_core_anomaly[f,:],'X',color=exp_colors[f],label=exp_names[f]+' r='+str(np.round(r[0,1],2)) ,markeredgecolor='k',markersize=7)
plt.legend()
plt.title(storm_num +'L '+cycle,fontsize=16)
plt.ylabel('Mean Warm Core Anomaly ($^oC$)',fontsize = 14)
plt.xlabel('Intensity (kt)',fontsize=14)
plt.savefig(dir_figures + 'int_vs_warm_core_anomaly_'+cycle,bbox_inches = 'tight',pad_inches = 0.1)


'''
fig,ax = plt.subplots(figsize=(10, 5))
for f,ff in enumerate(exp_folders):
    plt.plot(lead_time,height_max_warm_core_anomaly[f,:],'X-',color=exp_colors[f],label=exp_names[f] ,markeredgecolor='k',markersize=7)
plt.legend()
plt.title('Height Max. Warm Core Anomaly '+storm_num +'L '+cycle,fontsize=16)
plt.ylabel('($Km$)',fontsize = 14)
plt.xlabel('Forecast Lead Time',fontsize=14)
ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,126,12))
plt.savefig(dir_figures + 'height_max_warm_core_anomaly_'+cycle,bbox_inches = 'tight',pad_inches = 0.1)
'''
