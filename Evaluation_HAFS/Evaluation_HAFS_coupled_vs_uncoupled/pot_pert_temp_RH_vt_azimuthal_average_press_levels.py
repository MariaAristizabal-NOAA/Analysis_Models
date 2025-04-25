
#%% User input
#cycle = '2021063018'
cycle = '2021090206'

scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/AOML/aoml-hafs1/Lew.Gramer/ocean/HAFS_A/'

exp_names = ['HAFSv0p2a_coupled','HAFSv0p2a_uncoupled']
exp_labels = ['Coupled HAFS-A','Uncoupled HAFS-A']
exp_colors = ['blue','red']
#fhours = ['00','66']
fhours = ['75','93']

exp_folders = [scratch_folder1 + exp_names[0], scratch_folder1 + exp_names[1]]
name_trak_file = 'natl00l.' + cycle + '.trak.hafs.atcfunix.all'
files_adeck = [scratch_folder2 + 'coupled/adeck/' + name_trak_file, scratch_folder2 + '/uncoupled/adeck/' + name_trak_file] 

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

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

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_storm_track_and_int_high_resolution

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
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i],tan_theta[i],cos_theta[i] = Haversine(lat_track,lon_track,lat_lon_matrix[1,i],lat_lon_matrix[0,i])
   
    return oklat, oklon, R, tan_theta, cos_theta

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

#################################################################################
#%% Get storm track from trak atcf files
'''
lon_forec_track = np.empty((len(exp_names),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(exp_names),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(exp_names),43))
lead_time[:] = np.nan
int_track = np.empty((len(exp_names),43))
int_track[:] = np.nan
for f,ff in enumerate(files_adeck):
    lon_forec_track[f,:], lat_forec_track[f,:], lead_time[f,:], int_track[f,:] = get_storm_track_and_int(ff)
'''
#################################################################################
# Initiate figures profile to compare coupled vs uncoupled case

delta_r = 20 # Kilometers
ndr = int(560/delta_r) + 1

fig,apt1 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+ ' km) Perturbation Temperature ',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
apt1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

fig,apt2 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+ ' km) Perturbation Temperature ',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
apt2.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

fig,arh1 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+ ' km) Relative Humidity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Relative Humidity ($\%$)',fontsize=14)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
apt2.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

fig,arh2 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+ ' km) Relative Humidity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Relative Humidity ($\%$)',fontsize=14)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
apt2.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

fig,atv1 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+' km) Tangential Velocity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Tangential Velocity ($m/s$)',fontsize=14)
plt.legend()
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
atv1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

fig,atv2 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+' km) Tangential Velocity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Tangential Velocity ($m/s$)',fontsize=14)
plt.legend()
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
atv2.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

#########################################################################################
#%% Loop though hafs files of interest
for h,hh in enumerate(fhours):
    print('time #', h)

    for f,ff in enumerate(exp_folders):

        file_hafs_TMP = sorted(glob.glob(os.path.join(ff,cycle,'*0' + hh + '*TMP.nc')))[0]
        file_hafs_RH = sorted(glob.glob(os.path.join(ff,cycle,'*0' + hh + '*RH.nc')))[0]
        file_hafs_UGRD = sorted(glob.glob(os.path.join(ff,cycle,'*0' + hh + '*UGRD.nc')))[0]
        file_hafs_VGRD = sorted(glob.glob(os.path.join(ff,cycle,'*0' + hh + '*VGRD.nc')))[0]

        #################################################################################
        #%% Get storm track from trak atcf files
        file_adeck = files_adeck[f]
        lon_forec_track, lat_forec_track, lead_time, int_track = get_storm_track_and_int(file_adeck)
        
        #################################################################################
        #%% Get list name of variables to read
        fl = file_hafs_TMP
        hafs_tmp = xr.open_dataset(fl)
        sname = str(hafs_tmp.variables.items()).split('short_name')
        vars_list_TMP = []
        for n,sn in enumerate(sname[1:]): 
    	    if sn.split()[1].split('_')[0] == 'TMP':
                if sn.split()[1].split('_')[1][-2:] == 'mb':
            	    vars_list_TMP.append(sn.split()[1])

        fl = file_hafs_RH
        hafs_tmp = xr.open_dataset(fl)
        sname = str(hafs_tmp.variables.items()).split('short_name')
        vars_list_RH = []
        for n,sn in enumerate(sname[1:]): 
    	    if sn.split()[1].split('_')[0] == 'RH': 
                if sn.split()[1].split('_')[1][-2:] == 'mb':
            	    vars_list_RH.append(sn.split()[1])

        fl = file_hafs_UGRD
        hafs_tmp = xr.open_dataset(fl)
        sname = str(hafs_tmp.variables.items()).split('short_name')
        vars_list_UGRD = []
        for n,sn in enumerate(sname[1:]): 
    	    if sn.split()[1].split('_')[0] == 'UGRD': 
                if sn.split()[1].split('_')[1][-2:] == 'mb':
            	    vars_list_UGRD.append(sn.split()[1])

        fl = file_hafs_VGRD
        hafs_tmp = xr.open_dataset(fl)
        sname = str(hafs_tmp.variables.items()).split('short_name')
        vars_list_VGRD = []
        for n,sn in enumerate(sname[1:]): 
    	    if sn.split()[1].split('_')[0] == 'VGRD': 
                if sn.split()[1].split('_')[1][-2:] == 'mb':
            	    vars_list_VGRD.append(sn.split()[1])

        #################################################################################
        #%% Read profile of temperature and humidity from hafs files
        tmp_eye = np.empty((len(exp_names),len(vars_list_TMP)))
        tmp_eye[:] = np.nan
        tmp_mean550_650 = np.empty((len(exp_names),len(vars_list_TMP)))
        tmp_mean550_650[:] = np.nan
        rh_eye = np.empty((len(exp_names),len(vars_list_TMP)))
        rh_eye[:] = np.nan
        press = np.empty((len(exp_names),len(vars_list_TMP)))
        press[:] = np.nan
        time = np.empty((len(exp_names)))
        time[:] = np.nan

        tmp_mean_dr = np.empty((len(exp_names),len(vars_list_TMP),ndr))
        tmp_mean_dr[:] = np.nan
        pert_temp_dr = np.empty((len(exp_names),len(vars_list_TMP),ndr))
        pert_temp_dr[:] = np.nan
        rh_mean_dr = np.empty((len(exp_names),len(vars_list_RH),ndr))
        rh_mean_dr[:] = np.nan
        vt_mean_dr = np.empty((len(exp_names),len(vars_list_RH),ndr))
        vt_mean_dr[:] = np.nan

        hafs = xr.open_dataset(file_hafs_TMP)
        time = np.asarray(hafs.variables['time'][:])
        lat = np.asarray(hafs.variables['latitude'][:])
        lon = np.asarray(hafs.variables['longitude'][:])
        
        ind = np.where(lead_time == int(hh))[0][0]
        lon_track = lon_forec_track[ind]
        lat_track = lat_forec_track[ind]

        oklat, oklon, R, tan_theta, cos_theta = get_R_and_tan_theta_around_8degrees_from_eye(lat, lon, lon_track, lat_track)

    	# Obtaining points between 550 and 650 km of radius from eye
        okR550_650 = np.logical_and(R <= 650,R >= 550)

    	# Interpolating lon_track and lat_track into hafs grid
        oklon_track = np.int(np.round(np.interp(lon_track,lon,np.arange(len(lon)))))
        oklat_track = np.int(np.round(np.interp(lat_track,lat,np.arange(len(lat)))))

        #################################################################################
        #%% Perturbation temperature
        for z,var_name in enumerate(vars_list_TMP):
            hafs = xr.open_dataset(file_hafs_TMP)
            tmp_eye[f,z] = np.asarray(hafs.variables[var_name][0,oklat_track,oklon_track]) - 273.15
            tmp_mean550_650[f,z] = np.nanmean(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR550_650] - 273.15)
            press[f,z] = int(var_name.split('_')[1].split('mb')[0])

            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                tmp_mean_dr[f,z,dr] = np.nanmean(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR] - 273.15)  

        pert_temp_eye = tmp_eye - tmp_mean550_650

        pert_temp_dr = np.empty((len(exp_names),len(vars_list_TMP),ndr))
        pert_temp_dr[:] = np.nan
        for r in np.arange(ndr):
    	    pert_temp_dr[:,:,r] = tmp_mean_dr[:,:,r] - tmp_mean550_650

        #################################################################################
        #%% Humidity profiles
        for z,var_name in enumerate(vars_list_RH):
            hafs = xr.open_dataset(file_hafs_RH)
            rh_eye[f,z] = np.asarray(hafs.variables[var_name][0,oklat_track,oklon_track])

            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                rh_mean_dr[f,z,dr] = np.nanmean(np.ravel(np.asarray(hafs.variables[var_name][0,oklat,oklon]))[okR])  

	#################################################################################
	#%% Tangential velocity profiles
        for z in np.arange(len(vars_list_UGRD)):
            hafs_u = xr.open_dataset(file_hafs_UGRD)
            hafs_v = xr.open_dataset(file_hafs_VGRD)
            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                u = np.ravel(np.asarray(hafs_u.variables[vars_list_UGRD[z]][0,oklat,oklon]))[okR]     
                v = np.ravel(np.asarray(hafs_v.variables[vars_list_VGRD[z]][0,oklat,oklon]))[okR]
                vt = (v - tan_theta[okR] * u)/((1 + tan_theta[okR]**2) * cos_theta[okR])
                vt_mean_dr[f,z,dr] = np.nanmean(vt)  

    	#####################################################################################
        kw = dict(levels=np.arange(-14,14.1,1))
        fig,ax1 = plt.subplots(figsize=(6, 4))
        plt.contourf(np.arange(0,561,delta_r),press[f,:],pert_temp_dr[f,:,:],cmap='seismic',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$^oC$',fontsize = 14)
        plt.contour(np.arange(0,561,delta_r),press[f,:],pert_temp_dr[f,:,:],[0],colors='grey',alpha=0.2)
        plt.yscale('log')
        plt.ylim([50,1000])
        plt.yticks(np.arange(100,1001,100))
        ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
        plt.xticks(np.arange(0,551,50))
        plt.gca().invert_yaxis()
        #plt.colorbar()
        plt.title('Pertubation Temperature \n' + exp_names[f] + ' ' + str(time)[2:15])
        plt.ion()

    	#####################################################################################

        kw = dict(levels=np.arange(0,101,10))
        fig,ax1 = plt.subplots(figsize=(6, 4))
        plt.contourf(np.arange(0,561,delta_r),press[f,:],rh_mean_dr[f,:,:],cmap='Spectral_r',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$\%$',fontsize = 14)
        plt.contour(np.arange(0,561,delta_r),press[f,:],rh_mean_dr[f,:,:],[20],colors='grey')
        plt.yscale('log')
        plt.ylim([50,1000])
        plt.yticks(np.arange(100,1001,100))
        ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
        plt.xticks(np.arange(0,551,50))
        plt.gca().invert_yaxis()
        plt.title('Relative Humidity \n' + exp_names[f] + ' ' + str(time)[2:15])

        kw = dict(levels=np.arange(-10,71,5))
        fig,ax1 = plt.subplots(figsize=(6, 4))
        plt.contour(np.arange(0,561,delta_r),press[f,:],vt_mean_dr[f,:,:],colors='k',**kw)
        plt.contourf(np.arange(0,561,delta_r),press[f,:],vt_mean_dr[f,:,:],cmap='Spectral_r',**kw)
        plt.yscale('log')
        plt.ylim([50,1000])
        plt.yticks(np.arange(100,1001,100))
        ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
        plt.xticks(np.arange(0,551,50))
        plt.gca().invert_yaxis()
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$m/s$',fontsize = 14)
        plt.title('Tangential Velocity \n' +  exp_names[f] + ' ' + str(time)[2:15])

    	#####################################################################################
        if h == 0:
            apt1.plot(pert_temp_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            arh1.plot(rh_mean_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            atv1.plot(vt_mean_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
    
        if h == 1:
            apt2.plot(pert_temp_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            arh2.plot(rh_mean_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            atv2.plot(vt_mean_dr[f,:,0],press[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)

    apt1.legend()
    arh1.legend()
    atv1.legend()
    apt2.legend()
    arh2.legend()
    atv2.legend()

#####################################################################################
