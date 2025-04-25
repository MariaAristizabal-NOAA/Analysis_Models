
#%% User input

# Earl
#model = 'hfsb'
#storm_num = '06'
#storm_id = '06l'
#cycle = '2022090618'
#fhours = ['000','036','060']

# Fiona 
model = 'hfsb'
storm_num = '07'
storm_id = '07l'
cycle = '2022091900'
fhours = ['024','066']
#fhours = ['000','054','066']

# Larry 
#model = 'hfsb'
#storm_num = '12'
#storm_id = '12l'
#cycle = '2021090206'
#fhours = ['000','060','096']

# Ian
#model = 'hfsb'
#storm_num = '09'
#storm_id = '09l'
#cycle = '2022092706'
#fhours = ['000','030','042']

# Ida
#model = 'hfsb'
#storm_num = '09'
#storm_id = '09l'
#cycle = '2021082718'
#fhours = ['000','042','054']

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

#################################################################################
# Initiate figures profile to compare coupled vs uncoupled case

#delta_r = 20 # Kilometers
delta_r = 15 # Kilometers
ndr = int(560/delta_r) + 1

'''
fig,apt1 = plt.subplots(figsize=(6, 4))
plt.title('Mean (0-'+str(delta_r)+ ' km) Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.yticks(np.arange(0,20,2))
plt.ylim([0,20])

fig,apt2 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,20])
'''
'''
fig_arh1,arh1 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Relative Humidity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Relative Humidity ($\%$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])

fig_arh2,arh2 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Relative Humidity',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Relative Humidity ($\%$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
'''
'''
fig_atv1,atv1 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+' km) Tangential Velocity',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Tangential Velocity ($m/s$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])

fig_atv2,atv2 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+' km) Tangential Velocity',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Tangential Velocity ($m/s$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
'''

fig_appt1,appt1 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
plt.xlim([0,20])

fig_appt2,appt2 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
plt.xlim([0,20])

fig_appt3,appt3 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
plt.xlim([0,20])

'''
fig_apptt1,apptt1 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,7,1))
plt.ylim([0,6])

fig_apptt2,apptt2 = plt.subplots(figsize=(7,5))
plt.title('Mean (0-'+str(delta_r)+ ' km) Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,7,1))
plt.ylim([0,6])

fig_apptm1,apptm1 = plt.subplots(figsize=(7,5))
plt.title('Mean 200-300 km Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
plt.xlim([295,365])

fig_apptm2,apptm2 = plt.subplots(figsize=(7,5))
plt.title('Mean 200-300 km Potential Perturbation Temperature ',fontsize=16)
plt.ylabel('Height ($Km$)',fontsize = 14)
plt.xlabel('Temperature ($K$)',fontsize=14)
plt.yticks(np.arange(0,21,2))
plt.ylim([0,14])
plt.xlim([295,365])
'''

#########################################################################################
#%% Loop though hafs files of interest
for h,hh in enumerate(fhours):
    print('time #', h)

    #################################################################################
    #%% Read profile of temperature and humidity from hafs files
    nvert_lev = 45
    press = np.empty((len(exp_names),nvert_lev))
    press[:] = np.nan
    height = np.empty((len(exp_names),nvert_lev))
    height[:] = np.nan
    #time = np.empty((len(exp_names)))
    #time[:] = np.nan
    pot_pert_temp_dr = np.empty((len(exp_names),nvert_lev,ndr))
    pot_pert_temp_dr[:] = np.nan
    rh_mean_dr = np.empty((len(exp_names),nvert_lev,ndr))
    rh_mean_dr[:] = np.nan
    vt_mean_dr = np.empty((len(exp_names),nvert_lev,ndr))
    vt_mean_dr[:] = np.nan
    vr_mean_dr = np.empty((len(exp_names),nvert_lev,ndr))
    vr_mean_dr[:] = np.nan
    w_mean_dr = np.empty((len(exp_names),nvert_lev,ndr))
    w_mean_dr[:] = np.nan

    for f,ff in enumerate(exp_folders):
        file_hafs_fv3 = ff + 'grb2/' + cycle + '/' + storm_id + '.' + cycle + '.hfsb.parent.atm.f' + hh + '.grb2'

        #################################################################################
        #%% Get storm track from trak atcf files
        file_adeck = ff + 'adeck/' + track_file
        lon_forec_track, lat_forec_track, lead_time, int_track, rmw_track = get_storm_track_and_int(file_adeck,storm_num)
        
        #################################################################################
        #%% Read profile of temperature and humidity from hafs files
        hafs_fv3 = xr.open_dataset(file_hafs_fv3,engine="pynio")

        t0 = hafs_fv3['TMP_P0_L100_GLL0'].attrs['initial_time']
        dt = hafs_fv3['TMP_P0_L100_GLL0'].attrs['forecast_time'][0]
        t = datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt))
        timestamp = mdates.date2num(t)

        lon_hafs_fv3 = np.asarray(hafs_fv3.lon_0) - 360
        lat_hafs_fv3 = np.asarray(hafs_fv3.lat_0)
        press_hafs_fv3 = np.asarray(hafs_fv3.lv_ISBL0)/100 
        press[f,:] = press_hafs_fv3

        ind = np.where(lead_time == int(hh))[0][0]
        lon_track = lon_forec_track[ind]
        lat_track = lat_forec_track[ind]

        oklat, oklon, R, tan_theta, cos_theta, sin_theta = get_R_and_tan_theta_around_8degrees_from_eye(lat_hafs_fv3, lon_hafs_fv3, lon_track, lat_track)

    	# Obtaining points between 550 and 650 km of radius from eye
        okR550_650 = np.logical_and(R <= 650,R >= 550)
        okR200_300 = np.logical_and(R <= 300,R >= 200)

    	# Obtaining points at rmw
        okn = np.where(lead_time == int(hh))[0][0]
        ok_rmw = np.logical_and(R <= rmw_track[okn]+1,R >= rmw_track[okn]-1)

    	# Interpolating lon_track and lat_track into hafs grid
        oklon_track = int(np.round(np.interp(lon_track,lon_hafs_fv3,np.arange(len(lon_hafs_fv3)))))
        oklat_track = int(np.round(np.interp(lat_track,lat_hafs_fv3,np.arange(len(lat_hafs_fv3)))))

        #################################################################################
        #%% Perturbation temperature
        height[f,:] = metpy.calc.pressure_to_height_std(press_hafs_fv3*units.millibar).magnitude #in Km
        temp_eye = np.asarray(hafs_fv3['TMP_P0_L100_GLL0'][:,oklat_track,oklon_track])
        pot_temp_eye = metpy.calc.potential_temperature(press_hafs_fv3 * units.mbar, temp_eye * units.kelvin).magnitude 

        for z in np.arange(len(press_hafs_fv3)):
            temp_mean550_650 = np.nanmean(np.ravel(np.asarray(hafs_fv3['TMP_P0_L100_GLL0'][z,oklat,oklon]))[okR550_650])
            temp_mean200_300 = np.nanmean(np.ravel(np.asarray(hafs_fv3['TMP_P0_L100_GLL0'][z,oklat,oklon]))[okR200_300])
            temp_rmw = np.nanmean(np.ravel(np.asarray(hafs_fv3['TMP_P0_L100_GLL0'][z,oklat,oklon]))[ok_rmw])
            pot_temp_mean550_650 = metpy.calc.potential_temperature(press_hafs_fv3[z] * units.mbar, temp_mean550_650 * units.kelvin).magnitude  
            pot_temp_mean200_300 = metpy.calc.potential_temperature(press_hafs_fv3[z] * units.mbar, temp_mean200_300 * units.kelvin).magnitude 
            pot_temp_rmw = metpy.calc.potential_temperature(press_hafs_fv3[z] * units.mbar, temp_rmw * units.kelvin).magnitude  
            pert_temp_eye = temp_eye[z] - temp_mean200_300
            pot_pert_temp_eye = pot_temp_eye[z] - pot_temp_mean200_300

            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                temp_mean_dr = np.nanmean(np.ravel(np.asarray(hafs_fv3['TMP_P0_L100_GLL0'][z,oklat,oklon]))[okR])
                pot_temp_mean_dr = metpy.calc.potential_temperature(press_hafs_fv3[z] * units.mbar, temp_mean_dr * units.kelvin).magnitude
                pert_temp_dr = temp_mean_dr - temp_mean200_300
                pot_pert_temp_dr[f,z,dr] = pot_temp_mean_dr - pot_temp_mean200_300

        #################################################################################
        #%% Humidity profiles
        for z in np.arange(len(press_hafs_fv3)):
            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                rh_mean_dr[f,z,dr] = np.nanmean(np.ravel(hafs_fv3['RH_P0_L100_GLL0'][z,oklat,oklon])[okR])  

	#################################################################################
	#%% Tangential velocity profiles
        for z in np.arange(len(press_hafs_fv3)):
            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                u = np.ravel(np.asarray(hafs_fv3['UGRD_P0_L100_GLL0'][z,oklat,oklon]))[okR]     
                v = np.ravel(np.asarray(hafs_fv3['VGRD_P0_L100_GLL0'][z,oklat,oklon]))[okR]
                vt = (v - tan_theta[okR] * u)/((1 + tan_theta[okR]**2) * cos_theta[okR])
                vr = (u * cos_theta[okR] + v * sin_theta[okR]) 
                vt_mean_dr[f,z,dr] = np.nanmean(vt)  
                vr_mean_dr[f,z,dr] = np.nanmean(vr)  

	#################################################################################
	#%% Vertical velocity profiles
        for z in np.arange(len(press_hafs_fv3)):
            for dr in np.arange(ndr):
                okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
                w = np.ravel(np.asarray(hafs_fv3['DZDT_P0_L100_GLL0'][z,oklat,oklon]))[okR]
                w_mean_dr[f,z,dr] = np.nanmean(w)  

	#################################################################################
        xx = (np.arange(0,561,delta_r) + np.arange(20,561+20,delta_r))/2
        yy = height[f,:]
        X = np.tile(xx,(len(yy),1))
        Y = np.tile(yy,(len(xx),1)).T
        Z = np.tile(np.zeros(len(xx)),(len(yy),1))

        kw = dict(levels=np.arange(-20,20.1,1))
        fig,ax1 = plt.subplots(figsize=(7, 5))
        plt.contourf(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],cmap='seismic',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$K$',fontsize = 14)
        ct = plt.contour(np.arange(0,561,delta_r),height[f,:],vr_mean_dr[f,:,:],[0,2,4],colors='grey',alpha=0.5)
        ct = plt.contour(np.arange(0,561,delta_r),height[f,:],vr_mean_dr[f,:,:],[-8,-4,-2,-1],colors='k',alpha=0.5)
        plt.clabel(ct, ct.levels, inline=True, fontsize=12)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[0],colors='blue',alpha=0.1)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[15],colors='green',alpha=0.2)
        plt.yticks(np.arange(0,21,2))
        plt.ylim([0,14])
        plt.xticks(np.arange(0,551,50))
        plt.xlim([0,500])
        plt.ylabel('Height ($Km$)',fontsize = 14)
        plt.xlabel('Radius ($Km$)',fontsize = 14)
        plt.title('Potential Pertubation Temperature \n' + exp_names[f] + ' ' + cycle + ' ' + 'f' + str(hh))
        plt.savefig(dir_figures + 'PPT_'+cycle+'_'+ff.split('/')[-1]+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)

        '''
        kw = dict(levels=np.arange(-20,20.1,1))
        fig,ax1 = plt.subplots(figsize=(7,5))
        plt.contourf(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],cmap='seismic',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$K$',fontsize = 14)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[0],colors='blue',alpha=0.1)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[15],colors='green',alpha=0.2)
        Q = plt.quiver(X[::2,::2],Y[::2,::2],vr_mean_dr[f,::2,::2],w_mean_dr[f,::2,::2],scale=70) #,**kw)
        plt.quiverkey(Q,X=1.05,Y=1.05,U=10,label='10 m/s')
        plt.yticks(np.arange(0,21,2))
        #plt.ylim([0,20])
        plt.ylim([0,14])
        plt.xticks(np.arange(0,551,50))
        plt.xlim([0,500])
        plt.ylabel('Height ($Km$)',fontsize = 14)
        plt.xlabel('Radius ($Km$)',fontsize = 14)
        plt.title('Potential Pertubation Temperature \n' + exp_names[f] + ' ' + str(time)[2:15])
        plt.ion()
        plt.savefig(dir_figures + 'PPT2_'+cycle+'_'+ff.split('/')[-1]+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)

        kw = dict(levels=np.arange(-20,20.1,1))
        fig,ax1 = plt.subplots(figsize=(7,5))
        plt.contourf(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],cmap='seismic',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$K$',fontsize = 14)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[0],colors='blue',alpha=0.1)
        plt.contour(np.arange(0,561,delta_r),height[f,:],pot_pert_temp_dr[f,:,:],[15],colors='green',alpha=0.2)
        Q = plt.quiver(X[::2,::2],Y[::2,::2],vr_mean_dr[f,::2,::2],w_mean_dr[f,::2,::2],scale=70) #,**kw)
        plt.quiverkey(Q,X=1.05,Y=1.05,U=10,label='10 m/s')
        plt.yticks(np.arange(0,7,1))
        plt.ylim([0,6])
        plt.xticks(np.arange(0,551,50))
        plt.xlim([0,500])
        plt.ylabel('Height ($Km$)',fontsize = 14)
        plt.xlabel('Radius ($Km$)',fontsize = 14)
        plt.title('Potential Pertubation Temperature \n' + exp_names[f] + ' ' + str(time)[2:15])
        plt.savefig(dir_figures + 'PPT3_'+cycle+'_'+ff.split('/')[-1]+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        '''
    	#####################################################################################
        '''
        kw = dict(levels=np.arange(0,101,10))
        fig,ax1 = plt.subplots(figsize=(7,5))
        plt.contourf(np.arange(0,561,delta_r),height[f,:],rh_mean_dr[f,:,:],cmap='Spectral_r',**kw)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$\%$',fontsize = 14)
        plt.contour(np.arange(0,561,delta_r),height[f,:],rh_mean_dr[f,:,:],[20],colors='grey')
        plt.yticks(np.arange(0,21,2))
        #plt.ylim([0,20])
        plt.ylim([0,14])
        plt.xticks(np.arange(0,551,50))
        plt.xlim([0,500])
        plt.ylabel('Height ($Km$)',fontsize = 14)
        plt.xlabel('Radius ($Km$)',fontsize = 14)
        plt.title('Relative Humidity \n' + exp_names[f] + ' ' + cycle + ' ' + 'f' + str(hh))
        #plt.savefig(dir_figures + 'RH_'+cycle+'_'+ff.split('/')[-1]+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        '''

    	#####################################################################################

        '''
        kw = dict(levels=np.arange(-10,71,5))
        fig,ax1 = plt.subplots(figsize=(7,5))
        plt.contour(np.arange(0,561,delta_r),height[f,:],vt_mean_dr[f,:,:],colors='k',**kw)
        plt.contourf(np.arange(0,561,delta_r),height[f,:],vt_mean_dr[f,:,:],cmap='Spectral_r',**kw)
        plt.yticks(np.arange(0,21,2))
        #plt.ylim([0,20])
        plt.ylim([0,14])
        plt.xticks(np.arange(0,551,50))
        plt.xlim([0,500])
        plt.ylabel('Height ($Km$)',fontsize = 14)
        plt.xlabel('Radius ($Km$)',fontsize = 14)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('$m/s$',fontsize = 14)
        plt.title('Tangential Velocity \n' +  exp_names[f] + ' ' + str(time)[2:15])
        plt.savefig(dir_figures + 'VT_'+cycle+'_'+ff.split('/')[-1]+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        '''

    	#####################################################################################
        if h == 0:
            #apt1.plot(pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            appt1.plot(pot_pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + ' ' + cycle + ' ' + 'f' + str(hh),markeredgecolor='k',markersize=7)
            #apptt1.plot(pot_pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #apptm1.plot(pot_temp_mean200_300[f,:],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #arh1.plot(rh_mean_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #atv1.plot(vt_mean_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
    
        if h == 1:
            #apt2.plot(pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + str(time)[2:15],markeredgecolor='k',markersize=7)
            appt2.plot(pot_pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + ' ' + cycle + ' ' + 'f' + str(hh),markeredgecolor='k',markersize=7)
            #apptt2.plot(pot_pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #apptm2.plot(pot_temp_mean200_300[f,:],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #arh2.plot(rh_mean_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)
            #atv2.plot(vt_mean_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + '_' + str(time)[2:15],markeredgecolor='k',markersize=7)

        if h == 2:
            appt3.plot(pot_pert_temp_dr[f,:,0],height[f,:],'X-',color=exp_colors[f],label=exp_names[f] + ' ' + cycle + ' ' + 'f' + str(hh),markeredgecolor='k',markersize=7)

    if h == 0:
        appt1.legend()
        fig_appt1.savefig(dir_figures + 'PPT_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #apptt1.legend()
        #fig_apptt1.savefig(dir_figures + 'PPT2_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #apptm1.legend()
        #fig_apptm1.savefig(dir_figures + 'PPTM_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #arh1.legend()
        #fig_arh1.savefig(dir_figures + 'RH_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #atv1.legend()
        #fig_atv1.savefig(dir_figures + 'VT_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
    if h == 1:
        appt2.legend()
        fig_appt2.savefig(dir_figures + 'PPT_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #apptt2.legend()
        #fig_apptt2.savefig(dir_figures + 'PPT2_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #apptm2.legend()
        #fig_apptm2.savefig(dir_figures + 'PPTM_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #arh2.legend()
        #fig_arh2.savefig(dir_figures + 'RH_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
        #atv2.legend()
        #fig_atv2.savefig(dir_figures + 'VT_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)

    if h == 2:
        appt3.legend()
        fig_appt3.savefig(dir_figures + 'PPT_'+cycle+'_f'+hh,bbox_inches = 'tight',pad_inches = 0.1)
#####################################################################################
