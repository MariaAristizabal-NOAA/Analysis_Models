
#%% User input
#cycle = '2019083112'
#cycle = '2019083100'
cycle = '2019082800'


Dir_HWRF_POM_oper = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.' + cycle + '_grb2_to_nc_oper/' 
Dir_HWRF_POM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.' + cycle +'_grb2_to_nc_exp/'
Dir_HWRF_HYCOM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.' + cycle +'_grb2_to_nc_exp/'

file_atcf_HWRF_POM_oper = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.2019082800_oper/dorian05l.2019082800.trak.hwrf.atcfunix'
#dorian05l.2019082800.track_d03.patcf
file_atcf_HWRF_POM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.2019082800_exp/dorian05l.2019082800.trak.hwrf.atcfunix'
#dorian05l.2019082800.track_d03.patcf
file_atcf_HWRF_HYCOM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.2019082800_exp/dorian05l.2019082800.trak.hwrf.atcfunix'
#dorian05l.2019082800.track_d03.patcf

file_track_HWRF_POM_oper = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.2019082800_oper/dorian05l.2019082800.track_d03.patcf'
file_track_HWRF_POM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.2019082800_exp/dorian05l.2019082800.track_d03.patcf'
file_track_HWRF_HYCOM_exp = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.2019082800_exp/dorian05l.2019082800.track_d03.patcf'

scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

dir_figs = '/home/Maria.Aristizabal/Dorian_2019/Figures/'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

#%%
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import netCDF4
import cmocean
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
#file_track_HWRF_POM_oper = sorted(glob.glob(os.path.join(Dir_HWRF_POM_oper,'*trak*')))[0]
#file_track_HWRF_POM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_POM_exp,'*trak*')))[0]
#file_track_HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_HYCOM_exp,'*trak*')))[0]

lon_forec_track_HWRF_POM_oper, lat_forec_track_HWRF_POM_oper, lead_time_HWRF_POM_oper, int_track_HWRF_POM_oper = get_storm_track_and_int(file_atcf_HWRF_POM_oper)
lon_forec_track_HWRF_POM_exp, lat_forec_track_HWRF_POM_exp, lead_time_HWRF_POM_exp, int_track_HWRF_POM_exp = get_storm_track_and_int(file_atcf_HWRF_POM_exp)
lon_forec_track_HWRF_HYCOM_exp, lat_forec_track_HWRF_HYCOM_exp, lead_time_HWRF_HYCOM_exp, int_track_HWRF_HYCOM_exp = get_storm_track_and_int(file_atcf_HWRF_HYCOM_exp)

lon_forec_track_HWRF_POM_oper_hr, lat_forec_track_HWRF_POM_oper_hr, lead_time_HWRF_POM_oper_hr, int_track_HWRF_POM_oper_hr,_,_ = get_storm_track_and_int_high_resolution(file_track_HWRF_POM_oper)
lon_forec_track_HWRF_POM_exp_hr, lat_forec_track_HWRF_POM_exp_hr, lead_time_HWRF_POM_exp_hr, int_track_HWRF_POM_exp_hr,_,_ = get_storm_track_and_int_high_resolution(file_track_HWRF_POM_exp)
lon_forec_track_HWRF_HYCOM_exp_hr, lat_forec_track_HWRF_HYCOM_exp_hr, lead_time_HWRF_HYCOM_exp_hr, int_track_HWRF_HYCOM_exp_hr,_,_ = get_storm_track_and_int_high_resolution(file_track_HWRF_HYCOM_exp)

#################################################################################

#%% Get list HWRF files
HWRF_POM_oper = sorted(glob.glob(os.path.join(Dir_HWRF_POM_oper,'*.nc')))
HWRF_POM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_POM_exp,'*.nc')))
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_HYCOM_exp,'*.nc')))

#################################################################################
#%% Get list name of variables to read
fl = HWRF_HYCOM_exp[0]
HWRF = xr.open_dataset(fl)
sname = str(HWRF.variables.items()).split('short_name')
vars_list_TMP = []
vars_list_RH = []
vars_list_UGRD = []
vars_list_VGRD = []
for n,sn in enumerate(sname[1:]): 
    if sn.split()[1].split('_')[0] == 'TMP': 
        vars_list_TMP.append(sn.split()[1])
    if sn.split()[1].split('_')[0] == 'RH': 
        vars_list_RH.append(sn.split()[1])
    if sn.split()[1].split('_')[0] == 'UGRD': 
        vars_list_UGRD.append(sn.split()[1])
    if sn.split()[1].split('_')[0] == 'VGRD': 
        vars_list_VGRD.append(sn.split()[1])

#################################################################################
#%% Read profile of temperature and humidity from HWRF files
tmp_HWRF_POM_oper_eye = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_POM_oper_eye[:] = np.nan
tmp_HWRF_POM_oper_mean550_650 = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_POM_oper_mean550_650[:] = np.nan
rh_HWRF_POM_oper_eye = np.empty((3,len(vars_list_TMP[0:-3])))
rh_HWRF_POM_oper_eye[:] = np.nan
press_HWRF_POM_oper = np.empty((3,len(vars_list_TMP[0:-3])))
press_HWRF_POM_oper[:] = np.nan
time_HWRF_POM_oper = np.empty((3))
time_HWRF_POM_oper[:] = np.nan

delta_r = 20 # Kilometers
ndr = int(560/delta_r) + 1
tmp_HWRF_POM_oper_mean_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
tmp_HWRF_POM_oper_mean_dr[:] = np.nan
rh_HWRF_POM_oper_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
rh_HWRF_POM_oper_mean_dr[:] = np.nan
vt_HWRF_POM_oper_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
vt_HWRF_POM_oper_mean_dr[:] = np.nan
#for t,indx in enumerate(np.asarray([2,4])):
for t,indx in enumerate(np.asarray([0])):
    print(HWRF_POM_oper[indx])
    HWRF = xr.open_dataset(HWRF_POM_oper[indx])
    time_HWRF_POM_oper[t] = np.asarray(HWRF.variables['time'][:])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])

    lon_track = lon_forec_track_HWRF_POM_oper[indx]
    lat_track = lat_forec_track_HWRF_POM_oper[indx]

    xlim = [lon_track-8,lon_track+8]
    ylim = [lat_track-8,lat_track+8]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_track,meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_track,meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    delta_x = np.diff(lon_hwrf[oklon])[0]
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    tan_theta = np.empty(lat_lon_matrix.shape[1])
    tan_theta[:] = np.nan
    cos_theta = np.empty(lat_lon_matrix.shape[1])
    cos_theta[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i],tan_theta[i],cos_theta[i] = Haversine(lat_track,lon_track,lat_lon_matrix[1,i],lat_lon_matrix[0,i])

    okR550_650 = np.logical_and(R <= 650,R >= 550)

    # Interpolating lon_track and lat_track into HYCOM grid
    oklon_track = np.int(np.round(np.interp(lon_track,lon_hwrf,np.arange(len(lon_hwrf)))))
    oklat_track = np.int(np.round(np.interp(lat_track,lat_hwrf,np.arange(len(lat_hwrf)))))

    for z,var_name in enumerate(vars_list_TMP[0:-3]):
        tmp_HWRF_POM_oper_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track]) - 275.15
        tmp_HWRF_POM_oper_mean550_650[t,z] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR550_650] - 275.15)
        press_HWRF_POM_oper[t,z] = int(var_name.split('_')[1].split('mb')[0])

        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            tmp_HWRF_POM_oper_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR] - 275.15)  

    for z,var_name in enumerate(vars_list_RH[0:-1]):
        rh_HWRF_POM_oper_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track])

        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            rh_HWRF_POM_oper_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR])  

    for z in np.arange(len(vars_list_UGRD[0:-2])):
        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            u = np.ravel(np.asarray(HWRF.variables[vars_list_UGRD[0:-2][z]][0,oklat,oklon]))[okR]     
            v = np.ravel(np.asarray(HWRF.variables[vars_list_VGRD[0:-2][z]][0,oklat,oklon]))[okR]
            vt = (v - tan_theta[okR] * u)/((1 + tan_theta[okR]**2) * cos_theta[okR])
            vt_HWRF_POM_oper_mean_dr[t,z,dr] = np.nanmean(vt)  

#################################################################################
#%% Read profile of temperature and humidity from HWRF files
tmp_HWRF_POM_exp_eye = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_POM_exp_eye[:] = np.nan
tmp_HWRF_POM_exp_mean550_650 = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_POM_exp_mean550_650[:] = np.nan
rh_HWRF_POM_exp_eye = np.empty((3,len(vars_list_TMP[0:-3])))
rh_HWRF_POM_exp_eye[:] = np.nan
press_HWRF_POM_exp = np.empty((3,len(vars_list_TMP[0:-3])))
press_HWRF_POM_exp [:] = np.nan
time_HWRF_POM_exp = np.empty((3))
time_HWRF_POM_exp[:] = np.nan
tmp_HWRF_POM_exp_mean_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
tmp_HWRF_POM_exp_mean_dr[:] = np.nan
rh_HWRF_POM_exp_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
rh_HWRF_POM_exp_mean_dr[:] = np.nan
vt_HWRF_POM_exp_mean_dr = np.empty((3,len(vars_list_UGRD[0:-2]),ndr))
vt_HWRF_POM_exp_mean_dr[:] = np.nan
vt_HWRF_POM_exp_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
vt_HWRF_POM_exp_mean_dr[:] = np.nan
#for t,indx in enumerate(np.asarray([2,4])):
for t,indx in enumerate(np.asarray([0])):
    print(HWRF_POM_exp[indx])
    HWRF = xr.open_dataset(HWRF_POM_exp[indx])
    time_HWRF_POM_exp[t] = np.asarray(HWRF.variables['time'][:])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])

    lon_track = lon_forec_track_HWRF_POM_exp[indx]
    lat_track = lat_forec_track_HWRF_POM_exp[indx]

    xlim = [lon_track-8,lon_track+8]
    ylim = [lat_track-8,lat_track+8]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_track,meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_track,meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    delta_x = np.diff(lon_hwrf[oklon])[0]
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    tan_theta = np.empty(lat_lon_matrix.shape[1])
    tan_theta[:] = np.nan
    cos_theta = np.empty(lat_lon_matrix.shape[1])
    cos_theta[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i],tan_theta[i],cos_theta[i] = Haversine(lat_track,lon_track,lat_lon_matrix[1,i],lat_lon_matrix[0,i])

    okR550_650 = np.logical_and(R <= 650,R >= 550)

    # Interpolating lon_track and lat_track into HYCOM grid
    oklon_track = np.int(np.round(np.interp(lon_track,lon_hwrf,np.arange(len(lon_hwrf)))))
    oklat_track = np.int(np.round(np.interp(lat_track,lat_hwrf,np.arange(len(lat_hwrf)))))

    for z,var_name in enumerate(vars_list_TMP[0:-3]):
        tmp_HWRF_POM_exp_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track]) - 275.15
        tmp_HWRF_POM_exp_mean550_650[t,z] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR550_650] - 275.15)
        press_HWRF_POM_exp[t,z] = int(var_name.split('_')[1].split('mb')[0])

        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            tmp_HWRF_POM_exp_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR] - 275.15)

    for z,var_name in enumerate(vars_list_RH[0:-1]):
        rh_HWRF_POM_exp_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track])

        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            rh_HWRF_POM_exp_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR]) 

    for z in np.arange(len(vars_list_UGRD[0:-2])):
        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            u = np.ravel(np.asarray(HWRF.variables[vars_list_UGRD[0:-2][z]][0,oklat,oklon]))[okR]
            v = np.ravel(np.asarray(HWRF.variables[vars_list_VGRD[0:-2][z]][0,oklat,oklon]))[okR]
            vt = (v - tan_theta[okR] * u)/((1 + tan_theta[okR]**2) * cos_theta[okR])
            vt_HWRF_POM_exp_mean_dr[t,z,dr] = np.nanmean(vt)

#################################################################################
#%% Read profile of temperature and humidity from HWRF file
tmp_HWRF_HYCOM_exp_eye = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_HYCOM_exp_eye[:] = np.nan
tmp_HWRF_HYCOM_exp_mean550_650 = np.empty((3,len(vars_list_TMP[0:-3])))
tmp_HWRF_HYCOM_exp_mean550_650[:] = np.nan
rh_HWRF_HYCOM_exp_eye = np.empty((3,len(vars_list_TMP[0:-3])))
rh_HWRF_HYCOM_exp_eye[:] = np.nan
press_HWRF_HYCOM_exp = np.empty((3,len(vars_list_TMP[0:-3])))
press_HWRF_HYCOM_exp[:] = np.nan
time_HWRF_HYCOM_exp = np.empty((3))
time_HWRF_HYCOM_exp[:] = np.nan
tmp_HWRF_HYCOM_exp_mean_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
tmp_HWRF_HYCOM_exp_mean_dr[:] = np.nan
rh_HWRF_HYCOM_exp_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
rh_HWRF_HYCOM_exp_mean_dr[:] = np.nan
vt_HWRF_HYCOM_exp_mean_dr = np.empty((3,len(vars_list_RH[0:-1]),ndr))
vt_HWRF_HYCOM_exp_mean_dr[:] = np.nan
#for t,indx in enumerate(np.asarray([2,4])):
for t,indx in enumerate(np.asarray([0])):
    print(HWRF_HYCOM_exp[indx])
    HWRF = xr.open_dataset(HWRF_HYCOM_exp[indx])
    time_HWRF_HYCOM_exp[t] = np.asarray(HWRF.variables['time'][:])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])

    lon_track = lon_forec_track_HWRF_HYCOM_exp[indx]
    lat_track = lat_forec_track_HWRF_HYCOM_exp[indx]

    xlim = [lon_track-8,lon_track+8]
    ylim = [lat_track-8,lat_track+8]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_track,meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_track,meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    delta_x = np.diff(lon_hwrf[oklon])[0]
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    tan_theta = np.empty(lat_lon_matrix.shape[1])
    tan_theta[:] = np.nan
    cos_theta = np.empty(lat_lon_matrix.shape[1])
    cos_theta[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i],tan_theta[i],cos_theta[i] = Haversine(lat_track,lon_track,lat_lon_matrix[1,i],lat_lon_matrix[0,i])

    okR550_650 = np.logical_and(R <= 650,R >= 550)

    # Interpolating lon_track and lat_track into HYCOM grid
    oklon_track = np.int(np.round(np.interp(lon_track,lon_hwrf,np.arange(len(lon_hwrf)))))
    oklat_track = np.int(np.round(np.interp(lat_track,lat_hwrf,np.arange(len(lat_hwrf)))))
        
    for z,var_name in enumerate(vars_list_TMP[0:-3]):
        tmp_HWRF_HYCOM_exp_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track]) - 275.15
        tmp_HWRF_HYCOM_exp_mean550_650[t,z] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR550_650] - 275.15)
        press_HWRF_HYCOM_exp[t,z] = int(var_name.split('_')[1].split('mb')[0])
        
        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            tmp_HWRF_HYCOM_exp_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR] - 275.15)

    for z,var_name in enumerate(vars_list_RH[0:-1]):
        rh_HWRF_HYCOM_exp_eye[t,z] = np.asarray(HWRF.variables[var_name][0,oklat_track,oklon_track]) 

        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            rh_HWRF_HYCOM_exp_mean_dr[t,z,dr] = np.nanmean(np.ravel(np.asarray(HWRF.variables[var_name][0,oklat,oklon]))[okR])

    for z in np.arange(len(vars_list_UGRD[0:-2])):
        for dr in np.arange(ndr):
            okR = np.logical_and(R >= delta_r*dr,R <= delta_r*(dr+1))
            u = np.ravel(np.asarray(HWRF.variables[vars_list_UGRD[0:-2][z]][0,oklat,oklon]))[okR]
            v = np.ravel(np.asarray(HWRF.variables[vars_list_VGRD[0:-2][z]][0,oklat,oklon]))[okR]
            vt = (v - tan_theta[okR] * u)/((1 + tan_theta[okR]**2) * cos_theta[okR])
            vt_HWRF_HYCOM_exp_mean_dr[t,z,dr] = np.nanmean(vt)

#################################################################################
#%% Perturbation temperature
pert_temp_HWRF_POM_oper_eye = tmp_HWRF_POM_oper_eye - tmp_HWRF_POM_oper_mean550_650
pert_temp_HWRF_POM_exp_eye = tmp_HWRF_POM_exp_eye - tmp_HWRF_POM_exp_mean550_650
pert_temp_HWRF_HYCOM_exp_eye = tmp_HWRF_HYCOM_exp_eye - tmp_HWRF_HYCOM_exp_mean550_650

pert_temp_HWRF_POM_oper_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
pert_temp_HWRF_POM_oper_dr[:] = np.nan
pert_temp_HWRF_POM_exp_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
pert_temp_HWRF_POM_exp_dr[:] = np.nan
pert_temp_HWRF_HYCOM_exp_dr = np.empty((3,len(vars_list_TMP[0:-3]),ndr))
pert_temp_HWRF_HYCOM_exp_dr[:] = np.nan
for r in np.arange(ndr):
    pert_temp_HWRF_POM_oper_dr[:,:,r] = tmp_HWRF_POM_oper_mean_dr[:,:,r] - tmp_HWRF_POM_oper_mean550_650
    pert_temp_HWRF_POM_exp_dr[:,:,r] = tmp_HWRF_POM_exp_mean_dr[:,:,r] - tmp_HWRF_POM_exp_mean550_650
    pert_temp_HWRF_HYCOM_exp_dr[:,:,r] = tmp_HWRF_HYCOM_exp_mean_dr[:,:,r] - tmp_HWRF_HYCOM_exp_mean550_650

#################################################################################
#%% Figures
'''
t = 0
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(tmp_HWRF_POM_oper_eye[t,:],-press_HWRF_POM_oper[t,:],\
'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(tmp_HWRF_POM_exp_eye[t,:],-press_HWRF_POM_exp[t,:],\
'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(tmp_HWRF_HYCOM_exp_eye[t,:],-press_HWRF_HYCOM_exp[t,:],\
'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Atmosphere Temperature profile at Storm Eye ',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.legend()

####################################################################################
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(pert_temp_HWRF_POM_oper_eye[t,:],-press_HWRF_POM_oper[t,:],\
'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(pert_temp_HWRF_POM_exp_eye[t,:],-press_HWRF_POM_exp[t,:],\
'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(pert_temp_HWRF_HYCOM_exp_eye[t,:],-press_HWRF_HYCOM_exp[t,:],\
'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Perturbation Temperature at Storm Eye ',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.legend()
plt.xlim([-4,16])

####################################################################################
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(tmp_HWRF_POM_oper_mean550_650[t,:],-press_HWRF_POM_oper[t,:],\
'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(tmp_HWRF_POM_exp_mean550_650[t,:],-press_HWRF_POM_exp[t,:],\
'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(tmp_HWRF_HYCOM_exp_mean550_650[t,:],-press_HWRF_HYCOM_exp[t,:],\
'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Mean Temperature within 550 km < R < 650 km',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.legend()
'''
#####################################################################################
t=0
kw = dict(levels=np.arange(-5,6,0.5))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_oper[t,:],pert_temp_HWRF_POM_oper_dr[t,:,:],vmin=-5,vmax=5,cmap='seismic',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Pertubation Temperature HWRF POM oper')

#####################################################################################
kw = dict(levels=np.arange(-5,6,0.5))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],pert_temp_HWRF_POM_exp_dr[t,:,:],vmin=-5,vmax=5,cmap='seismic',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Pertubation Temperature HWRF POM exp')

#####################################################################################
kw = dict(levels=np.arange(-5,6,0.5))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_HYCOM_exp[t,:],pert_temp_HWRF_HYCOM_exp_dr[t,:,:],vmin=-5,vmax=6,cmap='seismic',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Pertubation Temperature HWRF HYCOM exp')

####################################################################################
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(pert_temp_HWRF_POM_oper_dr[t,:,0],press_HWRF_POM_oper[t,:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(pert_temp_HWRF_POM_exp_dr[t,:,0],press_HWRF_POM_exp[t,:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(pert_temp_HWRF_HYCOM_exp_dr[t,:,0],press_HWRF_HYCOM_exp[t,:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Mean (0-'+str(delta_r)+ ' km) Perturbation Temperature',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.legend()
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

#####################################################################################
kw = dict(levels=np.arange(0,101,10))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_oper[t,:],rh_HWRF_POM_oper_mean_dr[t,:,:],vmin=10,vmax=90,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Relative Humidity HWRF POM oper')

#####################################################################################
kw = dict(levels=np.arange(0,101,10))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],rh_HWRF_POM_exp_mean_dr[t,:,:],vmin=10,vmax=90,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Relative Humidity HWRF POM exp')

#####################################################################################
kw = dict(levels=np.arange(0,101,10))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contourf(np.arange(0,561,delta_r),press_HWRF_HYCOM_exp[t,:],rh_HWRF_HYCOM_exp_mean_dr[t,:,:],vmin=10,vmax=90,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Relative Humidity HWRF HYCOM exp')

######################################################################################
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(rh_HWRF_POM_oper_mean_dr[t,:,0],press_HWRF_POM_oper[t,:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(rh_HWRF_POM_exp_mean_dr[t,:,0],press_HWRF_POM_exp[t,:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(rh_HWRF_HYCOM_exp_mean_dr[t,:,0],press_HWRF_HYCOM_exp[t,:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Mean (0-'+str(delta_r)+ ' km) Relative Humidity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Relative Humidity ($\%$)',fontsize=14)
plt.legend()
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

######################################################################################
kw = dict(levels=np.arange(-24,25,2))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contour(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_POM_oper_mean_dr[t,:,:],colors='k',**kw)
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_POM_oper_mean_dr[t,:,:],vmin=-24,vmax=24,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Tangential Velocity HWRF POM oper')


######################################################################################
kw = dict(levels=np.arange(-24,25,2))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contour(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_POM_exp_mean_dr[t,:,:],colors='k',**kw)
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_POM_exp_mean_dr[t,:,:],vmin=-24,vmax=24,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Tangential Velocity HWRF POM exp')

######################################################################################
kw = dict(levels=np.arange(-24,25,2))
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.contour(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_HYCOM_exp_mean_dr[t,:,:],colors='k',**kw)
plt.contourf(np.arange(0,561,delta_r),press_HWRF_POM_exp[t,:],vt_HWRF_HYCOM_exp_mean_dr[t,:,:],vmin=-24,vmax=24,cmap='Spectral_r',**kw)
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.xticks(np.arange(0,551,50))
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('Tangential Velocity HWRF HYCOM exp')

######################################################################################
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(vt_HWRF_POM_oper_mean_dr[t,:,0],press_HWRF_POM_oper[t,:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(vt_HWRF_POM_exp_mean_dr[t,:,0],press_HWRF_POM_exp[t,:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(vt_HWRF_HYCOM_exp_mean_dr[t,:,0],press_HWRF_HYCOM_exp[t,:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
plt.title('Mean (0-'+str(delta_r)+' km) Tangential Velocity',fontsize=16)
plt.ylabel('Pressure (hPa)',fontsize=14)
plt.xlabel('Tangential Velocity ($m/s$)',fontsize=14)
plt.legend()
plt.yscale('log')
plt.ylim([50,1000])
plt.yticks(np.arange(100,1001,100))
ax1.set_yticklabels(['100','200','300','400','500','600','700','800','900','1000'])
plt.gca().invert_yaxis()

######################################################################################
