#%% User input
# forecasting cycle to be used

# Lee
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'
url_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2023.nc'

exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_MOM6_epbl']
exp_labels = ['HFSA-HYCOM (oper)','HWRF-POM (oper)','HFSA-MOM6 (exp)']
exp_colors = ['darkviolet','pink','forestgreen']

#exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_HYCOM','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_labels = ['HFSA_oper_HYCOM_kpp','HWRF_POM_MY2.5','HFSAv1p1_HYCOM_kpp','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_colors = ['darkviolet','pink','forestgreen','cyan','royalblue']

lon_lim = [-83,-30]
lat_lim = [10,45]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
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
                            get_var_from_model_following_trajectory


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#########################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder1 + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

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
# Loop the experiments to obtain forecasted track

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    #%% Get storm track from trak atcf files
    if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSAv1p1_HYCOM' or exp_names[i] == 'HFSAv1p1_MOM6_epbl' or exp_names[i] == 'HFSAv1p1_MOM6_kpp':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if exp_names[i] == 'HFSB':
        file_track = folder + + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
    if exp_names[i] == 'HWRF_2023':
        file_track = folder + storm_name + storm_id + '.' + cycle + '.trak.hwrf.atcfunix'

    print(file_track)
    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

#################################################################################
#%% Read drifter data

url = url_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
platform_code = np.asarray(gdata.platform_code)
wind_speed = np.asarray(gdata.windspd)
sea_level_press = np.asarray(gdata.slp)
sea_surface_temp = np.asarray(gdata.sst)
sea_surface_salt = np.asarray(gdata.salinity)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(times)
times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(times) >= mdates.date2num(tini),\
                         mdates.date2num(times) <= mdates.date2num(tend))

# Fields within time window
timeDr = times[oktimeg]
timestampDr = timestamps[oktimeg]
latDr = latitude[oktimeg]
lonDr = longitude[oktimeg]
platform_codeDr = platform_code[oktimeg]
windspdDr = wind_speed[oktimeg]
slpDr = sea_level_press[oktimeg]
sstDr = sea_surface_temp[oktimeg]
sssDr = sea_surface_salt[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
platform_codeD = platform_codeDr[oklat][oklon]
windspdD = windspdDr[oklat][oklon]
slpD = slpDr[oklat][oklon]
sstD = sstDr[oklat][oklon]
sssD = sssDr[oklat][oklon]

codes = np.unique(platform_codeD)

########################################################################

target_timeD_ocean = []
target_timeD_atm = []
target_temp = np.empty((len(folder_exps),43))
target_temp[:] = np.nan
target_salt = np.empty((len(folder_exps),43))
target_salt[:] = np.nan
target_press_surf = np.empty((len(folder_exps),43))
target_press_surf[:] = np.nan
#target_timed_atm = [[0 for i in range(43)] for j in range(len(folder_exps))]
#target_timed_ocean = [[0 for i in range(43)] for j in range(len(folder_exps))]

# Loop through nearby drifter
#for code in codes:
#for code in [1701715.,1301700.,1301783.]:
#for code in [3201837.]:
for code in [1301700.]:
    print(code)
    okcode = platform_codeD == code
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
    platform_coded = platform_codeD[okcode]
    windspdd = windspdD[okcode]
    slpd = slpD[okcode]
    sstdd = sstD[okcode]
    sssd = sssD[okcode]

    sstd = np.empty((len(sstdd)))
    sstd[:] = np.nan
    for i,sst in enumerate(sstdd):
        if sst == 0:
            sstd[i] = np.nan
        else:
            sstd[i] = float(sst)

    # Loop the experiments
    for i,folder in enumerate(folder_exps):
        print(folder)

        if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSB' or exp_names[i] == 'HFSAv1p1_HYCOM':
            #%% Get list files
            files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
            files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))

            # Reading HAFS/HYCOM grid
            hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
            lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
            lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
            depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

            #%% Reading HAFS/FV3 grid
            fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
            lon_hafs_fv3 = np.asarray(fv3.lon_0)
            lat_hafs_fv3 = np.asarray(fv3.lat_0)

            # Retrieve HAFS_HYCOM temp. following drifter trajectory
    
            # Conversion from glider longitude and latitude to HYCOM convention
            lond_hyc, latd_hyc = geo_coord_to_HYCOM_coord(lond,latd)
            lond_hyc = np.asarray(lond_hyc)
            latd_hyc = np.asarray(latd_hyc)

            files_model = files_hafs_hycom
            time_name = 'MT'
            lat_name = 'Latitude'
            lon_name = 'Longitude'
            depth_level = 0
            timestamp_obss = timestampd
            kwargs = dict(depth_level = 0)
            if np.min(lon_hafs_hycom) < 0:
                lon_obss = lond
            else: 
                lon_obss = lond_hyc
            lat_obss = latd_hyc

            oklo = np.isfinite(lon_obss)
            lon_obs = lon_obss[oklo]
            lat_obs = lat_obss[oklo]
            timestamp_obs = timestamp_obss[oklo]

            target_timedhyc, target_temp[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'temperature',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

            _, target_salt[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'salinity',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

            target_timeD_ocean.append(target_timedhyc)

            #################################################################
            # Retrieve HAFS_atm press following drifters trajectory

            files_model = files_hafs_fv3
            time_name = 'time'
            lat_name = 'latitude'
            lon_name = 'longitude'
            #depth_level = 0
            if np.min(lon_hafs_fv3) < 0:
                lon_obss = lond
            else:
                lon_obss = lond + 360
            lat_obss = latd
            timestamp_obss = timestampd
    
            oklo = np.isfinite(lon_obss)
            lon_obs = lon_obss[oklo]
            lat_obs = lat_obss[oklo]
            timestamp_obs = timestamp_obss[oklo]
    
            target_timedatm, target_press_surf[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)
            #_, target_tmp_surf[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)
    
            target_timeD_atm.append(target_timedatm)

        if exp_names[i] == 'HWRF_2023':
            #%% Get list files
            files_pom = sorted(glob.glob(os.path.join(folder,'*pom*00*.nc')))
    
            #%% Reading POM grid
            print('Retrieving coordinates from POM')
            grid_file = glob.glob(os.path.join(folder,'*pom.grid.nc'))[0]
            pom_grid = xr.open_dataset(grid_file)
            lon_pom = np.asarray(pom_grid['east_e'][:])
            lat_pom = np.asarray( pom_grid['north_e'][:])
            zlevc = np.asarray(pom_grid['zz'][:])
            topoz = np.asarray(pom_grid['h'][:])
    
            #%% Read POM time
            time_pom = []
            for n,file in enumerate(files_pom):
                print(file)
                pom = xr.open_dataset(file)
                #ocean = xr.open_dataset(files_ocean[0],decode_times=False)
                t = pom.variables['time'][:]
                timestamp = mdates.date2num(t)[0]
                time_pom.append(mdates.num2date(timestamp))
    
            time_pom = np.asarray(time_pom)
    
            #%% Get list HWRF files
            files_hwrf = sorted(glob.glob(os.path.join(folder,'*hwrfprs.synoptic.0p125*.grb2')))
    
            #%% Reading hwrf grid
            hwrf = xr.open_dataset(files_hwrf[0],engine="pynio")
            lon_hwrf = np.asarray(hwrf.lon_0)
            lat_hwrf = np.asarray(hwrf.lat_0)
    
            # Read track file
            okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
            lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
    
            ############################################################################
            #%% Retrieve POM temp. following saildrone trajectory
    
            files_model = files_pom
            time_name = 'time'
            lat_name = 'north_e'
            lon_name = 'east_e'
            timestamp_obss = timestampd
            kwargs = dict(depth_level = 0)
    
            lon_obss = lond
            lat_obss = latd
            oklo = np.isfinite(lon_obss)
            lon_obs = lon_obss[oklo]
            lat_obs = lat_obss[oklo]
            timestamp_obs = timestamp_obss[oklo]
    
            target_timedpom, target_temp[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'t',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)
    
            _, target_salt[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'s',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)
    
            target_timeD_ocean.append(target_timedpom)

            ############################################################################
            #%% Retrieve HWRF press following saildrone trajectory
    
            files_model = files_hwrf
            time_name = 'time'
            #lat_name = 'latitude'
            #lon_name = 'longitude'
            depth_level = 0
            if np.min(lon_hwrf) < 0:
                lon_obss = lond
            else:
                lon_obss = lond + 360
            lat_obss = latd
            timestamp_obss = timestampd
    
            oklo = np.isfinite(lon_obss)
            lon_obs = lon_obss[oklo]
            lat_obs = lat_obss[oklo]
            timestamp_obs = timestamp_obss[oklo]
    
            target_timedatm, target_press_surf[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)
            #_, target_tmp_surf[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

            target_timeD_atm.append(target_timedatm)

        if exp_names[i] == 'HFSAv1p1_MOM6_epbl' or exp_names[i] == 'HFSAv1p1_MOM6_kpp':
        
            #%% Get list files
            files_hafs_mom6 = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))

            #%% Reading MOM6 grid
            hafs_mom6_grid = xr.open_dataset(files_hafs_mom6[0],decode_times=False)
            lon_hafs_mom6 = np.asarray(hafs_mom6_grid['xh'][:])
            lat_hafs_mom6 = np.asarray(hafs_mom6_grid['yh'][:])
            depth_hafs_mom6 = np.asarray(hafs_mom6_grid['z_l'][:])

            #%% Get list files
            #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.storm.atm.*.grb2')))
            files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))

            #%% Reading FV3 grid
            fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
            lon_hafs_fv3 = np.asarray(fv3.lon_0)
            lat_hafs_fv3 = np.asarray(fv3.lat_0)

            #%% Read HAFS/HYCOM time
            time_mom6 = []
            timestamp_mom6 = []
            for n,file in enumerate(files_hafs_mom6):
                print(file)
                MOM6 = xr.open_dataset(file)
                t = MOM6['time'][:]
                timestamp = mdates.date2num(t)[0]
                time_mom6.append(mdates.num2date(timestamp))
                timestamp_mom6.append(timestamp)

            time_mom6 = np.asarray(time_mom6)
            timestamp_mom6 = np.asarray(timestamp_mom6)

            # Read track file
            okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
            lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

            #################################################################################
            #%% Retrieve HAFS_HYCOM temp. following saildrone trajectory

            files_model = files_hafs_mom6
            time_name = 'time'
            lat_name = 'yh'
            lon_name = 'xh'
            depth_level = 0
            timestamp_obsd = timestampd
            kwargs = dict(depth_level = 0)
            lon_obsd = lond
            lat_obsd = latd

            oklo = np.isfinite(lon_obsd)
            lon_obs = lon_obsd[oklo]
            lat_obs = lat_obsd[oklo]
            timestamp_obs = timestamp_obsd[oklo]

            target_timeDmom6, target_temp[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'temp',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

            _, target_salt[i,0:len(files_model)] = get_var_from_model_following_trajectory(files_model,'so',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

            target_timeD_ocean.append(target_timeDmom6)

            #################################################################################
            #%% Retrieve HAFS_atm press following saildrone trajectory

            files_model = files_hafs_fv3
            time_name = 'time'
            lat_name = 'latitude'
            lon_name = 'longitude'
            #depth_level = 0
            if np.min(lon_hafs_fv3) < 0:
                lon_obsd = lond
            else:
                lon_obsd = lond + 360
            lat_obsd = latd
            timestamp_obsd = timestampd

            oklo = np.isfinite(lon_obsd)
            lon_obs = lon_obsd[oklo]
            lat_obs = lat_obsd[oklo]
            timestamp_obs = timestamp_obsd[oklo]
 
            #target_timedatm, target_press_surf[i,:] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)
            target_timedatm, target_press_surf[i,:] = get_var_from_model_following_trajectory(files_model,'PRMSL_P0_L101_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

            target_timeD_atm.append(target_timedatm)
    
    # Figure track
    lev = np.arange(-9000,9100,100)
    fig,ax = plt.subplots(figsize=(12,4))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.plot(lond, latd,'.',color='orangered',label='Drifter id ' + str(code))
    plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
    plt.axis('scaled')
    plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
    plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])

    # Figure track
    lev = np.arange(-9000,9100,100)
    fig,ax = plt.subplots(figsize=(12,4))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.plot(lond, latd,'.',color='orangered',label='Drifter id ' + str(code))
    plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
    plt.axis('scaled')
    plt.xlim([np.nanmin(lond)-3,np.nanmax(lond)+3])
    plt.ylim([np.nanmin(latd)-3,np.nanmax(latd)+3])
    
    # Figure sea level pressure
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,slpd,'o-',color='orangered',label='Drifter '+ str(code),markersize=5,markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        plt.plot(target_timeD_atm[i],target_press_surf[i,0:len(target_timeD_atm[i])]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
    plt.ylabel('(hPa)',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    
    # Figure SST
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,sstd,'o-',color='orangered',label='Drifter '+ str(code),markersize=5,markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        plt.plot(target_timeD_ocean[i],target_temp[i,0:len(target_timeD_ocean[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.title('Sea Surface Temperature Cycle '+ cycle,fontsize=18)
    plt.ylabel('($^oC$)',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    
    '''
    # Figure SSS
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,sssd,'.-',color='blue',label='Drifter '+ str(code))
    for i in np.arange(len(exp_names)):
        plt.plot(target_timeD_ocean[i],target_salt[i,0:len(target_timeD_ocean[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend()
    plt.title('Sea Surface Salinity Cycle '+ cycle,fontsize=18)
    plt.ylabel(' ',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    '''
    ########################################################################

'''
for code in [4102649.0,4102651.0,4102654.0]:
    print(code)
    okcode = platform_codeD == code
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
    platform_coded = platform_codeD[okcode]
    windspdd = windspdD[okcode]
    slpd = slpD[okcode]
    sstd = sstD[okcode]
    sssd = sssD[okcode]
'''

# Figure track
lev = np.arange(-9000,9100,100)
fig,ax = plt.subplots(figsize=(12,4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonD, latD,'.',color='orangered',label='Drifters')
plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])


