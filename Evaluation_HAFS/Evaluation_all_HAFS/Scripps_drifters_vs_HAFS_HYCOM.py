#%% User input
# forecasting cycle to be used

# Ian
cycle = '2022092718'
storm_num = '09'
basin = 'al'
url_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2022.nc'

# Ian
#cycle = '2022092906'
#storm_num = '09'
#basin = 'al'

# Fiona
#cycle = '2022091806'
#storm_num = '07'
#basin = 'al'

exp_names = ['hafsv1_fnl_hfsa','hafsv1_baseline']
exp_labels = ['HFSA','HAFS_baseline']
exp_colors = ['c','orange']

lon_lim = [-75,-55.0]
lat_lim = [20.0,50.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

#folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',
#               scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

folder_exps = [scratch_folder2 + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',scratch_folder2 + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = abdeck_folder + 'aid/a' + basin + storm_num + cycle[0:4] + '.dat'

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
#import seawater as sw

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

    # Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
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

# Get drifter variables within lat and lon range
#oklat = np.logical_and(latDr >= np.nanmin(lat_best_track) - 4, latDr <= np.nanmax(lat_best_track) + 4)
#lonDD = lonDr[oklat]
#oklon = np.logical_and(lonDD >= np.nanmin(lon_best_track), lonDD <= np.nanmax(lon_best_track))

oklat = np.logical_and(latDr >= np.nanmin(lat_forec_track[0,:]) - 0.5, latDr <= np.nanmax(lat_forec_track[0,:]) + 0.5)
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= np.nanmin(lon_forec_track[0,:]) - 0.5, lonDD <= np.nanmax(lon_forec_track[0,:]) + 0.5)

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

# Find the different drifter within lat, lon and time window
codes = np.unique(platform_codeD)

########################################################################

target_tmp_surf = np.empty((len(folder_exps),22))
target_tmp_surf[:] = np.nan
target_salt_surf = np.empty((len(folder_exps),22))
target_salt_surf[:] = np.nan
target_press_surf = np.empty((len(folder_exps),43))
target_press_surf[:] = np.nan

# Loop through nearby drifter
#for code in codes:
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

    # Loop the experiments
    for i,folder in enumerate(folder_exps):
        print(folder)
        # Get list files
        files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfsa.parent.atm.*.grb2')))
        if len(files_hafs_fv3) == 0:
            files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs.grid01.*.grb2')))
        
        # Reading HAFS/HYCOM grid
        hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
        lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
        lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
        depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

        #%% Reading HAFS/FV3 grid
        fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
        lon_hafs_fv3 = np.asarray(fv3.lon_0)
        lat_hafs_fv3 = np.asarray(fv3.lat_0)

        '''
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
        '''
        #################################################################
        # Retrieve HAFS_HYCOM temp. following saildrone trajectory
    
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

        target_timedhyc, target_tmp_surf[i,:] = get_var_from_model_following_trajectory(files_model,'temperature',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

        _, target_salt_surf[i,:] = get_var_from_model_following_trajectory(files_model,'salinity',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

        #################################################################
        # Retrieve HAFS_atm press following drifters trajectory

        files_model = files_hafs_fv3
        time_name = 'time'
        lat_name = 'latitude'
        lon_name = 'longitude'
        depth_level = 0
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

        target_timedfv3, target_press_surf[i,:] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

        #_, target_tmp_surf[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

        #_, target_tmp_2mabove[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L103_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

        #_, target_rh_2mabove[i,:] = get_var_from_model_following_trajectory(files_model,'RH_P0_L103_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    # Figure track
    lev = np.arange(-9000,9100,100)
    fig,ax = plt.subplots()
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.plot(lond, latd,'.',color='purple',label='Drifter id ' + str(code))
    #plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
    plt.legend(loc='upper right')
    plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
    plt.axis('scaled')
    plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
    plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])

    # Figure sea level pressure
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,slpd,'.-',color='blue',label='Drifter '+ str(code))
    for i in np.arange(len(exp_names)):
        plt.plot(target_timedfv3,target_press_surf[i,:]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend()
    plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
    plt.ylabel('(hPa)',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)

    # Figure SST
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,sstd,'.-',color='blue',label='Drifter '+ str(code))
    for i in np.arange(len(exp_names)):
        plt.plot(target_timedhyc,target_tmp_surf[i,:],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend()
    plt.title('Sea Surface Temperature Cycle '+ cycle,fontsize=18)
    plt.ylabel('($^oC$)',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)

    # Figure SSS
    '''
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,sssd,'.-',color='blue',label='Drifter '+ str(code))
    for i in np.arange(len(exp_names)):
        plt.plot(target_timedhyc,target_salt_surf[i,:],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.legend()
    plt.title('Sea Surface Salinity Cycle '+ cycle,fontsize=18)
    plt.ylabel(' ',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    '''
########################################################################
#%% Figure track
'''
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonD, latD,'.',color='purple',label='Scrips Drifter id ' + str(platform_codeD[0]))
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
'''

########################################################################
