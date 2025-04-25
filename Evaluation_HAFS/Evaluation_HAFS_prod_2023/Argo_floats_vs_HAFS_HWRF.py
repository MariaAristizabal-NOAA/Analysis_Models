#%% User input

# forecasting cycle to be used
# Otis
cycle = '2023102218'
storm_num = '18'
basin = 'ep'
storm_id = '18e'
storm_name = 'otis'
url_argo = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Argo_floats/2023/argo_floats_Hurr_Otis_2023_qc.nc'

exp_names = ['HFSA_oper']
exp_labels = ['HFSA_oper_HYCOM_kpp']
exp_colors = ['darkviolet']

lon_lim = [-110,-80]
lat_lim = [0,20]

home_folder = '/home/Maria.Aristizabal/'
home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder1 + exp_names[0] + '/' + cycle + '/' + storm_num + basin[0] + '/']

#bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2021_sub_ice_topo.nc'

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
import seawater as sw

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
bath_lat = np.asarray(ncbath.variables['lat'][:])
bath_lon = np.asarray(ncbath.variables['lon'][:])
bath_elev = np.asarray(ncbath.variables['elevation'][:])

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

okt = time_best_track >= tini

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
#%% Read Argo Float data

url = url_argo

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
platform_number = np.asarray(gdata.platform_number)
depth = np.asarray(gdata.pres)
depth_qc = np.asarray(gdata.pres_qc)
temp = np.asarray(gdata.temp)
temp_qc = np.asarray(gdata.temp_qc)
salt = np.asarray(gdata.psal)
salt_qc = np.asarray(gdata.psal_qc)

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
platform_numberDr = platform_number[oktimeg]
depthDr = depth[oktimeg]
depthDr_qc = depth_qc[oktimeg]
tempDr = temp[oktimeg]
tempDr_qc = temp_qc[oktimeg]
saltDr = salt[oktimeg]
saltDr_qc = salt_qc[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
platform_numberD = platform_numberDr[oklat][oklon]
depthD = depthDr[oklat][oklon]
depthD_qc = depthDr_qc[oklat][oklon]
tempD = tempDr[oklat][oklon]
tempD_qc = tempDr_qc[oklat][oklon]
saltD = saltDr[oklat][oklon]
saltD_qc = saltDr_qc[oklat][oklon]

platf_numbers = np.unique(platform_numberD)

########################################################################

target_timeD_ocean = []
target_temp = np.empty((len(folder_exps),43))
target_temp[:] = np.nan
target_salt = np.empty((len(folder_exps),43))
target_salt[:] = np.nan
target_depth = np.empty((len(folder_exps),43))
target_depth[:] = np.nan

# Loop through nearby drifter
for number in platf_numbers:
    print(number)
    okcode = platform_numberD == number
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
    platf_numd = platform_numberD[okcode]
    depthd = depthD[okcode]
    depthd_qc = depthD_qc[okcode]
    tempd = tempD[okcode]
    tempd_qc = tempD_qc[okcode]
    saltd = saltD[okcode]
    saltd_qc = saltD_qc[okcode]

    densd = sw.dens(saltd,tempd,depthd)
    ohcd = OHC_from_profile(depthd,tempd,densd)

    # Loop the experiments
    for i,folder in enumerate(folder_exps):
        print(folder)

        if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSB' or exp_names[i] == 'HFSAv1p1_HYCOM':
            #%% Get list files
            files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))

            # Reading HAFS/HYCOM grid
            hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],engine="pynio")
            lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
            lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
            depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

            # Reading HAFS/HYCOM time
            time_hafs_hycom = []
            timestamp_hafs_hycom = []
            for file in files_hafs_hycom:
                hafs_hycom = xr.open_dataset(file,engine="pynio")
                time_hafs_hycom.append(np.asarray(hafs_hycom['MT'][:])[0])
            
            timestamp_hafs_hycom = mdates.date2num(time_hafs_hycom)

            # Conversion from Argo float longitude and latitude to HYCOM convention
            lond_hyc, latd_hyc = geo_coord_to_HYCOM_coord(lond,latd)
            lond_hyc = np.asarray(lond_hyc)
            latd_hyc = np.asarray(latd_hyc)

            timestamp_obss = timestampd

            if np.min(lon_hafs_hycom) < 0:
                lon_obss = lond
            else: 
                lon_obss = lond_hyc
            lat_obss = latd_hyc

            oklo = np.isfinite(lon_obss)
            lon_obs = np.unique(lon_obss[oklo])
            lat_obs = np.unique(lat_obss[oklo])
            timestamp_obs = np.unique(timestamp_obss[oklo])

            # Retrieve HAFS_HYCOM temp. for Argo Float time and position
            oktt = np.round(np.interp(timestamp_obs,timestamp_hafs_hycom,np.arange(len(timestamp_hafs_hycom)))).astype(int)[0]
            oklat = np.where(lat_hafs_hycom >= lat_obs[0])[0][0]
            oklon = np.where(lon_hafs_hycom >= lon_obs[0])[0][0]
            hafs_hycom = xr.open_dataset(files_hafs_hycom[oktt],engine="pynio")
            target_temp = np.asarray(hafs_hycom['temperature'][0,:,oklat,oklon])
            target_salt = np.asarray(hafs_hycom['salinity'][0,:,oklat,oklon])

            dens_prof = sw.dens(target_salt,target_temp,depth_hafs_hycom)
            ohc_hafs_hycom = OHC_from_profile(depth_hafs_hycom,target_temp,dens_prof)

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
    

            ############################################################################
        if exp_names[i] == 'HFSAv1p1_MOM6_epbl' or exp_names[i] == 'HFSAv1p1_MOM6_kpp':
        
            #%% Get list files
            files_hafs_mom6 = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))

            #%% Reading MOM6 grid
            hafs_mom6_grid = xr.open_dataset(files_hafs_mom6[0],decode_times=False)
            lon_hafs_mom6 = np.asarray(hafs_mom6_grid['xh'][:])
            lat_hafs_mom6 = np.asarray(hafs_mom6_grid['yh'][:])
            depth_hafs_mom6 = np.asarray(hafs_mom6_grid['z_l'][:])

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

            #################################################################################
    
    # Figure track
    lev = np.arange(-9000,9100,100)
    fig,ax = plt.subplots(figsize=(12,4))
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.plot(np.unique(lond),np.unique(latd),'*',markersize=10,color='orangered',label='Argo Float id ' + str(number))
    plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
    plt.axis('scaled')
    plt.xlim(lon_lim)
    plt.ylim(lat_lim)
    
    # Figure temperature profile
    fig,ax = plt.subplots(figsize=(6,9))
    plt.plot(tempd,-depthd,'o-',color='orangered',label='Argo Float id '+ str(number),markersize=5,markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        plt.plot(target_temp,-depth_hafs_hycom,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    #plt.legend(loc='lower right',bbox_to_anchor=[1.7,1.0])
    plt.legend(loc='lower right')
    plt.title('Temperature Cycle '+ cycle,fontsize=18)
    plt.xlabel('($^oC$)',fontsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.ylim([-300,0])
    plt.xlim([5,32])
    plt.grid()
    plt.text(0.45,0.4,'OHC Argo = '+str(np.round(ohcd,1))+' $kJ/cm^2$',fontsize=14,transform=ax.transAxes)
    plt.text(0.45,0.35,'OHC HFSA = '+str(np.round(ohc_hafs_hycom,1))+' $kJ/cm^2$',fontsize=14,transform=ax.transAxes)
    plt.text(0.45,0.3,'Mean pres_qc = '+str(np.round(np.nanmean(depthd_qc.astype(int)),1)),fontsize=14,transform=ax.transAxes)
    plt.text(0.45,0.25,'Mean temp_qc = '+str(np.round(np.nanmean(tempd_qc.astype(int)),1)),fontsize=14,transform=ax.transAxes)

    # Figure salinity profile
    fig,ax = plt.subplots(figsize=(6,9))
    plt.plot(saltd,-depthd,'o-',color='orangered',label='Argo Float id '+ str(number),markersize=5,markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        plt.plot(target_salt,-depth_hafs_hycom,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    #plt.legend(loc='lower right',bbox_to_anchor=[1.7,1.0])
    plt.legend(loc='lower left')
    plt.title('Salinity Cycle '+ cycle,fontsize=18)
    plt.xlabel(' ',fontsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.ylim([-300,0])
    plt.text(0.1,0.2,'Mean pres_qc = '+str(np.round(np.nanmean(depthd_qc.astype(int)),1)),fontsize=14,transform=ax.transAxes)
    plt.text(0.1,0.14,'Mean salt_qc = '+str(np.round(np.nanmean(saltd_qc.astype(int)),1)),fontsize=14,transform=ax.transAxes)
    #plt.xlim([5,32])
    plt.grid()

    ########################################################################

# Figure track
lev = np.arange(-9000,9100,100)
fig,ax = plt.subplots(figsize=(12,4))
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonD, latD,'*',color='orangered',label='Argo Floats')
plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim(lon_lim)
plt.ylim(lat_lim)

