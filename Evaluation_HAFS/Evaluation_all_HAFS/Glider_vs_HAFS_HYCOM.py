#%% User input

# forecasting cycle to be used

# Fiona
cycle = '2022091806'
storm_num = '07'
basin = 'al'
url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2022/SG610-20220621T1130.nc' 

exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start','hafsv0p3a_2022rt_natl']
exp_labels = ['HF3C','HF3A']
exp_colors = ['orange','c']

#exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start']
#exp_labels = ['HF3C']
#exp_colors = ['orange']

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',
               scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

#folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = scratch_folder + 'adeck/a' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'
#folder_glider_comp = '/home/Maria.Aristizabal/Repos/glider_model_comparisons_Python/'

folder_fig = scratch_folder + 'figures/'

#url_sg633 = scratch_folder + 'Data/Glider_data/SG663-20210615T1202_d95b_2e68_343c.nc'
#url_bio_jack = scratch_folder + 'Data/Glider_data/bios_jack-20210709T1945_c92d_71b9_1bdc.nc'
#url_ng645 = scratch_folder + 'Data/Glider_data/ng645-20210613T0000_3bdf_c72d_67da.nc'

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
#import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_HYCOM,\
                            figure_transect_temp, glider_data_vector_to_array,\
                            grid_glider_data


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
#%% Read glider data
'''
# Elsa
if cycle == '2021070100':
    url_glider = url_sg633

# Larry and Sam
if cycle == '2021090800' or cycle == '2021100100':
    url_glider = url_bios_jack

# Ida
if cycle == '2021082800' or cycle == '2021082706':
    url_glider = url_ng645
'''

gdata = xr.open_dataset(url_glider)#,decode_times=False)

dataset_id = gdata.id.split('_')[0]
temperature = np.asarray(gdata.variables['temperature'][:])
salinity = np.asarray(gdata.variables['salinity'][:])
density = np.asarray(gdata.variables['density'][:])
latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
depth = np.asarray(gdata.depth)

time = np.asarray(gdata.time)
time = np.asarray(mdates.num2date(mdates.date2num(time)))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
temperat =  temperature[oktimeg].T
salinit =  salinity[oktimeg].T
densit =  density[oktimeg].T
latitud = latitude[oktimeg]
longitud = longitude[oktimeg]
depthh = depth[oktimeg].T
timee = time[oktimeg]

contour_plot = 'no' # default value is 'yes'
delta_z = 0.2     # default value is 0.3

depthg, timeg, tempg, latg, long = glider_data_vector_to_array(depthh,timee,temperat,latitud,longitud)
depthg, timeg, saltg, latg, long = glider_data_vector_to_array(depthh,timee,salinit,latitud,longitud)

ok = np.where(np.isfinite(timeg[0,:]))[0]
timegg = timeg[0,ok]
tempgg = tempg[:,ok]
saltgg = saltg[:,ok]
depthgg = depthg[:,ok]
longg = long[0,ok]
latgg = latg[0,ok]

tempg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(tempgg,timegg,depthgg,delta_z)
saltg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(saltgg,timegg,depthgg,delta_z)

#tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
tstamp_glider = timeg_gridded

# Conversion from glider longitude and latitude to HYCOM convention
target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
lon_glider = target_lonG
lat_glider = target_latG
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
rmw_track = np.empty((len(folder_exps),43))
rmw_track[:] = np.nan
target_temp_10m = np.empty((len(folder_exps),22))
target_temp_10m[:] = np.nan
target_salt_10m = np.empty((len(folder_exps),22))
target_salt_10m[:] = np.nan
target_time = np.empty((len(folder_exps),22))
target_time[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    #if len(files_hafs_fv3) == 0:
    #    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], rmw_track[i,0:okn] = get_storm_track_and_int(file_track,storm_num)

    #%% Read HAFS/Fv3 time
    '''
    time_fv3 = []
    for n,file in enumerate(files_hafs_fv3):
        print(file)
        #FV3 = xr.open_dataset(file)
        #t = FV3.variables['time'][:]
        #timestamp = mdates.date2num(t)[0]
        #time_fv3.append(mdates.num2date(timestamp))

        FV3 = xr.open_dataset(file,engine="pynio")
        t0 = FV3.variables['TMP_P0_L1_GLL0'].attrs['initial_time']
        dt = FV3.variables['TMP_P0_L1_GLL0'].attrs['forecast_time'][0]
        time_fv3.append(datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt)))

    time_fv3 = np.asarray(time_fv3)
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

    #%% Retrieve glider transect from HAFS_HYCOM
    ncfiles = files_hafs_hycom
    lon = lon_hafs_hycom
    lat = lat_hafs_hycom
    depth = depth_hafs_hycom

    if np.min(lon) < 0:
        lon_glid = longg
    else: 
        lon_glid = lon_glider
    lat_glid = lat_glider

    target_t, target_temp_hafs_hycom, target_salt_hafs_hycom = \
    get_glider_transect_from_HAFS_HYCOM(ncfiles,lon,lat,depth,tstamp_glider,lon_glid,lat_glid)

    timestamp = mdates.date2num(target_t)
    
    max_depth = 200
    kw_temp = dict(levels = np.arange(11,32,1))
    figure_transect_temp(target_t,-depth,target_temp_hafs_hycom,date_ini,date_end,max_depth,kw_temp,'Spectral_r')
    plt.title(exp_labels[i],fontsize=16)

    max_depth = 200
    kw_salt = dict(levels = np.arange(34,37.5,0.3))
    figure_transect_temp(target_t,-depth,target_salt_hafs_hycom,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r')
    plt.title(exp_labels[i],fontsize=16)

    # Temp at 10 meters depth
    okd = np.where(depth_hafs_hycom <= 10)[0]
    target_temp_10m[i,:] = target_temp_hafs_hycom[okd[-1],:]
    target_salt_10m[i,:] = target_salt_hafs_hycom[okd[-1],:]
    target_time[i,:] = timestamp 
  
##################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='purple',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='purple',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(long)-1,np.nanmax(long)+1])
plt.ylim([np.nanmin(latg)-1,np.nanmax(latg)+1])

##################################################################
#%% Figure time series temp at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
tempg10 = tempg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,tempg10,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)): 
    plt.plot(target_time[i,:],target_temp_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Temperature ($^oC$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Temperature at 10 m depth',fontsize=16)
plt.grid(True)

#################################################################################
#%% Figure time series salt at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
saltg10 = saltg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,saltg10,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)):
    plt.plot(target_time[i,:],target_salt_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Salt',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Salinity at 10 m depth',fontsize=16)
plt.grid(True)

#################################################################################
#%% Figures transects

#%% Glider
max_depth = 200
kw_temp = dict(levels = np.arange(11,32,1))
figure_transect_temp(timegg,-depthg_gridded,tempg_gridded,date_ini,date_end,max_depth,kw_temp,'Spectral_r')
plt.title(dataset_id,fontsize=16)

max_depth = 200
kw_salt = dict(levels = np.arange(34,37.5,0.3))
figure_transect_temp(timegg,-depthg_gridded,saltg_gridded,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r')
plt.title(dataset_id,fontsize=16)

