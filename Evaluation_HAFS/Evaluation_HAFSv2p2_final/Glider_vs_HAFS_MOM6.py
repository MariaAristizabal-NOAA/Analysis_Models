#%% User input

# forecasting cycle to be used
# forecasting cycle to be used

# Erin
cycle = '2025081500'
storm_num = '05'
basin = 'al'
storm_id = '05l'
storm_name= 'erin'
models = ['hfsa','hfsa','hfsb','hfsb']

exp_names = ['HFSA_oper','HFSAv2p2_final_HAF3','HFSB_oper','HFSBv2p2_final_HBFZ']
exp_labels = ['HA21','HA22','HB21','HB22']
exp_colors = ['darkviolet','cyan','lime','orange']
ocean = ['mom6','mom6','hycom','mom6']

lon_lim = [-80,-50.0]
lat_lim = [10.0,40.0]

scratch_folder = '/scratch3/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch4/HFIP/hwrfv3/noscrub/input/abdeck/'
folder_myutils= '/home/Maria.Aristizabal/Maria_Utils/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

url_glider = scratch_folder + 'Data/Gliders/2025/sg668-20250625T0000.nc'
#url_glider = scratch_folder + 'Data/Gliders/2025/ng783-20250611T0000.nc'
#url_glider = scratch_folder + 'Data/Gliders/2025/echo-20250810T0000.nc'

################################################################################
import sys
import glob
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_glider_transect_from_HAFS_OCEAN,\
                            figure_transect_time_vs_depth, glider_data_vector_to_array,\
                            grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

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
#%% Read glider data
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
#target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
#lon_glider = target_lonG
#lat_glider = target_latG
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
target_temp_10m = np.empty((len(folder_exps),43))
target_temp_10m[:] = np.nan
target_salt_10m = np.empty((len(folder_exps),43))
target_salt_10m[:] = np.nan
target_time = np.empty((len(folder_exps),43))
target_time[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    #%% Get list files
    files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*storm.atm*.grb2')))

    hafs_ocean = xr.open_dataset(files_hafs_ocean[0])
    lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
    lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
    depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])

    #%% Get storm track from trak atcf files
    file_track = folder+storm_id+'.'+cycle+'.'+models[i]+'.trak.atcfunix'

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

    #%% Read HAFS/MOM6 time
    time_ocean = []
    timestamp_ocean = []
    for n,file in enumerate(files_hafs_ocean):
        print(file)
        OCEAN = xr.open_dataset(file)
        t = OCEAN['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_ocean.append(mdates.num2date(timestamp))
        timestamp_ocean.append(timestamp)

    time_ocean = np.asarray(time_ocean)
    timestamp_ocean = np.asarray(timestamp_ocean)

    #%% Retrieve glider transect from HAFS_HYCOM
    ncfiles = files_hafs_ocean
    lon = lon_hafs_ocean
    lat = lat_hafs_ocean
    depth = depth_hafs_ocean
    time_name = 'time'
    temp_name = 'temp'
    salt_name = 'so'

    if np.min(lon) < 0:
        lon_glid = longg
    else: 
        lon_glid = lon_glider
    lat_glid = latgg

    target_t, target_temp_hafs_ocean, target_salt_hafs_ocean = \
    get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

    timestamp = mdates.date2num(target_t)
    
    max_depth = 200
    kw_temp = dict(levels = np.arange(11,32,1))
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_temp_hafs_ocean,date_ini,date_end,max_depth,kw_temp,'Spectral_r','Degrees')
    plt.title(exp_labels[i],fontsize=16)

    max_depth = 200
    kw_salt = dict(levels = np.arange(34,37.5,0.3))
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_salt_hafs_ocean,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
    plt.title(exp_labels[i],fontsize=16)

    # Temp at 10 meters depth
    okd = np.where(depth_hafs_ocean <= 10)[0]
    target_temp_10m[i,:] = target_temp_hafs_ocean[okd[-1],:]
    target_salt_10m[i,:] = target_salt_hafs_ocean[okd[-1],:]
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
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')

plt.xlim(lon_lim)
plt.ylim(lat_lim)

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
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
plt.plot(timegg,tempg10,'o-',color='blue',label=dataset_id.split('-')[0],markeredgecolor='k')
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
plt.plot(timegg,saltg10,'o-',color='blue',label=dataset_id.split('-')[0],markeredgecolor='k')
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
figure_transect_time_vs_depth(timegg,-depthg_gridded,tempg_gridded,date_ini,date_end,max_depth,kw_temp,'Spectral_r','Degress')
plt.title(dataset_id,fontsize=16)

max_depth = 200
kw_salt = dict(levels = np.arange(34,37.5,0.3))
figure_transect_time_vs_depth(timegg,-depthg_gridded,saltg_gridded,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
plt.title(dataset_id,fontsize=16)

