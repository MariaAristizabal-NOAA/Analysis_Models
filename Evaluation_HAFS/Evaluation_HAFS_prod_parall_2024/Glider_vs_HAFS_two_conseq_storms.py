#%% User input

# Helena and Milton
cycles = ['2024092500','2024100800']
storm_nums = ['09','14']
basins = ['al','al']
storm_ids = ['09l','14l']
storm_names = ['Helena','Milton']
# time eye passage closest to glider
teyes = ['2024092618','2024100918']

exp_names = ['HFSA_oper','hafs_20241220_v2p1a_baseline','hafs_20250210_v2p1a_ha30']
exp_labels = ['HFSA_oper','HAFSv2.1A baseline','HAFSv2.1A final']
exp_colors = ['purple','dodgerblue','#00c8c8']
hafs_ab = ['hfsa','hfsa','hfsa']
ocean = ['mom6','mom6','mom6']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']

'''
exp_names = ['HFSA_oper','HFSB_oper','HAFSv2p0p1a_2024rt','HAFSv2p0p1b_2024rt']
exp_labels = ['HFSA_oper','HFSB_oper','HFSA_para','HFSB_para']
exp_colors = ['purple','lime','dodgerblue','olive']
hafs_ab = ['hfsa','hfsb','hfsa','hfsb']
ocean = ['mom6','hycom','mom6','mom6']
'''

url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2024/usf-jaialai-20240920T0000.nc' 

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# folder utils for Hycom 
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_OCEAN,\
                            figure_transect_time_vs_depth,\
                            glider_data_vector_to_array,grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder[i] + exp_names[i])

################################################################################
#%% Time window
date_ini = cycles[0][0:4]+'/'+cycles[0][4:6]+'/'+cycles[0][6:8]+'/'+cycles[0][8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = datetime.strptime(cycles[1][0:4]+'/'+cycles[1][4:6]+'/'+cycles[1][6:8]+'/'+cycles[1][8:]+'/00/00','%Y/%m/%d/%H/%M/%S') + timedelta(hours=126)
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
lon_best_track = np.empty((len(cycles),50))
lon_best_track[:] = np.nan
lat_best_track = np.empty((len(cycles),50))
lat_best_track[:] = np.nan
#time_best_track = np.empty((len(cycles),50))
#time_best_track[:] = np.nan
time_best_track = []
for c,cycle in enumerate(cycles):
    best_track_file = abdeck_folder + 'btk/b' + basins[c] + storm_nums[c] + cycles[c][0:4] + '.dat'
    lon,_,_,_,_ = get_best_track_and_int(best_track_file)

    lon_best_track[c,0:len(lon)], lat_best_track[c,0:len(lon)], t_best_track, _, _ = get_best_track_and_int(best_track_file)

    #time_best_track[c,0:len(lon)] = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])
    time_best_track.append([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

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
delta_z = 1     # default value is 0.3

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

tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
tstamp_glider = timeg_gridded

# Conversion from glider longitude and latitude to HYCOM convention
target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
lon_glider = target_lonG
lat_glider = target_latG

#################################################################################
#%% Loop the experiments
nc = len(cycles)

lon_forec_track = np.empty((len(folder_exps),nc,len(time_fv3)))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),nc,len(time_fv3)))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),nc,len(time_fv3)))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),nc,len(time_fv3)))
int_track[:] = np.nan
target_time = np.empty((len(folder_exps),nc,len(time_fv3)))
target_time[:] = np.nan
target_temp_hafs_ocean = np.empty((len(folder_exps),nc,55,len(time_fv3)))
target_temp_hafs_ocean[:] = np.nan
target_salt_hafs_ocean = np.empty((len(folder_exps),nc,55,len(time_fv3)))
target_salt_hafs_ocean[:] = np.nan
target_depth_hafs_ocean = np.empty((len(folder_exps),nc,55))
target_depth_hafs_ocean[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    for c,cycle in enumerate(cycles):
        folderc = folder + '/' + cycle + '/' + storm_nums[c] + basins[c][-1] + '/'
        #%% Get list files
        if ocean[i] == 'hycom':
            files_hafs_ocean = sorted(glob.glob(os.path.join(folderc,'*3z*.nc')))
            hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
            lon_hafs_ocean = np.asarray(hafs_ocean['Longitude'][:])
            lat_hafs_ocean = np.asarray(hafs_ocean['Latitude'][:])
            depth_hafs_ocean = np.asarray(hafs_ocean['Z'][:])
        if ocean[i] == 'mom6':
            files_hafs_ocean = sorted(glob.glob(os.path.join(folderc,'*mom6*.nc')))
            hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
            lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
            lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
            depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])
    
        #%% Get storm track from trak atcf files
        if hafs_ab[i] == 'hfsa':
            file_track = folderc + storm_ids[c]+'.' + cycle + '.hfsa.trak.atcfunix'
        if hafs_ab[i] == 'hfsb':
            file_track = folderc + storm_ids[c] +'.' + cycle + '.hfsb.trak.atcfunix'
        print(file_track)
    
        okn = get_storm_track_and_int(file_track,storm_nums[c])[0].shape[0]
        lon_forec_track[i,c,0:okn], lat_forec_track[i,c,0:okn], lead_time[i,c,0:okn], int_track[i,c,0:okn], _ = get_storm_track_and_int(file_track,storm_nums[c])
    
        #%% Read time
        time_ocean = []
        timestamp_ocean = []
        for n,file in enumerate(files_hafs_ocean):
            print(file)
            hafs_ocean = xr.open_dataset(file)
            if ocean[i] == 'hycom':
                t = hafs_ocean['MT'][:]
            if ocean[i] == 'mom6':
                t = hafs_ocean['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_ocean.append(mdates.num2date(timestamp))
        timestamp_ocean.append(timestamp)
    
        time_ocean = np.asarray(time_ocean)
        timestamp_ocean = np.asarray(timestamp_ocean)
    
        #%% Retrieve glider transect
        ncfiles = files_hafs_ocean
        lon = lon_hafs_ocean
        lat = lat_hafs_ocean
        depth = depth_hafs_ocean
        if ocean[i] == 'hycom':
            time_name = 'MT'
            temp_name = 'temperature'
            salt_name = 'salinity'
        if ocean[i] == 'mom6':
            time_name = 'time'
            temp_name = 'temp'
            salt_name = 'so'
    
        if np.min(lon) < 0:
            lon_glid = longg
        else: 
            lon_glid = lon_glider
        lat_glid = lat_glider
    
        target_t, target_temp_hafs_oc, target_salt_hafs_oc = \
        get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)
    
        target_temp_hafs_ocean[i,c,0:target_temp_hafs_oc.shape[0],0:target_temp_hafs_oc.shape[1]] = target_temp_hafs_oc
        target_salt_hafs_ocean[i,c,0:target_temp_hafs_oc.shape[0],0:target_temp_hafs_oc.shape[1]] = target_salt_hafs_oc
        target_depth_hafs_ocean[i,c,0:len(depth_hafs_ocean)] = depth_hafs_ocean
    
        target_time[i,c,0:target_temp_hafs_oc.shape[1]] = mdates.date2num(target_t)
  
##################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
#okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for c,cyc in enumerate(cycles):
    for i in np.arange(len(exp_names)): 
        if c == 0:
            plt.plot(lon_forec_track[i,c,::2], lat_forec_track[i,c,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
        else:
            plt.plot(lon_forec_track[i,c,::2], lat_forec_track[i,c,::2],'o-',color=exp_colors[i],markeredgecolor='k',markersize=7)
        if c==0 and i==0:
            plt.plot(lon_best_track[c], lat_best_track[c],'o-',color='k',label='Best Track')
        else:
            plt.plot(lon_best_track[c], lat_best_track[c],'o-',color='k')
plt.plot(long[0,:], latg[0,:],'.-',color='orange',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_nums[0] + ' cycle '+ cycles[0],fontsize=18)
plt.axis('scaled')

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
#okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for c,cyc in enumerate(cycles):
    for i in np.arange(len(exp_names)):
        if c == 0:
            plt.plot(lon_forec_track[i,c,::2], lat_forec_track[i,c,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
        else:
            plt.plot(lon_forec_track[i,c,::2], lat_forec_track[i,c,::2],'o-',color=exp_colors[i],markeredgecolor='k',markersize=7)
        if c==0 and i==0:
            plt.plot(lon_best_track[c], lat_best_track[c],'o-',color='k',label='Best Track')
        else:
            plt.plot(lon_best_track[c], lat_best_track[c],'o-',color='k')
plt.plot(long[0,:], latg[0,:],'.-',color='orange',label='Glider Track')
plt.legend(loc='upper right',bbox_to_anchor=[0.4,1.0])
#plt.legend()
plt.title('Track Forecast ' + storm_nums[0] + ' cycle '+ cycles[0],fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(long)-4,np.nanmax(long)+4])
plt.ylim([np.nanmin(latg)-4,np.nanmax(latg)+4])

##################################################################
# Figures vertical profiles
y = int(teyes[0][0:4])
m = int(teyes[0][4:6])
d = int(teyes[0][6:8])
h = int(teyes[0][8:10])
time_eye_passage1 = mdates.date2num(datetime(y,m,d,h,0,0)) 
time_eye_passage1_plus48 = mdates.date2num(datetime(y,m,d,h,0,0)+ timedelta(hours=72)) 
before_eye1 = timeg[0,:] < time_eye_passage1 
after_eye1 = np.logical_and(timeg[0,:] >= time_eye_passage1,timeg[0,:] <= time_eye_passage1_plus48) 

y = int(teyes[1][0:4])
m = int(teyes[1][4:6])
d = int(teyes[1][6:8])
h = int(teyes[1][8:10])
time_eye_passage2_minus48 = mdates.date2num(datetime(y,m,d,h,0,0)-timedelta(hours=72)) 
time_eye_passage2 = mdates.date2num(datetime(y,m,d,h,0,0)) 
before_eye2 = np.logical_and(timeg[0,:] >= time_eye_passage2_minus48,timeg[0,:] < time_eye_passage2) 
after_eye2 = timeg[0,:] >= time_eye_passage2 

fig, ax = plt.subplots(figsize=(8,11))
plt.plot(tempg[0,before_eye1],-depthg[0,before_eye1],'.',color='r',label='Before Helene')
plt.plot(tempg[0,after_eye1],-depthg[0,after_eye1],'.',color='b',label='After Helene')
plt.plot(tempg[:,before_eye1],-depthg[:,before_eye1],'.',color='r',alpha=0.5,markersize=1)
plt.plot(tempg[:,after_eye1],-depthg[:,after_eye1],'.',color='b',alpha=0.3,markersize=1)
plt.plot(tempg[0,before_eye2],-depthg[0,before_eye2],'.',color='orange',label='Before Milton')
plt.plot(tempg[0,after_eye2],-depthg[0,after_eye2],'.',color='cyan',label='After Milton')
plt.plot(tempg[:,before_eye2],-depthg[:,before_eye2],'.',color='orange',alpha=0.5,markersize=1)
plt.plot(tempg[:,after_eye2],-depthg[:,after_eye2],'.',color='cyan',alpha=0.3,markersize=1)
plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Temperature ($^oC$)',fontsize=16)
plt.legend()
plt.title(label=dataset_id,fontsize=16)

# Profilesn before and after storm passage
for i,ff in enumerate(folder_exps):
    fig, ax = plt.subplots(figsize=(8,11))
    before_eye = target_time[i,0,:] < time_eye_passage1
    after_eye = target_time[i,0,:] > time_eye_passage1
    plt.plot(target_temp_hafs_ocean[i,0,:,before_eye][0,:].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,before_eye].shape[0],1))[0,:].T,'o-',color='red',label='Before Helene',markeredgecolor='k')
    plt.plot(target_temp_hafs_ocean[i,0,:,before_eye].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,before_eye].shape[0],1)).T,'o-',color='red',markeredgecolor='k',alpha=0.5)
    plt.plot(target_temp_hafs_ocean[i,0,:,after_eye][0,:].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,after_eye].shape[0],1))[0,:].T,'o-',color='blue',label='After Helene',markeredgecolor='k')
    plt.plot(target_temp_hafs_ocean[i,0,:,after_eye].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,after_eye].shape[0],1)).T,'o-',color='blue',markeredgecolor='k',alpha=0.5)

    before_eye = target_time[i,1,:] < time_eye_passage2
    after_eye = target_time[i,1,:] > time_eye_passage2
    plt.plot(target_temp_hafs_ocean[i,1,:,before_eye][0,:].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,before_eye].shape[0],1))[0,:].T,'o-',color='orange',label='Before Milton',markeredgecolor='k')
    plt.plot(target_temp_hafs_ocean[i,1,:,before_eye].T,-np.tile(target_depth_hafs_ocean[i,0,:],(target_temp_hafs_ocean[i,0,:,before_eye].shape[0],1)).T,'o-',color='orange',markeredgecolor='k',alpha=0.5)
    plt.plot(target_temp_hafs_ocean[i,1,:,after_eye][0,:].T,-np.tile(target_depth_hafs_ocean[i,1,:],(target_temp_hafs_ocean[i,0,:,after_eye].shape[0],1))[0,:].T,'o-',label='After Milton',color='cyan',markeredgecolor='k')
    plt.plot(target_temp_hafs_ocean[i,1,:,after_eye].T,-np.tile(target_depth_hafs_ocean[i,1,:],(target_temp_hafs_ocean[i,0,:,after_eye].shape[0],1)).T,'o-',color='cyan',markeredgecolor='k',alpha=0.5)

    plt.ylabel('Depth (m)',fontsize=16)
    plt.xlabel('Temperature ($^oC$)',fontsize=16)
    plt.title(exp_labels[i] + ' ' + str(tini)+ '-' + str(tend),fontsize=16)
    plt.legend()
    plt.ylim([-160,0])
    plt.xlim([16,31])


