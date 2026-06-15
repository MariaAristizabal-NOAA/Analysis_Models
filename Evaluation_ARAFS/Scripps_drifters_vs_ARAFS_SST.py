#%% User input

# forecasting cycle to be used
cycle = '2026010300'
url_drifter = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/Data/Lagrangian_Drifters_Scripps/LDL_SST_Jan_2026.nc'

exp_labels = ['Coupled']
exp_colors = ['orange']

lon_lim = [-180,-80]
lat_lim = [0,70]

folder_exps = ['/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/']

# folder utils for Hycom 
#folder_myutils= '/home/Maria.Aristizabal/Utils/'

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

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#sys.path.append(folder_myutils)
#from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
#                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
#                            get_var_from_model_following_trajectory

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
def get_var_from_model_following_trajectory_nc_file(files_model,var_name,time_name,lon_model,lat_model,lon_obs,lat_obs,timestamp_obs,**kwargs):
    # depth_level is an optional argument
    # kwargs = dict(depth_level = 0)

    depth_level = kwargs.get('depth_level', None)

    target_var = np.empty((len(files_model)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        #model = xr.open_dataset(file,engine="pynio")
        model = xr.open_dataset(file)
        t = model[time_name][:]
        timestamp = mdates.date2num(t)[0]
        target_time.append(mdates.num2date(timestamp))

        # Interpolating lat_obs and lon_obs into model grid
        sublon = np.interp(timestamp,timestamp_obs,lon_obs)
        sublat = np.interp(timestamp,timestamp_obs,lat_obs)
        if lon_model.ndim == 1:
            oksort_lon = np.argsort(lon_model)
            oksort_lat = np.argsort(lat_model)
            lonn = lon_model[oksort_lon]
            latt = lat_model[oksort_lat]
        if lon_model.ndim == 2:
            oksort_lon = np.argsort(lon_model[0,:])
            oksort_lat = np.argsort(lat_model[:,0])
            lonn = lon_model[0,oksort_lon]
            latt = lat_model[oksort_lat,0]
        oklon = int(np.round(np.interp(sublon,lonn,np.arange(len(lonn)))))
        oklat = int(np.round(np.interp(sublat,latt,np.arange(len(latt)))))
        if model[var_name].ndim == 4:
            var = np.asarray(model[var_name][0,depth_level,oksort_lat,oksort_lon])
        if model[var_name].ndim == 3:
            var = np.asarray(model[var_name][0,oksort_lat,oksort_lon])
        target_var[x] = var[oklat,oklon]

    return target_time, target_var

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=120)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

################################################################################
#%% Time fv3
#time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,121,3)]

###############################################################################
#%% Read drifter data
url = url_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
wmo_id = np.asarray(gdata.wmo_ID)
sea_surface_temp = np.asarray(gdata.sea_surface_temperature)

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
wmo_idDr = wmo_id[oktimeg]
sstDr = sea_surface_temp[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
sstD = sstDr[oklat][oklon]
wmo_idD = wmo_idDr[oklat][oklon]

codes = np.unique(wmo_idD)

'''
oklon = np.logical_and(lonD <= -157,lonD >= -160)
oklat = np.logical_and(latD[oklon] <= 50,latD[oklon] >= 49)
np.unique(wmo_idD[oklon][oklat])

oklon = np.logical_and(lonD <= -160,lonD >= -170)
oklat = np.logical_and(latD[oklon] <= 20,latD[oklon] >= 16)
np.unique(wmo_idD[oklon][oklat])
array([1801727., 2802018., 7810256.])

oklon = np.logical_and(lonD <= -123,lonD >= -127)
oklat = np.logical_and(latD[oklon] <= 50,latD[oklon] >= 46)
np.unique(wmo_idD[oklon][oklat])
[7801738.]
'''

########################################################################
target_timeD_ocean = []
target_timeD_atm = []
target_temp = np.empty((len(folder_exps),43))
target_temp[:] = np.nan

okc = codes == [1801723.0]
#for code in codes[okc]:
#for code in [7810738., 7810749., 7810750.]:
for code in [7801738.]:
    print(code)
    okcode = wmo_idD == code
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
    platform_coded = wmo_idD[okcode]
    sstdd = sstD[okcode]

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

        #%% Get list files
        files_hafs_mom6 = sorted(glob.glob(os.path.join(folder+cycle+'/00E/','*mom6*.nc')))

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

        ########################################################################
        #%% Retrieve model temp. following saildrone trajectory
        files_model = files_hafs_mom6
        time_name = 'time'
        lat_model = lat_hafs_mom6
        lon_model = lon_hafs_mom6
        depth_level = 0
        timestamp_obsd = timestampd
        kwargs = dict(depth_level = 0)
        lon_obsd = lond
        lat_obsd = latd

        oklo = np.isfinite(lon_obsd)
        lon_obs = lon_obsd[oklo]
        lat_obs = lat_obsd[oklo]
        timestamp_obs = timestamp_obsd[oklo]

        target_timeDmom6, target_temp[i,0:len(files_model)] = get_var_from_model_following_trajectory_nc_file(files_model,'temp',time_name,lon_model,lat_model,lon_obs,lat_obs,timestamp_obs,**kwargs)

        target_timeD_ocean.append(target_timeDmom6)

        #######################################################################
    
        # Figure SST
        fig,ax = plt.subplots(figsize=(10, 4))
        plt.plot(timed,sstd,'o-',color='k',label='Drifter '+ str(code),markersize=5,markeredgecolor='k')
        for i in np.arange(len(exp_labels)):
            plt.plot(target_timeD_ocean[i],target_temp[i,0:len(target_timeD_ocean[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
        #plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
        plt.legend(loc='upper left')
        plt.title('Sea Surface Temperature Cycle '+ cycle,fontsize=18)
        plt.ylabel('($^oC$)',fontsize=14)
        date_form = DateFormatter("%m-%d")
        ax.xaxis.set_major_formatter(date_form)
        
        fig = plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax.axis('scaled')
        ax.plot(lonD, latD,'.',color='k',label='Lagrangian Drifters: '+str(len(lonD))+' Data Points')
        ax.plot(lon_obsd, lat_obsd,'.',color='red',label='Lagrangian Drifters '+str(code)+': '+str(len(lon_obsd))+' Data Points')
        plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
        plt.title('Time Window: '+ str(tini)[0:13] + ' - ' + str(tend)[0:13],fontsize=18)
        plt.axis('scaled')
        plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
        plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])
        
        ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    
    
    ########################################################################

fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')
ax.plot(lonD, latD,'.',color='k',label='Lagrangian Drifters: '+str(len(lonD))+' Data Points')
plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
plt.title('Time Window: '+ str(tini)[0:13] + ' - ' + str(tend)[0:13],fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')


