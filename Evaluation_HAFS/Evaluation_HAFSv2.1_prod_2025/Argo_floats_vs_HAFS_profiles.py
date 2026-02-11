#!/usr/bin/env python3

"""This script is to compare ARAFS with Argo floats"""

import os
import glob
import sys
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

# Parse the yaml config file
print('Parse the config file: plot_atmos4.yml:')
with open('plot_atmos4.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
storm_id = conf['Stormid'] 
storm_num = conf['Stormid'][0:2] 
basin = conf['Basin'] 
cycle = conf['ymdh']
experiments = conf['Experiments']
exp_markers = conf['Exp_markers']
exp_colors = conf['Exp_colors']
data_prod = conf['Data_prod']
fhhhs = conf['fhhhs']
commodels = conf['COMmodels']
comobs = conf['COMobs']
fileobs = conf['Fileobs']
folder_myutils = conf['folder_myutils']
abdeck_folder = conf['abdeck_folder']

xlim = conf['xlim']
ylim = conf['ylim']

print(conf)

cartopy.config['data_dir'] = conf['cartopyDataDir']

sys.path.append(folder_myutils)
from my_models_utils import get_best_track_and_int, get_storm_track_and_int

################################################################################
def get_var_from_model_for_Argo_float_profile(files_model,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs):

    target_var = np.empty((len(depth_levels)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        model = xr.open_dataset(file)
        if file.split('.')[-1] == 'nc':
            t = model[time_name][:]
            timestamp = mdates.date2num(t)[0]
            target_time.append(mdates.num2date(timestamp))
            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])

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
        if file.split('.')[-1] == 'nc':
            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name][0,depth_levels,oksort_lat,oksort_lon])
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name][0,oksort_lat,oksort_lon])
        if len(depth_levels) > 1:
            target_var[:,x] = var[:,oklat,oklon]
        else:
            target_var[x] = var[oklat,oklon]

    return target_time, target_var

#####################################################################
# Read  temp from observations products
ncfile = comobs[0] + '/' + fileobs[0]
OBS = xr.open_dataset(ncfile)

time_obs = np.asarray(OBS['time'])
lat_obs = np.asarray(OBS['latitude'])
lon_obs = np.asarray(OBS['longitude'])
temp_obs = np.asarray(OBS['temp']) 
temp_qc_obs = np.asarray(OBS['temp_qc']) 
salt_obs = np.asarray(OBS['psal']) 
salt_qc_obs = np.asarray(OBS['psal_qc']) 
depth_obs = np.asarray(OBS['pres']) 
depth_qc_obs = np.asarray(OBS['pres_qc']) 

depth_qc_obs = np.asarray(OBS['pres_qc']) 

#timestamp_obs = time_obs.astype('int')/1e9
timestamp_obs = mdates.date2num(time_obs)

total_nobs = len(np.asarray(OBS['platform_number']))
total_plat = len(np.unique(np.asarray(OBS['platform_number'])))
plat_ids = np.unique(np.asarray(OBS['platform_number']))

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')
plt.plot(lon_obs,lat_obs,'.',label='ARGO = '+ str(total_nobs))
plt.legend()
plt.title(str(time_obs[0]).split('T')[0])
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

####################################################################
# Read best track
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

####################################################################
# Read forecast track
lon_forec_track = np.empty((len(experiments),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(experiments),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(experiments),43))
lead_time[:] = np.nan
int_track = np.empty((len(experiments),43))
int_track[:] = np.nan
for ex,exp in enumerate(experiments):
    #print(exp)
    #%% Get storm track from trak atcf files
    if experiments[ex] == 'HFSA' or experiments[ex] == 'HFXA':
        file_track = commodels[ex] + '/' + storm_id + '.' + cycle + '.hfsa.trak.atcfunix'
    if experiments[ex] == 'HFSB' or experiments[ex] == 'HFXB':
        file_track = commodels[ex] + '/' + storm_id + '.' + cycle + '.hfsb.trak.atcfunix'

    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[ex,0:okn], lat_forec_track[ex,0:okn], lead_time[ex,0:okn], int_track[ex,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

okt = np.where(t_best_track == cycle)[0][0] 

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')
plt.plot(lon_best_track[okt:],lat_best_track[okt:],'o-',color='k',label='Best Track')
for i in np.arange(len(experiments)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=experiments[i],markersize=7)
plt.plot(lon_obs,lat_obs,'o',color='green',label='ARGOS')
plt.legend()
plt.title(str(time_obs[0]).split('T')[0])
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
plt.xlim(xlim)
plt.ylim(ylim)

########################################################################
for plat_id in plat_ids:
    print(plat_id)
    okid = np.asarray(OBS['platform_number']) == plat_id
    timestamp_obs_unique = np.unique(timestamp_obs[okid])

    fig = plt.figure(figsize=(6,8))
    
    for tt in timestamp_obs_unique:
        oktime = timestamp_obs[okid] == tt 
        timepp = time_obs[okid][oktime]
        timestamppp = timestamp_obs[okid][oktime]
        latpp = lat_obs[okid][oktime]
        lonpp = lon_obs[okid][oktime]
        depthpp = depth_obs[okid][oktime]
        depthpp_qc = depth_qc_obs[okid][oktime]
        temppp = temp_obs[okid][oktime]
        temppp_qc = temp_qc_obs[okid][oktime]
        saltpp = salt_obs[okid][oktime]
        saltpp_qc = salt_qc_obs[okid][oktime]
        
        # Use only observations with qc flag = '1', good data
        good = temppp_qc == '1'
        timep = timepp[good] 
        timestampp = timestamppp[good]
        latp = latpp[good] 
        lonp = lonpp[good] 
        depthp = depthpp[good]
        tempp = temppp[good] 
        saltp = saltpp[good] 

        if len(timep) == 0:
            print('Bad data')
            break

    for ex,exp in enumerate(experiments):
        print(exp)
        files_model = np.sort(glob.glob(conf['COMmodels'][ex]+'/*mom6*.nc'))

        times_model = []
        timestamps_model = []

        for x,file in enumerate(files_model):
            model = xr.open_dataset(file)
            t = np.asarray(model['time'][:])
            timestamp = mdates.date2num(t)[0]
            times_model.append(mdates.num2date(timestamp))
            timestamps_model.append(timestamp)

        ok_file_model = int(np.round(np.interp(tt,timestamps_model,np.arange(len(timestamps_model)))))

        file_model = files_model[ok_file_model]
        nc_model = xr.open_dataset(file_model)
        time_model = np.asarray(nc_model['time'])
        lon_model = np.asarray(nc_model['xh'])
        lat_model = np.asarray(nc_model['yh'])

        latp_target = np.unique(latp)[0]
        lonp_target = np.unique(lonp)[0]
        oklat = int(np.round(np.interp(latp_target,lat_model,np.arange(len(lat_model)))))
        oklon = int(np.round(np.interp(lonp_target,lon_model,np.arange(len(lon_model)))))
        temp_model =  np.asarray(nc_model['temp'][0,:,oklat,oklon])
        depth_model = np.asarray(nc_model['z_l'])

        #################################################################

        plt.plot(temp_model,-depth_model,'o-',label=exp,color=exp_colors[ex],alpha=0.5)

    plt.plot(tempp,-depthp,'.-',color='green',label='ARGO '+plat_id)
    plt.legend()
    plt.xlabel('Temperature $^oC$')
    plt.ylabel('Depth (m)')
    plt.title('Lon = '+str(np.round(lonp[0],2))+' Lat = '+str(np.round(latp[0],2))+' time =  '+str(timep[0])[0:16])       #plt.ylim([np.min(-depthp),0])
    plt.ylim([-250,0])
    plt.xlim([15,30])
     
##################################################################




