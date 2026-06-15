#%% User input

# forecasting cycle to be used
cycle = '2026010300'
url_drifter = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/Data/Lagrangian_Drifters_Scripps/LDL_sea_level_press_Jan_2026.nc'

exp_labels = ['uncoupled','Coupled']
exp_colors = ['blue','orange']

lon_lim = [-180,-80]
lat_lim = [0,70]

folder_exps = ['/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_uncoupled/','/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/']

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
#folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
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
import grib2io

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#sys.path.append(folder_myutils)
#from my_models_utils get_var_from_model_following_trajectory

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
def get_var_from_model_following_trajectory_grb2_file(files_model,var_name,lon_model,lat_model,lon_obs,lat_obs,timestamp_obs):

    target_var = np.empty((len(files_model)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        model = grib2io.open(file,mode='r')
        t = model.select(shortName=var_name)[0].validDate
        timestamp = t.timestamp()
        target_time.append(t)

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
        var_sort = model.select(shortName=var_name)[0].data[oksort_lat,:][:,oksort_lon]
    
        target_var[x] = var_sort[oklat,oklon]

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
sea_level_pressure = np.asarray(gdata.sea_level_pressure)

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
slpDr = sea_level_pressure[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
slpD = slpDr[oklat][oklon]
wmo_idD = wmo_idDr[oklat][oklon]

codes = np.unique(wmo_idD)

'''
oklon = np.logical_and(lonD <= -157,lonD >= -160)
oklat = np.logical_and(latD[oklon] <= 50,latD[oklon] >= 49)
np.unique(wmo_idD[oklon][oklat])

oklon = np.logical_and(lonD <= -156,lonD >= -160)
oklat = np.logical_and(latD[oklon] <= 20,latD[oklon] >= 16)
np.unique(wmo_idD[oklon][oklat])
'''

########################################################################
# Get lon and lat from model
files_fv3 = sorted(glob.glob(os.path.join(folder_exps[0]+cycle+'/00E/','arafs*.grb2')))
model = grib2io.open(files_fv3[0],mode='r')

lon_fv3 = model.select(shortName='ELON')[0].data
lat_fv3 = model.select(shortName='NLAT')[0].data

########################################################################
target_time = []
target_slp = np.empty((len(folder_exps),43))
target_slp[:] = np.nan

for code in [7801738.]:
    print(code)
    okcode = wmo_idD == code
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
    platform_coded = wmo_idD[okcode]
    slpdd = slpD[okcode]

    slpd = np.empty((len(slpdd)))
    slpd[:] = np.nan
    for i,slp in enumerate(slpdd):
        if slp == 0:
            slpd[i] = np.nan
        else:
            slpd[i] = float(slp)

    # Loop the experiments
    for i,folder in enumerate(folder_exps):
        print(folder)

        #%% Get list files
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder+cycle+'/00E/','arafs*.grb2')))

        ########################################################################
        #%% Retrieve model variable following saildrone trajectory

        #with grib2io.open(file, mode='r') as f:
        #    unique_vars = set(msg.shortName for msg in f)

        files_model = files_hafs_fv3
        lon_model = lon_fv3-360
        lat_model = lat_fv3

        timestamp_obsd = timestampd
        lon_obsd = lond
        lat_obsd = latd

        oklo = np.isfinite(lon_obsd)
        lon_obs = lon_obsd[oklo]
        lat_obs = lat_obsd[oklo]
        timestamp_obs = timestamp_obsd[oklo]

        target_timeD, target_slp[i,0:len(files_model)] = get_var_from_model_following_trajectory_grb2_file(files_model,'PRMSL',lon_model,lat_model,lon_obs,lat_obs,timestamp_obs)

        target_time.append(target_timeD)

        #######################################################################
    
    # Figure SST
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,slpd,'o-',color='k',label='Drifter '+ str(code),markersize=5,markeredgecolor='k')
    for i in np.arange(len(exp_labels)):
        plt.plot(target_timeD,target_slp[i,0:len(target_timeD)]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    #plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.legend(loc='upper left')
    plt.title('Sea Level Pressure Cycle '+ cycle,fontsize=18)
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


