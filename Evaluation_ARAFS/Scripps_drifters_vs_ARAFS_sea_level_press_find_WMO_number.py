#%% User input

# forecasting cycle to be used
cycle = '2025120900'

# 2023-2026
#WMO_AR_recon = [7810735,7810736,7810737,7810738,7810739,7810740,7810741,7810742,7810743,7810744,7810745,7810746,7810747,7810748,7810749,7810750,7810751,7810752,7810753,7810754,7810755,7810756,7810757,7810758,7810759,7810760,7810761,7810762,7810763,7810764,7810847,7810848,7810849,7810850,7810851,7810852,7810853,7810854,7810855,7810856,7810857,7810858,7810859,7810860,7810861,7810862,7810863,7810864,7810865,7810866,7810837,7810838,7810839,7810840,7810841,7810842,7810843,7810844,7810845,7810846,7810555,2802108,1801806,2802107,6801916,6801917,6801914,5802101,2802104,6801915,3801707,7801736,2802106,7801735,5802104,2802110,7801738,1801808,3801709,7801737,2802109,1801809,7801734,4804122,2802111,3801708,2802105,7801733,5802102,5802103,1801807,6801936,4804132,6801937,1801818,7801759,5802124,2802129,2802128,7801758,1801821,7801761,6801934,4804134,1801820,4804135,3801716,6801938,2802126,5802122,4804133,6801939,2802127,1801819,6801935,5802125,5802123,1801822,3801717,7801760,1801823,5802111,3801715,7801746,2802121,7801750,4804142,1801827,6801948,7801768,6801949,6801947,3801723,1801828,4804140,2802133,3801722,3801725,3801726,3801724,4804141,7801770,6801946,7801769]

# 2025-2026
WMO_AR_recon = [7810735,7810736,7810737,7810738,7810739,7810740,7810741,7810742,7810743,7810744,7810745,7810746,7810747,7810748,7810749,7810750,7810751,7810752,7810753,7810754,7810755,7810756,7810757,7810758,7810759,7810760,7810761,7810762,7810763,7810764,7810847,7810848,7810849,7810850,7810851,7810852,7810853,7810854,7810855,7810856,7810857,7810858,7810859,7810860,7810861,7810862,7810863,7810864,7810865,7810866,7810837,7810838,7810839,7810840,7810841,7810842,7810843,7810844,7810845,7810846]

url_drifter = '/scratch3/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_sea_level_press_Dec_1_2025_to_Mar_31_2026.nc'
#url_drifter = '/work/noaa/hwrf/noscrub/maristiz//Data/Scripts_lagrang_drifters/LDL_sea_level_press_Dec_1_2025_to_Mar_31_2026.nc'

exp_labels = ['uncoupled','Coupled']
exp_colors = ['blue','orange']

lon_lim = [-180,-80]
lat_lim = [0,70]

#folder_exps = ['/work/noaa/hurricane/malasala/Maria_Murali/ARAFSv1_coupled/']
folder_exps = ['/scratch4/HFIP/hafs-west/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_uncoupled/','/scratch4/HFIP/hafs-west/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/']

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import sys
import os
import glob
import grib2io
import random

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
def generate_unique_rgb_colors(n):
    if n > 256**3:
        raise ValueError("Requested more colors than possible unique RGB values.")

    colors = np.ones((n,3))
    for i in np.arange(n):
        # Generate a random RGB tuple
        color = [random.randint(0, 255)/255,
                 random.randint(0, 255)/255,
                 random.randint(0, 255)/255]
        colors[i,:] = color

    return colors

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

###############################################################################
#%% Read drifter data
url = url_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
wmo_id = np.asarray(gdata.wmo_ID)
wmo_id_unique = np.unique(wmo_id)
sea_level_pressure = np.asarray(gdata.sea_level_pressure)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(times)
#oktimeg = np.logical_and(mdates.date2num(times) >= mdates.date2num(tini),\
#                         mdates.date2num(times) <= mdates.date2num(tend))

# figure of all lagragian drifter in the publicably available file
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.plot(longitude, latitude,'.',color='k',markersize=1)
#plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
plt.title('Time Window: '+ str(times[0])[0:10] + ' - ' + str(times[-1])[0:10]+'\n Total number drifters = '+str(len(wmo_id_unique)),fontsize=18)
plt.axis('scaled')
#plt.xlim(lon_lim)
#plt.ylim(lat_lim)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latitude >= lat_lim[0], latitude <= lat_lim[1])
lonDD = longitude[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = times[oklat][oklon]
timestampD = timestamps[oklat][oklon]
latD = latitude[oklat][oklon]
lonD = longitude[oklat][oklon]
slpD = sea_level_pressure[oklat][oklon]
wmo_idD = wmo_id[oklat][oklon]

codes = np.unique(wmo_idD)

# figure of all lagragian drifter in the publicably available file
# within lat and lon limits
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.plot(lonD, latD,'.',color='k',markersize=1)
plt.title('Time Window: '+ str(times[0])[0:10] + ' - ' + str(times[-1])[0:10]+'\n Total number drifters = '+str(len(codes)),fontsize=18)
plt.axis('scaled')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False

# figure of all lagragian drifter in the publicably available file
# within lat and lon limits
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
plt.title('Time Window: '+ str(times[0])[0:10] + ' - ' + str(times[-1])[0:10]+'\n Total number drifters = '+str(len(codes)),fontsize=18)
plt.axis('scaled')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
plt.xlim(lon_lim)
plt.ylim(lat_lim)

colors = generate_unique_rgb_colors(len(codes))

for c,code in enumerate(codes):
    print(code)
    ok_id = wmo_idD == code
    ax.plot(lonD[ok_id][1:-1],latD[ok_id][1:-1],'.',markersize=0.02,color=colors[c])
    ax.plot(lonD[ok_id][0],latD[ok_id][0],'o',color=colors[c],markersize=4,markeredgecolor='k')
    ax.plot(lonD[ok_id][-1],latD[ok_id][-1],'o',color=colors[c],markersize=4,markerfacecolor='none')

# Fields within time window
oktimeg = np.logical_and(mdates.date2num(timeD) >= mdates.date2num(tini),\
                         mdates.date2num(timeD) <= mdates.date2num(tend))

timeDr = timeD[oktimeg]
timestampDr = timestampD[oktimeg]
latDr = latD[oktimeg]
lonDr = lonD[oktimeg]
wmo_idDr = wmo_idD[oktimeg]
slpDr = slpD[oktimeg]
codes = np.unique(wmo_idDr)

########################################################################
# Find drifters that are part of the AR recon campaing
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
#plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
plt.title('Time Window: '+ str(times[0])[0:10] + ' - ' + str(times[-1])[0:10],fontsize=18)
plt.axis('scaled')
plt.xlim(lon_lim)
plt.ylim(lat_lim)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False

wmo_ar_recon = []
wmo_ar_recon_wc = []
for c,wmo_num in enumerate(WMO_AR_recon):
    ok_id = wmo_id == wmo_num
    if len(wmo_id[ok_id]) != 0:
        #print(np.unique(wmo_id[ok_id])[0])
        wmo_ar_recon.append(np.unique(wmo_id[ok_id])[0])

        # Fields for this specific drifter
        latdd = latitude[ok_id]
        londd = longitude[ok_id]
        wmo_iddd = wmo_id[ok_id]

        # Find drifter within lon limits
        oklon = np.logical_and(londd >= lon_lim[0],londd < lon_lim[1])

        # Fields within the lon and lat limits
        if len(londd[oklon]) > 0:
            latld = latdd[oklon]
            lonld = londd[oklon]
            wmo_ar_recon_wc.append(np.unique(wmo_iddd[oklon])[0])
            print(np.unique(wmo_iddd[oklon])[0])

            ax.plot(lonld[1:-1],latld[1:-1],'.',markersize=0.02,color=colors[c])
            ax.plot(lonld[0],latld[0],'o',color=colors[c],markersize=4,markeredgecolor='k')
            ax.plot(lonld[-1],latld[-1],'o',color=colors[c],markersize=4,markerfacecolor='none')

#for wmo_num in wmo_ar_recon_wc:
#for wmo_num in [7801736.,2802110.,2802121.,3801723.]:
for wmo_num in [7810735.]:
    ok_id = wmo_id == wmo_num
    print(np.unique(wmo_id[ok_id])[0])

    # Fields for this specific drifter
    timedd = times[ok_id]
    timestampdd = timestamps[ok_id]
    latdd = latitude[ok_id]
    londd = longitude[ok_id]
    wmo_iddd = wmo_id[ok_id]
    code_id = np.unique(wmo_iddd)[0]
    slpdd = sea_level_pressure[ok_id]

    print(timedd[0])

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.plot(londd, latdd,'.',color='k',label='Lagrangian Drifter: '+str(code_id))
    plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
    plt.axis('scaled')
    plt.xlim(lon_lim)
    plt.ylim(lat_lim)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

    plt.figure(figsize=(12,5))
    plt.plot(timedd,slpdd,'.-',label='Lagrangian Drifter: '+str(code_id))
    plt.title('Sea Level Pressure')
    plt.legend()

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

#for code in [7801736.,2802110.,2802121.,3801723.]:
for code in [7810735.]:
    print(code)
    # Fields within time window and lat and lon limits
    okcode = wmo_idD == code
    timed = timeD[okcode]
    timestampd = timestampD[okcode]
    latd = latD[okcode]
    lond = lonD[okcode]
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
    
    # Figure SLP
    fig,ax = plt.subplots(figsize=(10, 4))
    plt.plot(timed,slpd,'o-',color='k',label='Drifter '+ str(code),markersize=5,markeredgecolor='k')
    for i in np.arange(len(folder_exps)):
        plt.plot(target_timeD,target_slp[i,0:len(target_timeD)]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    #plt.legend(loc='upper right',bbox_to_anchor=[1.7,1.0])
    plt.legend(loc='upper right')
    plt.title('Sea Level Pressure Cycle '+ cycle,fontsize=18)
    plt.ylabel('($^oC$)',fontsize=14)
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    
    fig = plt.figure(figsize=(6,7))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.axis('scaled')
    ax.plot(lonD, latD,'.',color='k',label='Lagrangian Drifters: '+str(len(lonD))+' Data Points')
    ax.plot(lon_obsd, lat_obsd,'.',color='red',label='Lagrangian Drifter '+str(code)+': '+str(len(lon_obsd))+' Data Points')
    plt.legend(loc='lower right',bbox_to_anchor=[1.0,-0.2])
    plt.title('Time Window: '+ str(tini)[0:13] + ' - ' + str(tend)[0:13],fontsize=18)
    plt.axis('scaled')
    plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
    plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])
    
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False


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

