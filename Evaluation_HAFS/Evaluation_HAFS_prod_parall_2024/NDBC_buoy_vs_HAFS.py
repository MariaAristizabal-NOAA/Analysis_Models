#%% User input
# forecasting cycle to be used

# Milton
storm_num = '09'
basin = 'al'
storm_id = '09l'
cycle = '2024092500'
year = '2024'
hafs = ['hfsa','hfsa']

# Milton
'''
storm_num = '14'
basin = 'al'
storm_id = '14l'
cycle = '2024100800'
year = '2024'
hafs = ['hfsa','hfsa']
'''

#url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2024/NDBC_buoy_42098.nc'
#url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2024/NDBC_buoy_42097.nc'
url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2024/NDBC_buoy_42036.nc'

exp_names = ['HFSA_oper','HAFSv2p0p1a_2024rt']
exp_labels = ['HFSA','HFXA']
exp_colors = ['purple','cadetblue']

'''
# Lee
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'
url_NDBC = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/2023/NDBC_buoy_Lee_2023.nc'

exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_HYCOM','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
exp_labels = ['HFSA_oper_MOM6_epbl','HWRF','HFSAv1p1_HYCOM_kpp','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
exp_colors = ['darkviolet','pink','forestgreen','cyan','royalblue']
'''

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

lon_lim = [-75,-55.0]
lat_lim = [20.0,50.0]

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import pygrib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            get_var_from_model_following_trajectory,\
                            get_var_from_ww3_grb2_following_trajectory

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
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
#%% Read Saildrone data
# Sam

url = url_NDBC

gdata = xr.open_dataset(url)#,decode_times=False)

stationB = gdata.station.values
latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
time = np.asarray(gdata.time)
wd = np.asarray(gdata.wd) # wind direction (the direction the wind is coming from in degrees clockwise from true N)
wspd = np.asarray(gdata.wspd) # Average wind speed (m/s)
wvht = np.asarray(gdata.wvht) # Significant wave height (meters) is calculated as the average of the highest one-third of all of the wave heights during the 20-minute sampling period.
dpd = np.asarray(gdata.dpd) # Dominant wave period (seconds) is the period with the maximum wave energy.
apd = np.asarray(gdata.apd) # Average wave period (seconds) of all waves during the 20-minute period.
mwd = np.asarray(gdata.mwd) # Mean wave direction corresponding to energy of the dominant period (DOMPD).
bar = np.asarray(gdata.bar) # Air pressure (hPa).
wtmp = np.asarray(gdata.wtmp) # Sea surface temperature (Celsius).
tide = np.asarray(gdata.tide) # The water level in meters.

times = np.asarray(gdata.time)
timestamp = mdates.date2num(time)
#times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
timesB = times[oktimeg]
timestampB = timestamp[oktimeg]
latB = latitude[oktimeg]
lonB = longitude[oktimeg]
wdB = wd[oktimeg]
wspdB = wspd[oktimeg]
wvhtB = wvht[oktimeg]
dpdB = dpd[oktimeg]
apdB = apd[oktimeg]
mwdB = mwd[oktimeg]
barB = bar[oktimeg]
wtmpB = wtmp[oktimeg]
tideB = tide[oktimeg]

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

target_timeB_fv3 = []
target_timeB_ocean = []
target_timeB_ww3 = []
target_wdB = np.empty((len(folder_exps),43))
target_wdB[:] = np.nan
target_wspdB = np.empty((len(folder_exps),43))
target_wspdB[:] = np.nan
target_wvhtB = np.empty((len(folder_exps),43))
target_wvhtB[:] = np.nan
target_dpdB = np.empty((len(folder_exps),43))
target_dpdB[:] = np.nan
target_apdB = np.empty((len(folder_exps),43))
target_apdB[:] = np.nan
target_mwdB = np.empty((len(folder_exps),43))
target_mwdB[:] = np.nan
target_barB = np.empty((len(folder_exps),43))
target_barB[:] = np.nan
target_wtmpB = np.empty((len(folder_exps),43))
target_wtmpB[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get storm track from trak atcf files
    if hafs[i] == 'hfsa':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if hafs[i] == 'hfsb':
        file_track = folder + + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'

    # Read track file
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

    #%% Get list files
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hfs*.parent.atm.*.grb2')))
    files_hafs_mom6 = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
    files_hafs_ww3 = sorted(glob.glob(os.path.join(folder,'*ww3*.grb2')))

    #################################################################################
    #%% Retrieve HAFS_atm press
    
    '''
    #%% Reading FV3 grid
    grbindx = pygrib.index(files_hafs_fv3[0],'shortName','typeOfLevel','level')
    selected_grbs = grbindx.select(shortName='t',typeOfLevel='surface',level=0)
    lat_hafs_fv3, lon_hafs_fv3 = selected_grbs[0].latlons()

    grib2files = files_hafs_fv3
    time_name = 'time'
    lat_name = 'latitude'
    lon_name = 'longitude'
    #var_name = 'prmsl'
    if np.min(lon_hafs_fv3) < 0:
        lon_obss = lonB
    else:
        lon_obss = lonB + 360
    lat_obss = latB
    timestamp_obss = timestampB

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    target_timeBfv, target_barB[i,0:len(grib2files)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2files,'sp',typeoflevel='surface',level='0')

    target_timeB_fv3.append(target_timeBfv)
    '''
    #################################################################################
    #%% Retrieve wave parameters

    #%% Reading WW3 grid
    grbs = pygrib.open(files_hafs_ww3[0])
    grbs.seek(0)
    for grb in grbs[1:30]:
        print(grb.name + '=' +  grb.shortName +',' + grb.typeOfLevel + ',' + str(grb.level))

    grbindx_ww3 = pygrib.index(files_hafs_ww3[0],'shortName','typeOfLevel','level')
    selected_grbs_ww3 = grbindx_ww3.select(shortName='ws',typeOfLevel='surface',level=1)
    lat_hafs_ww3, lon_hafs_ww3 = selected_grbs_ww3[0].latlons()

    grib2file = files_hafs_ww3[0]
    if np.min(lon_hafs_ww3) < 0:
        lon_obss = lonB
    else:
        lon_obss = lonB + 360
    lat_obss = latB
    timestamp_obss = timestampB

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    target_timeBww3, target_wspdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'ws',typeoflevel='surface',level=1)

    target_timeBww3, target_wdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'wdir',typeoflevel='surface',level=1)

    target_timeBww3, target_wvhtB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'swh',typeoflevel='surface',level=1)

    target_timeBww3, target_dpdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'perpw',typeoflevel='surface',level=1)

    target_timeBww3, target_apdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'mwp',typeoflevel='surface',level=1)

    target_timeBww3, target_mwdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'dirpw',typeoflevel='surface',level=1)

    target_timeB_ww3.append(target_timeBww3)

    #################################################################################
    #%% Retrieve ocean output

    #%% Reading oceam
    files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
    hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
    lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
    lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
    depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])
    time_name = 'time'
    temp_name = 'temp'
    salt_name = 'so'
    lon_name = 'xh'
    lat_name = 'yh'

    #%% Read time
    time_ocean = []
    timestamp_ocean = []
    for n,file in enumerate(files_hafs_ocean):
        print(file)
        hafs_ocean = xr.open_dataset(file)
        t = hafs_ocean[time_name][:]
    timestamp = mdates.date2num(t)[0]
    time_ocean.append(mdates.num2date(timestamp))
    timestamp_ocean.append(timestamp)
    time_ocean = np.asarray(time_ocean)
    timestamp_ocean = np.asarray(timestamp_ocean)

    ncfiles = files_hafs_ocean
    lon = lon_hafs_ocean
    lat = lat_hafs_ocean
    depth_level = 0
    timestamp_obss = timestampB
    kwargs = dict(depth_level = 0)

    oklo = np.isfinite(lon_obss)
    lon_obs = lonB[oklo]
    lat_obs = latB[oklo]
    timestamp_obs = timestamp_obss[oklo]

    target_timeBocean, target_wtmpB[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,ncfiles,temp_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name,depth_level=0)
        
    target_timeB_ocean.append(target_timeBocean)


#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonB, latB,'*',color='blue',label='NDBC Buoy '+np.unique(stationB)[0],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonB, latB,'*',color='blue',label='NDBC Buoy '+np.unique(stationB)[0],markersize=7)
plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonB)-3,np.nanmax(lonB)+3])
plt.ylim([np.nanmin(latB)-3,np.nanmax(latB)+3])


#################################################################################

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,wtmpB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ocean[i],target_wtmpB[i,0:len(target_timeB_ocean[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Water Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

'''
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,barB,'.-',color='blue',label='NDBC Buoy')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_fv3,target_barB[i,:]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.legend()
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
#plt.ylim([950,1050])
plt.ylim([1000,1020])
'''

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,wspdB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_wspdB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Wind Speed Cycle '+ cycle,fontsize=18)
plt.ylabel('($m/s$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,wvhtB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_wvhtB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Significant Wave Height Cycle '+ cycle,fontsize=18)
plt.ylabel('(meters)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,wdB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_wdB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Wind Direction Cycle '+ cycle,fontsize=18)
plt.ylabel('(degrees)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,dpdB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_dpdB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Dominant Wave Period Cycle '+ cycle,fontsize=18)
plt.ylabel('(Seconds)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,apdB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_apdB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Average Wave Period Cycle '+ cycle,fontsize=18)
plt.ylabel('(Seconds)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timesB,mwdB,'.-',color='blue',label='NDBC Buoy '+np.unique(stationB)[0])
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ww3[i],target_mwdB[i,0:len(target_timeB_ww3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
plt.legend(loc='upper right')
plt.title('Dominant Wave Direction Cycle '+ cycle,fontsize=18)
plt.ylabel('(Degrees)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

