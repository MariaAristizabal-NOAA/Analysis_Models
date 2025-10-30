#%% User input
# forecasting cycle to be used

# Erin
storm_num = '98'
basin = 'al'
storm_id = '98l'
cycle = '2025101818'
year = '2025'

'''
storm_num = '13'
basin = 'al'
storm_id = '13l'
cycle = '2025102112'
year = '2025'
'''

exp_names = ['HFSA_oper']
exp_labels = ['HFSA']
exp_colors = ['purple']
hafs = ['hfsa']

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

scratch_folder = '/scratch3/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch4/HFIP/hwrfv3/noscrub/input/abdeck/'
folder_myutils= '/home/Maria.Aristizabal/Maria_Utils/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

url_buoy = '/scratch3/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Buoy_Hurr_Melissa_2025/buoy_1.nc'

lon_lim = [-80,-50.0]
lat_lim = [10.0,40.0]

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
#%% Read  data
# Sam

url = url_buoy

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.lat)
longitude = np.asarray(gdata.lon)
time = np.asarray(gdata.time)
sst = np.asarray(gdata.sst)-273.15 
timestampp = mdates.date2num(time)
#times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
timeB = time[oktimeg]
timestampB = timestampp[oktimeg]
latB = latitude[oktimeg]
lonB = longitude[oktimeg]
sstB = sst[oktimeg]

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
target_sstB = np.empty((len(folder_exps),43))
target_sstB[:] = np.nan

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
    selected_grbs_ww3 = grbindx_ww3.select(shortName='ws',typeOfLevel='surface',level='1')
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

    target_timeBww3, target_wspdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'ws',typeoflevel='surface',level='1')

    target_timeBww3, target_wdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'wdir',typeoflevel='surface',level='1')

    target_timeBww3, target_wvhtB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'swh',typeoflevel='surface',level='1')

    target_timeBww3, target_dpdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'perpw',typeoflevel='surface',level='1')

    target_timeBww3, target_apdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'mwp',typeoflevel='surface',level='1')

    target_timeBww3, target_mwdB[i,:] = get_var_from_ww3_grb2_following_trajectory(lon_obs,lat_obs,timestamp_obs,grib2file,'dirpw',typeoflevel='surface',level='1')

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

    target_timeBocean, target_sstB[i,0:len(ncfiles)] = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,ncfiles,temp_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name,depth_level=0)
        
    target_timeB_ocean.append(target_timeBocean)


#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

okt = np.where(t_best_track == cycle)[0][0]

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track[okt:], lat_best_track[okt:],'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonB, latB,'*',color='blue',label='Buoy ',markersize=7)
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
plt.plot(lon_best_track[okt:], lat_best_track[okt:],'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonB, latB,'*',color='blue',label='Buoy',markersize=7)
plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
#plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonB)-8,np.nanmax(lonB)+8])
plt.ylim([np.nanmin(latB)-8,np.nanmax(latB)+8])

#################################################################################

fig,ax = plt.subplots(figsize=(10, 4))
#plt.plot(timestampB,sstB,'.-',color='blue',label='Buoy')
plt.plot(timestampp,sst,'.-',color='blue',label='Buoy')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeB_ocean[i],target_sstB[i,0:len(target_timeB_ocean[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
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


