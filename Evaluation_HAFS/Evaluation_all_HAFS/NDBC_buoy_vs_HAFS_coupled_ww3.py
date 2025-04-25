#%% User input

# forecasting cycle to be used
cycle = '2020082506'

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

list_exp = ['a2o2a_a2w2a','a2o_a2w2a']

folder_hafs_phase3_wave = scratch_folder + 'HAFSv0p2a_phase3_wave/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = scratch_folder + 'bdeck/bal132020.dat'
GFS_track_file = scratch_folder + 'adeck/aal132020.dat'

# NDBC buoy location
folder_ndbc = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NDBC_buoys/'
ndbc_buoys = ['42360h2020.nc','42395h2020.nc']

#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = scratch_folder + 'figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
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
                            Haversine, get_glider_transect_from_HAFS_HYCOM,\
                            figure_transect_temp, glider_data_vector_to_array,\
                            grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
def angle_0_360(u,v):
    ang = np.rad2deg(np.arctan(v/u))
    if u<0:
        ang = 180 + ang 
    else:
        if v>0:
            ang = ang
        else: 
            ang = 360 + ang
    return ang

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

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
lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)

#################################################################################
#%% Read ww3 
i = 0
files_hafs_phase3_wave = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'/'+cycle+'/00L/'+'hafs_202108_'+list_exp[i]+'/com/'+cycle+'/00L/','*ww3_ounf.nc')))
WW3 = xr.open_dataset(files_hafs_phase3_wave[0],decode_times=True)

# read WW3 time
t = np.asarray(WW3.variables['time'][:])
timestamp_ww3 = mdates.date2num(t)
time_ww3 = np.asarray(mdates.num2date(timestamp_ww3))

lat_ww3 = np.asarray(WW3.variables['latitude'])
lon_ww3 = np.asarray(WW3.variables['longitude'])

uwind_ww3_a2o2a_a2w2a = np.asarray(WW3.variables['uwnd'][:])
vwind_ww3_a2o2a_a2w2a = np.asarray(WW3.variables['vwnd'][:])
hs_ww3_a2o2a_a2w2a = np.asarray(WW3.variables['hs'][:])
fp_ww3_a2o2a_a2w2a = np.asarray(WW3.variables['fp'][:])

#################################################################################
#%% Read ww3
i = 1
files_hafs_phase3_wave = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'/'+cycle+'/00L/'+'hafs_202108_'+list_exp[i]+'/com/'+cycle+'/00L/','*ww3_ounf.nc')))
WW3 = xr.open_dataset(files_hafs_phase3_wave[0],decode_times=True)

# read WW3 time
t = np.asarray(WW3.variables['time'][:])
timestamp_ww3 = mdates.date2num(t)
time_ww3 = np.asarray(mdates.num2date(timestamp_ww3))

lat_ww3 = np.asarray(WW3.variables['latitude'])
lon_ww3 = np.asarray(WW3.variables['longitude'])

uwind_ww3_a2o_a2w2a = np.asarray(WW3.variables['uwnd'][:])
vwind_ww3_a2o_a2w2a = np.asarray(WW3.variables['vwnd'][:])
hs_ww3_a2o_a2w2a = np.asarray(WW3.variables['hs'][:])
fp_ww3_a2o_a2w2a = np.asarray(WW3.variables['fp'][:])


#################################################################################
#%% Read NDBC buoy
i = 1
file_ndbc = folder_ndbc + ndbc_buoys[i]
ndbc = xr.open_dataset(file_ndbc,decode_times=True)

tt = np.asarray(ndbc['time'][:])
timestamp_ndbc = mdates.date2num(tt)
tt_ndbc = np.asarray(mdates.num2date(timestamp_ndbc))

okt = np.where(np.logical_and(tt_ndbc >= time_ww3[0],tt_ndbc <= time_ww3[-1]))[0] 
time_ndbc = tt_ndbc[okt]
lat_ndbc = np.asarray(ndbc['latitude'][:])
lon_ndbc = np.asarray(ndbc['longitude'][:])
hs_ndbc = np.asarray(ndbc['wave_height'][okt,0,0])
fp_ndbc = (np.asarray(ndbc['dominant_wpd'][okt,0,0])/np.timedelta64(1000000000, 'ns'))*10**(-1)
mwd_ndbc = np.asarray(ndbc['mean_wave_dir'][okt,0,0])
wind_spd_ndbc = np.asarray(ndbc['wind_spd'][okt,0,0])
wind_dir_ndbc = np.asarray(ndbc['wind_dir'][okt,0,0])

#################################################################################

oklon = np.where(lon_ww3 >= lon_ndbc)[0][0]
oklat = np.where(lat_ww3 >= lat_ndbc)[0][0]

wind_spd_ww3_a2o2a_a2w2a = np.sqrt(uwind_ww3_a2o2a_a2w2a[:,oklat,oklon]**2 + vwind_ww3_a2o2a_a2w2a[:,oklat,oklon]**2) 
wind_spd_ww3_a2o_a2w2a = np.sqrt(uwind_ww3_a2o_a2w2a[:,oklat,oklon]**2 + vwind_ww3_a2o_a2w2a[:,oklat,oklon]**2) 

kw = dict(levels=np.arange(-32,33,2))
plt.figure()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track,lat_best_track,'.-',color='grey',label='Best Track')
plt.plot(lon_ndbc,lat_ndbc,'.-',color='red',label=ndbc.station)
plt.axis('scaled')
plt.ylim([0,45])
plt.xlim([-100,-10])
plt.legend()

fig,ax = plt.subplots(figsize=(7,4))
plt.plot(time_ww3,hs_ww3_a2o2a_a2w2a[:,oklat,oklon],'.-',color='indianred',label=list_exp[0])
plt.plot(time_ww3,hs_ww3_a2o_a2w2a[:,oklat,oklon],'o-',color='chocolate',label=list_exp[1])
plt.plot(time_ndbc,hs_ndbc,'.-',color='royalblue',label=ndbc.station)
plt.title('Significant Wave Height',fontsize=16)
plt.ylabel('m',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.plot(np.tile(datetime(2020, 8, 26, 12, 0),len(np.arange(11))),np.arange(11),'--k')

fig,ax = plt.subplots(figsize=(7,4))
plt.plot(time_ww3,fp_ww3_a2o2a_a2w2a[:,oklat,oklon],'.-',color='indianred',label=list_exp[0])
plt.plot(time_ww3,fp_ww3_a2o_a2w2a[:,oklat,oklon],'o-',color='chocolate',label=list_exp[1])
plt.plot(time_ndbc,fp_ndbc,'*',color='royalblue',label=ndbc.station)
plt.title('Dominant Wave Frequency',fontsize=16)
plt.ylabel('1/s',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.plot(np.tile(datetime(2020, 8, 26, 12, 0),len(np.arange(3))),np.arange(3),'--k')

fig,ax = plt.subplots(figsize=(7,4))
plt.plot(time_ww3,wind_spd_ww3_a2o2a_a2w2a,'.-',color='indianred',label=list_exp[0])
plt.plot(time_ww3,wind_spd_ww3_a2o2a_a2w2a,'o-',color='chocolate',label=list_exp[1])
plt.plot(time_ndbc,wind_spd_ndbc,'*',color='royalblue',label=ndbc.station)
plt.title('Wind Speed',fontsize=16)
plt.ylabel('m/s',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.plot(np.tile(datetime(2020, 8, 26, 12, 0),len(np.arange(30))),np.arange(30),'--k')

'''
u = uwind_ww3[:,oklat,oklon]
u[u==0] = np.nan
v = vwind_ww3[:,oklat,oklon]
v[v==0] = np.nan
V = [u,v]
wind_dir_ww3 = np.asarray([angle_0_360(V[0][i],V[1][i]) for i in np.arange(len(u))])
fig,ax = plt.subplots(figsize=(7,4))
plt.plot(time_ww3,wind_dir_ww3,'.-',color='indianred',label='WW3')
plt.plot(time_ndbc,wind_dir_ndbc,'*',color='royalblue',label=ndbc.station)
plt.title('Wind Direction',fontsize=16)
plt.ylabel('m/s',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.plot(np.tile(datetime(2020, 8, 26, 12, 0),len(np.arange(360))),np.arange(360),'--k')
'''



