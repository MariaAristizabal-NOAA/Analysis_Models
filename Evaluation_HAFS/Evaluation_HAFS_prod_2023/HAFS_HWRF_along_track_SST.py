#%% User input

# forecasting cycle to be used
# Lee
cycle = '2023090706'
storm_id = '13l'
fhour = 'f006'
storm_name = 'lee'
basin = 'al'
storm_num = '13'

exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1']
exp_labels = ['HFSA','HWRF','HFSAv1p1']
exp_colors = ['darkviolet','pink','steelblue']

# Transect lon and lat limits
lon_lim = [-59.8,-59.8]
lat_lim = [15.8,22.4]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder + exp_names[0] + '/' + cycle + '/' + storm_id + '/',scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_id + '/',scratch_folder + exp_names[2] + '/' + cycle + '/' + storm_id + '/']

best_track_file = abdeck_folder + 'btk/b' + basin + storm_id + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int

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

#################################################################################
#%% Loop the experiments

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
target_sst = np.empty((len(folder_exps),100))
target_sst[:] = np.nan
target_sst0 = np.empty((len(folder_exps),100))
target_sst0[:] = np.nan
target_lat = np.empty((len(folder_exps),100))
target_lat[:] = np.nan
for i,folder in enumerate(folder_exps):
    print(folder)    

    if exp_labels[i] == 'HFSA':
        file0 = folder + storm_id + '.' + cycle + '.hfsa.hycom.3z.f000.nc'
        file = folder + storm_id + '.' + cycle + '.hfsa.hycom.3z.' + fhour + '.nc'
        file_track = glob.glob(os.path.join(folder,'*hfsa.trak.atcfunix'))[0]

        #okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
        #lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], _, _, _ = get_storm_track_and_int(file_track,storm_num)
        lon_forec_track, lat_forec_track, lead_time, int_track, _ = get_storm_track_and_int(file_track,storm_id[0:-1])

        #%% Read HAFS/HYCOM time
        MODEL0 = xr.open_dataset(file0)
        MODEL = xr.open_dataset(file)
        temp0 = MODEL0['temperature']
        temp = MODEL['temperature']
        t = MODEL['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading HAFS/HYCOM grid)
        lon = np.asarray(MODEL['Longitude'][:])-360
        lat = np.asarray(MODEL['Latitude'][:])
        depth = np.asarray(MODEL['Z'][:])

    if exp_labels[i] == 'HFSAv1p1':
        file003 = folder + storm_id + '.' + cycle + '.hfsa.mom6.f003.nc'
        file = folder + storm_id + '.' + cycle + '.hfsa.mom6.' + fhour + '.nc'
        file_track = os.path.join(folder,storm_id+'.'+cycle+'.hfsa.trak.atcfunix')

        #okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
        #lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], _, _, _ = get_storm_track_and_int(file_track,storm_num)
        lon_forec_track, lat_forec_track, lead_time, int_track, _ = get_storm_track_and_int(file_track,storm_id[0:-1])
    
        #%% Read HAFS/HYCOM time
        MODEL0 = xr.open_dataset(file003)
        MODEL = xr.open_dataset(file)
        temp0 = MODEL0['temp']
        temp = MODEL['temp']
        t = MODEL['time'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading grid
        lon = np.asarray(MODEL['xh'][:]) 
        lat = np.asarray(MODEL['yh'][:])
    
    if exp_labels[i] == 'HWRF':
        fnu = int(int(fhour[1:])/6)
        if fnu <= 10:
            fnum = '0' + str(fnu)
        else:
            fnum = str(fnu)  
        file0 = folder + storm_name + storm_id + '.' + cycle + '.pom.0000.nc'
        file = folder + storm_name + storm_id + '.' + cycle + '.pom.00' + fnum + '.nc'
        grid = folder + storm_name + storm_id + '.' + cycle + '.pom.grid.nc'
        file_track = glob.glob(os.path.join(folder,'*'+storm_id+'.'+cycle+'*trak.hwrf.atcfunix'))[0]

        #okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
        #lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], _, _, _ = get_storm_track_and_int(file_track,storm_num)
        lon_forec_track, lat_forec_track, lead_time, int_track, _ = get_storm_track_and_int(file_track,storm_id[0:-1])

        #%% Read time
        MODEL0 = xr.open_dataset(file0)
        MODEL = xr.open_dataset(file)
        temp0 = MODEL0['t']
        temp = MODEL['t']
        t = MODEL['time'][:]
        timestamp = mdates.date2num(t)[0]
        time = mdates.num2date(timestamp)

        #%% Reading grid
        GRID = xr.open_dataset(grid)
        lonp = np.asarray(GRID['east_e'][:])
        latp = np.asarray(GRID['north_e'][:])
        zlev_pom = np.asarray(GRID['zz'][:])
        hpom = np.asarray(GRID['h'][:])
        zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
        zmatrix = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])     
        lon = lonp[0,:]
        lat = latp[:,0]

    #%% Along-track transect
    #lon_forec_track_interp = np.interp(lat,lat_forec_track[i,:],lon_forec_track[i,:],left=np.nan,right=np.nan)
    lon_forec_track_interp = np.interp(lat,lat_forec_track,lon_forec_track,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat)
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan 
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    oklon = np.round(np.interp(lon_forec_track_int,lon,np.arange(len(lon)))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat,np.arange(len(lat)))).astype(int)
               
    for x in np.arange(len(lon_forec_track_int)):
        #target_sst0[i,x] = np.asarray(MODEL0['temperature'][0,0,oklat[x],oklon[x]]) # HFSA
        #target_sst[i,x] = np.asarray(MODEL['temperature'][0,0,oklat[x],oklon[x]])
        target_sst0[i,x] = np.asarray(temp0[0,0,oklat[x],oklon[x]])
        target_sst[i,x] = np.asarray(temp[0,0,oklat[x],oklon[x]])

    target_lat[i,0:len(lat[oklat])] = lat[oklat]

fig,ax = plt.subplots(figsize=(10,5))
for i,label in enumerate(exp_labels):
    ax.plot(target_lat[i,:],target_sst[i,:],'o-',color = exp_colors[i],markeredgecolor='k',label=label)
ax.legend()
ax.set_ylabel('SST',fontsize=14)
ax.set_xlabel('Latitude',fontsize=14)
plt.savefig('along_track_sst.png')

