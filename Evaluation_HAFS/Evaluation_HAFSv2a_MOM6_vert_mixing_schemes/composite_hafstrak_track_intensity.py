#%% User input

# forecasting cycle to be used

# Lee
storm_num = '13'
basin = 'al'
storm_id = '13l'
year = '2023'

# Fiona
#storm_num = '07'
#basin = 'al'
#storm_id = '07l'
#year = '2022'

# Phillipe
#storm_num = '17'
#basin = 'al'
#storm_id = '17l'
#year = '2023'

exp_names = ['hafsv2a_final','hafsv2_h3a0_rtofsv2','hafsv2_h3a0_rtofsv2_mom6epbl']
exp_labels = ['HV2A','H2AR','HEBL']
exp_colors = ['cyan','deeppink','lawngreen']


#exp_names = ['hafsv2a_final','hafsv2_h3a0_rtofsv2_mom6ri','hafsv2_h3a0_rtofsv2_mom6epbl']
#exp_labels = ['HV2A','HR25','HEBL']
#exp_colors = ['cyan','blue','lawngreen']

#exp_names = ['hafsv2_20231212_h2ab','hafsv2a_hkpp','hafsv2a_hkp2']
#exp_labels = ['hafsv2_baseline','hafsv2a_hkpp','hafsv2a_hkp2']
#exp_colors = ['orange','forestgreen','royalblue']

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

lat_lim = [-100,-60]
lon_lim = [10,55]

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
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

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder1 + 'hafstrak/' + exp_names[i] + '/')

###################################################################
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

###################################################################
#%% Read best track
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + year + '.dat'

lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

###################################################################
#%% Loop the experiments

# Find all cycles
Cycles = []
for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get storm track from trak atcf files
    if exp_labels[i] == 'HV2A' or exp_labels[i] == 'H2AR' or exp_names[i] == 'HEBL':
        file_tracks = np.sort(glob.glob(os.path.join(folder,storm_id+'*.hfsa.trak.atcfunix')))
    if exp_names[i] == 'HFSB':
        file_tracks = np.sort(glob.glob(os.path.join(folder,storm_id+'*.hfsb.trak.atcfunix')))
    
    print(file_tracks)
    cycles = [f.split('/')[-1].split('.')[1] for f in file_tracks]
    print(cycles)
    Cycles.append(cycles)

all_cycles = []
for cyc in Cycles:
    for cy in cyc:
        print(cy)
        if cy not in all_cycles:
            print(cy)
            all_cycles.append(cy)

num_cycles = len(all_cycles)

#%% Loop the experiments
lon_forec_track = np.empty((len(folder_exps),num_cycles,43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),num_cycles,43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),num_cycles,43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),num_cycles,43))
int_track[:] = np.nan
for i,folder in enumerate(folder_exps):
    print(folder)
    for c,cycle in enumerate(all_cycles):
        print(cycle)
        #%% Get storm track from trak atcf files
        if exp_names[i] == 'hafsv2a_final' or exp_names[i] == 'hafsv2_h3a0_rtofsv2' or exp_names[i] == 'hafsv2_h3a0_rtofsv2_mom6epbl':
            file_track = folder + storm_id + '.' + cycle + '.hfsa.trak.atcfunix'
        if exp_names[i] == 'HFSB':
            file_track = glob.glob(os.path.join(folder,storm_id+'*.hfsb.trak.atcfunix'))[0]
    
        print(file_tracks)

        if os.path.isfile(file_track):
            okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
            lon_forec_track[i,c,0:okn], lat_forec_track[i,c,0:okn], lead_time[i,c,0:okn], int_track[i,c,0:okn], _ = get_storm_track_and_int(file_track,storm_num)
        else:
            lon_forec_track[i,c,:] = np.nan
            lat_forec_track[i,c,:] = np.nan
            lead_time[i,c,:] = np.nan
            int_track[i,c,:] = np.nan
 
###################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

for f,folder in enumerate(folder_exps):
    fig,ax = plt.subplots(figsize=(8, 4))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.plot(lon_forec_track[f,0,::2], lat_forec_track[f,0,::2],'-',color=exp_colors[f],label=exp_labels[f])
    for c in np.arange(len(all_cycles)):
        plt.plot(lon_forec_track[f,c,::2], lat_forec_track[f,c,::2],'-',color=exp_colors[f])
    plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
    plt.legend(loc='upper right')
    #plt.legend(loc='upper right',bbox_to_anchor=[1.4,1.0])
    plt.title('Track Forecast ' + storm_num ,fontsize=18)
    plt.axis('scaled')
    plt.xlim([np.min(lon_best_track)-5,np.max(lon_best_track)+5])
    plt.ylim([np.min(lat_best_track)-5,np.max(lat_best_track)+5])
    #plt.savefig('track'+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

########################################################################
#%% Figure intensity

for f,folder in enumerate(folder_exps):
    fig,ax1 = plt.subplots(figsize=(7, 4))
    plt.plot(time_best_track,int_best_track,'o-k',label='Best')
    #plt.plot(lead_time[f,cc,::2][oktt],int_best_track[okt],'o-k',label='Best')
    for c,cycle in enumerate(all_cycles):
        #%% Time window
        date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
        tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
        tend = tini + timedelta(hours=126)
        date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')
        date_vec = [tini + timedelta(hours=3*int(dt)) for dt in np.arange(0,127,3)]

        cc = [cc for cc,cyc in enumerate(all_cycles) if cyc == cycle][0]
        #plt.plot(lead_time[f,cc,::2],int_track[f,cc,::2],'o-',color=exp_colors[f],label=exp_labels[f],markeredgecolor='k',markersize=7)
        plt.plot(date_vec[::2],int_track[f,cc,::2],'-',color=exp_colors[f])
        if c == 0:
            plt.plot(date_vec[::2],int_track[f,0,::2],'-',color=exp_colors[f],label=exp_labels[f])

    plt.legend(loc='upper right')

    ax1.tick_params(which='major', width=2)
    ax1.tick_params(which='major', length=7)
    #ax1.tick_params(which='minor', length=4, color='k')

    #ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    #ax1.xaxis.set_minor_locator(MultipleLocator(3))
    #ax1.xaxis.set_ticks(np.arange(0,126,12))
    #plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
    plt.ylim([10,190])
    #plt.xlim([0,126])
    #plt.xticks(np.arange(0,126,12))
    plt.title('Intensity Forecast',fontsize=18)
    plt.ylabel('Max 10m Wind (kt)',fontsize=14)

    ax2 = ax1.twinx()
    plt.ylim([10,190])
    yticks = [64,83,96,113,137]
    plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
    plt.grid(True)
        
################################################################################
