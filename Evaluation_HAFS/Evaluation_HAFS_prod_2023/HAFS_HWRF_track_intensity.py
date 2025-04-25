#%% User input

# forecasting cycle to be used

# Lee
cycle = '2023090706'
#cycle = '2023090718'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'

exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_MOM6_epbl','HWRF_HYCOM']
exp_labels = ['HFSA-HYCOM (oper)','HWRF-POM (oper)','HFSA-MOM6 (exp)','HWRF-HYCOM (exp)']
exp_colors = ['darkviolet','pink','forestgreen','royalblue']

#exp_names = ['HFSA_oper','HWRF_2023','HFSAv1p1_HYCOM','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_labels = ['HFSA_oper_HYCOM_epbl','HWRF_POM_MY2.5','HFSAv1p1_HYCOM_kpp','HFSAv1p1_MOM6_epbl','HFSAv1p1_MOM6_kpp']
#exp_colors = ['darkviolet','pink','forestgreen','cyan','royalblue']

#exp_names = ['HFSA_oper','HWRF_2023','HMON_2023','HFSAv1p1']
#exp_labels = ['HFSA','HWRF','HMON','HFSAv1p1']
#exp_colors = ['darkviolet','pink','c','steelblue']

lon_lim = [-80,-55]
lat_lim = [10.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

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

#########################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder1 + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

###################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

###################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

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
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

###################################################################
#%% Loop the experiments
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan
for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get storm track from trak atcf files
    if exp_labels[i] == 'HMON':
        file_track = glob.glob(os.path.join(folder,'*'+storm_id+'.'+cycle+'*trak.hmon.atcfunix'))[0]
        print(file_track)
    if exp_names[i] == 'HFSA_oper' or exp_names[i] == 'HFSAv1p1_HYCOM' or exp_names[i] == 'HFSAv1p1_MOM6_epbl' or exp_names[i] == 'HFSAv1p1_MOM6_kpp':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
        print(file_track)
    if exp_names[i] == 'HFSB':
        file_track = folder + + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
        print(file_track)
    if exp_names[i] == 'HWRF_2023' or exp_names[i] == 'HWRF_HYCOM':
        file_track = folder + storm_name + storm_id + '.' + cycle + '.trak.hwrf.atcfunix'
        print(file_track)

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

###################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots(figsize=(8, 4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
#plt.legend(loc='upper right')
plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
#plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.min(lon_best_track[okt])-5,np.max(lon_best_track[okt])+5])
plt.ylim([np.min(lat_best_track[okt])-5,np.max(lat_best_track[okt])+5])
plt.savefig('track'+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

########################################################################
#%% Figure intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)
tt_fv3 = np.asarray([t.replace(tzinfo=None) for t in time_fv3])
oktt = np.logical_and(tt_fv3[::2] >= time_best_track[0],tt_fv3[::2] <= time_best_track[-1])

fig,ax1 = plt.subplots(figsize=(7, 4))
for ii,i in enumerate(np.arange(len(exp_names))): 
    plt.plot(lead_time[i,::2],int_track[i,::2],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)

plt.plot(lead_time[i,::2][oktt],int_best_track[okt],'o-k',label='Best')
#plt.legend(loc='upper right',bbox_to_anchor=[1.1,1.15])
#plt.legend(loc='lower right')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
#ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.ylim([10,190])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Cycle '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([10,190])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)
plt.savefig('intensity'+'.png',format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
