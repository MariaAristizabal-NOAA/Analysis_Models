#%% User input

# forecasting cycle to be used
cycle = '20220902'

# Transect lon and lat limits
lon_lim = [-44,-44]
lat_lim = [34,42]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

folder_rtofs = scratch_folder + 'RTOFS2022_oper/rtofs.' + cycle +'/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# folder utils for Hycom
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = scratch_folder + 'figures/'

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
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Time window
'''
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]
tini = datetime.strptime(date_ini,'%Y/%m/%d')
tend = tini + timedelta(hours=192)
date_end = tend.strftime('%Y/%m/%d')
'''

#################################################################################
#%% Reading bathymetry data
'''
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
'''

#################################################################################
#%% Get list files
files_rtofs = sorted(glob.glob(folder_rtofs +'*.nc'))
#files_hafs_hycom = sorted(glob.glob(os.path.join(folder_hafs,'*3z*.nc')))

################################################################################
#%% Reading HAFS/HYCOM grid
rtofs_grid = xr.open_dataset(files_rtofs[0],decode_times=False)
lon_rtofs = np.asarray(rtofs_grid['Longitude'][:])
lat_rtofs = np.asarray(rtofs_grid['Latitude'][:])
depth_rtofs = np.asarray(rtofs_grid['Depth'][:])

#################################################################################
#%% Read rtofs_da time
time_rtofs = []
timestamp_rtofs = []
for n,file in enumerate(files_rtofs):
    print(file)
    RTOFS = xr.open_dataset(file)
    t = RTOFS['MT'][:]
    timestamp = mdates.date2num(t)[0]
    time_rtofs.append(mdates.num2date(timestamp))
    timestamp_rtofs.append(timestamp)

time_rtofss = np.asarray(time_rtofs)
okt = np.argsort(time_rtofss)
time_rtofs = time_rtofss[okt]
timestamp_rtofs = np.asarray(timestamp_rtofs)[okt]

files_rtofs = np.asarray(files_rtofs)[okt]

##################################################################
#%% Longitudinal transect
lon = lon_rtofs
lat = lat_rtofs

xlim = lon_lim
ylim = lat_lim

#xmin = int(np.round(np.interp(xlim[0],lon[0,:],np.arange(len(lon[0,:])))))
#xmax = int(np.round(np.interp(xlim[1],lon[0,:],np.arange(len(lon[0,:])))))
xmin = np.where(lon[0,:] < xlim[0])[0][-1]
xmax = np.where(lon[0,:] < xlim[1])[0][-1]
ymin = int(np.round(np.interp(ylim[0],lat[:,0],np.arange(len(lat[:,0])))))
ymax = int(np.round(np.interp(ylim[1],lat[:,0],np.arange(len(lat[:,0])))))
if xmin == xmax:
    xmax = xmax + 1

    latt = lat[ymin:ymax,:][:,xmin:xmax][:,0]

    MODEL0 = xr.open_dataset(files_rtofs[0])
    temp0 = np.asarray(MODEL0['pot_temp'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
    temp0[temp0 == 0] = np.nan

    for n,file in enumerate(files_rtofs):
        print(file)
        f = file.split('/')[-1].split('.')[-3]

        rtofs = xr.open_dataset(file)

        MODEL = xr.open_dataset(file)
        temp = np.asarray(MODEL['pot_temp'][0,:,ymin:ymax,:][:,:,xmin:xmax])[:,:,0]
        temp[temp == 0] = np.nan

        print(np.nanmin(temp-temp0))
        #############################################################
        # Temp

        kw = dict(levels=np.arange(15,31.1,0.5))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(latt,-depth_rtofs,temp,cmap='Spectral_r',**kw,extend='both')

        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(latt,-depth_rtofs,temp,[26],colors='k')
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_title('RTOFS Forecast for ' + ' Init: ' + cycle + ' ' + f + '\n' + 'Temperature')
        fname = cycle + '.rtofs.temp_trans.' + f + '.png'
        fig.savefig(fname,bbox_inches='tight')
        #plt.close()
        #####################################################################
        # Figure temp-temp_f000

        kw = dict(levels=np.arange(-4,4.1,0.2))
        fig,ax = plt.subplots(figsize=(8,4))
        ctr = ax.contourf(latt,-depth_rtofs,temp-temp0,cmap='seismic',**kw,extend='both')
        cbar = fig.colorbar(ctr,extendrect=True)
        cbar.set_label('$^oC$',fontsize=14)
        cs = ax.contour(latt,-depth_rtofs,temp-temp0,[0],colors='k',alpha=0.3)
        ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
        ax.set_ylim([-300,0])
        ax.set_title('RTOFS Forecast for ' + ' Init: ' + cycle + ' ' + f + '\n' + 'Temperature difference')
        ax.set_xlabel('dt min = ' + str(np.round(np.nanmin(temp-temp0),2)),fontsize=14)
        fname = cycle + '.pom.temp_diff_trans.f' + f + '.png'
        fig.savefig(fname,bbox_inches='tight')
        #plt.close()

    ############################################################################




################################################################################
'''
#%% Figure SST all domain HAFS
n = 3
RTOFS_DA0 = xr.open_dataset(files_rtofs_da[n])
SST_rtofs_da = np.asarray(RTOFS_DA0['temperature'][0,0,:,:])

if cycle == '20211016':
    kw = dict(levels=np.arange(-5,32.1,0.5))
    lev = [22,23,24,25,26,27,28,29,30,31]

if cycle == '2021082800' or '2021093006':
    kw = dict(levels=np.arange(22,32.5,0.5))
    lev = [22,23,24,25,26,27,28,29,30,31,32]

plt.figure()
#plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
#plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contour(lon_rtofs_da,lat_rtofs_da,SST_rtofs_da,lev,colors='grey',alpha=0.5)
plt.contourf(lon_rtofs_da,lat_rtofs_da,SST_rtofs_da,cmap='Spectral_r',**kw) #,vmin=-3.0,vmax=3.0)
plt.colorbar()
plt.axis('scaled')
#plt.ylim([2,40])
#plt.xlim([-98,-40])
plt.title('SST RTOFS-DA '+ str(time_rtofs_da[n])[0:13])
plt.text(-140,20,file.split('/')[-1].split('_')[3],fontsize=14)
#plt.legend()
file_name = folder_fig + 'SST_rtofs_da_' + file.split('/')[-1].split('.')[-2]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
'''
################################################################################
'''
#%% find grid point
target_lat = 10.18
target_lon = -119.22

# search in xi_rho direction
oklatmm = []
oklonmm = []

lat = lat_rtofs_da
lon = lon_rtofs_da

for pos_xi in np.arange(lat.shape[1]):
    pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),\
                                             left=np.nan,right=np.nan))
    if np.isfinite(pos_eta):
        oklatmm.append((pos_eta).astype(int))
        oklonmm.append(pos_xi)

pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
oklat = oklatmm[pos]
oklon = oklonmm[pos]
'''
'''
#search in eta-rho direction
oklatmm = []
oklonmm = []
for pos_eta in np.arange(lat.shape[0]):
    pos_xi = np.round(np.interp(target_lon,lon[pos_eta,:],np.arange(len(lon[pos_eta,:])),\
                                            left=np.nan,right=np.nan))
    if np.isfinite(pos_xi):
        oklatmm.append(pos_eta)
        oklonmm.append(pos_xi.astype(int))

pos_lat = np.round(np.interp(target_lat,lat[oklatmm,oklonmm],np.arange(len(lat[oklatmm,oklonmm])))).astype(int)
oklat2 = oklatmm[pos_lat]
oklon2 = oklonmm[pos_lat]

#check for minimum distance
dist1 = np.sqrt((oklon1-target_lon)**2 + (oklat1-target_lat)**2)
dist2 = np.sqrt((oklon2-target_lon)**2 + (oklat2-target_lat)**2)
if dist1 >= dist2:
    oklat = oklat1
    oklon = oklon1
else:
    oklat = oklat2
    oklon= oklon2

oklat = oklat.astype(int)
oklon = oklon.astype(int)
'''

