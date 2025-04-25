#%% User input
cycle = '20211016'

target_lat = 3.77
target_lon = -95.56

# folder ab files HYCOM
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
rtofs_folder_par3 = scratch_folder + 'RTOFS_DA/parallel3/' + 'rtofs.' + cycle + '_ab/'
rtofs_folder_run4 = scratch_folder + 'RTOFS_DA/rtofs_run4/' + 'rtofs.' + cycle + '/'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_DA/GRID_DEPTH/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

##################################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
def find_grid_position_hycom(lon,lat,target_lon,target_lat):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),\
                                             left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklon,oklat

##################################################################################
def get_profile_rtofs_ab_file_desn_layers(nz,oklon,oklat,afile,bfile):

    #nz = 41
    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz))
    target_temp_rtofs[:] = np.nan
    target_z_rtofs = np.empty((nz))
    target_z_rtofs[:] = np.nan
    timeRTOFS = []

    lines = [line.rstrip() for line in open(bfile)]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days)-1)
    timeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT)

    ztmp = readVar(afile[:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklat,oklon]
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(afile[:-2],'archive','temp',[lyr+1])
        target_temp_rtofs[lyr] = temp_RTOFS[oklat,oklon]
        dp = readVar(afile[:-2],'archive','thknss',[lyr+1])/2/9806
        target_ztmp = np.append(target_ztmp,dp[oklat,oklon])

    target_z3d = np.cumsum(target_ztmp)              # [idm,jdm,kdm+1]
    target_z3d = np.squeeze(target_z3d[1:])             # [idm,jdm,kdm]
    target_z3d = np.asarray(target_z3d)
    target_z3d[target_z3d > 10**8] = np.nan
    target_z_rtofs = target_z3d

    time_rtofs = np.asarray(timeRTOFS)

    return target_temp_rtofs, target_z_rtofs, time_rtofs

##################################################################################
#%% 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import glob
import sys

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

sys.path.append(folder_myutils)
from my_models_utils import geo_coord_to_HYCOM_coord

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

##################################################################################
#%% finding list of hycom files
artofs_files_par3 = sorted(glob.glob(os.path.join(rtofs_folder_par3,'*n00.archv.a')))
brtofs_files_par3 = sorted(glob.glob(os.path.join(rtofs_folder_par3,'*n00.archv.b')))

artofs_files_run4 = sorted(glob.glob(os.path.join(rtofs_folder_run4,'*n00.archv.a')))
brtofs_files_run4 = sorted(glob.glob(os.path.join(rtofs_folder_run4,'*n00.archv.b')))
#artofs_files_run4 = sorted(glob.glob(os.path.join(rtofs_folder_run4,'*f06.archv.a')))
#brtofs_files_run4 = sorted(glob.glob(os.path.join(rtofs_folder_run4,'*f06.archv.b')))

##################################################################################
#%% Reading HYCOM grid from ab files
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(rtofs_grid_file+'.b')]
lon_rtofs_ab = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
lat_rtofs_ab = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth_rtofs_ab = np.asarray(readdepth(rtofs_depth_file,'depth'))

##################################################################################
#%% find grid point
target_lonH,target_latH = geo_coord_to_HYCOM_coord(target_lon,target_lat)

lat = lat_rtofs_ab
lon = lon_rtofs_ab

oklon_ab, oklat_ab = find_grid_position_hycom(lon,lat,target_lonH,target_latH) 

##################################################################################
#%% get target temp prof from ab files
nz = 41
target_temp_rtofs_run4, target_z_rtofs_run4, time_rtofs_run4 = get_profile_rtofs_ab_file_desn_layers(nz,oklon_ab,oklat_ab,artofs_files_run4[0],brtofs_files_run4[0])

target_temp_rtofs_par3, target_z_rtofs_par3, time_rtofs_par3 = get_profile_rtofs_ab_file_desn_layers(nz,oklon_ab,oklat_ab,artofs_files_par3[0],brtofs_files_par3[0])

##################################################################################
#%% Plot profile

plt.figure(figsize=(5,7))
plt.plot(target_temp_rtofs_par3,-target_z_rtofs_par3,'.-',label = 'Parallel3 '+str(time_rtofs_par3[0])[0:13])
plt.plot(target_temp_rtofs_run4,-target_z_rtofs_run4,'.-',label = 'Run4 '+str(time_rtofs_run4[0])[0:13])
plt.ylim([-1800,0])
plt.title('Temperature RTOFS '+ '\n Position ' + '(' + str(np.round(lon_rtofs_ab[oklat_ab,oklon_ab]-360,4)) + ',' + str(np.round(lat_rtofs_ab[oklat_ab,oklon_ab],4)) + ')')
plt.legend()

##################################################################################
