#%% User input
cycle = '20230907'

# folder ab files HYCOM
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
rtofs_folder = scratch_folder + 'RTOFS_2023/20230907/'

#rtofs_folder_par3 = scratch_folder + 'RTOFS_DA/parallel3/' + 'rtofs.' + cycle + '_ab/'
#rtofs_folder_run4 = scratch_folder + 'RTOFS_DA/rtofs_run4/' + 'rtofs.' + cycle + '/'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

lon_lim = [-55,-55]
lat_lim = [15,25]


##################################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
def find_grid_position_hycom(lon,lat,target_lon,target_lat):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),left=np.nan,right=np.nan))
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
    timeRT = tzero+timedelta(float(hycom_days))
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
def get_latitudinal_transect_rtofs_ab_file_desn_layers(nz,oklon,oklats,afile,bfile):

    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz,len(oklats)))
    target_temp_rtofs[:] = np.nan
    target_z_2d = np.empty((nz+1,len(oklats)))
    target_z_2d[:] = np.nan
    timeRTOFS = []

    lines = [line.rstrip() for line in open(bfile)]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days))
    timeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT)

    ztmp = readVar(afile[:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklats,oklon]
    target_z_2d[0,:] = target_ztmp
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(afile[:-2],'archive','temp',[lyr+1])
        target_temp_rtofs[lyr,:] = temp_RTOFS[oklats,oklon]
        dp = readVar(afile[:-2],'archive','thknss',[lyr+1])/2/9806
        target_z_2d[lyr+1,:] = dp[oklats,oklon]

    target_z3d = np.cumsum(target_z_2d,0)              # [idm,jdm,kdm+1]
    target_z3d = target_z3d[1:]             # [idm,jdm,kdm]
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
artofs_files = sorted(glob.glob(os.path.join(rtofs_folder,'*archv.a')))

##################################################################################
#%% Reading HYCOM grid from ab files
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(rtofs_grid_file+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth_hycom = np.asarray(readdepth(rtofs_depth_file,'depth'))

#######################################################################
#%% find grid point
lon_limH,lat_limH = geo_coord_to_HYCOM_coord(lon_lim,lat_lim)

lat = lat_hycom
lon = lon_hycom

okloni, oklati = find_grid_position_hycom(lon,lat,lon_limH[0],lat_limH[0]) 
oklonf, oklatf = find_grid_position_hycom(lon,lat,lon_limH[1],lat_limH[1]) 

#######################################################################
#%% get target temp latitudinal profile from ab files
nz = 41
oklats = np.arange(oklati,oklatf)
oklon = okloni

# reference rtofs file
#n0 = 4
n0 = -1
afile = artofs_files[n0]
bfile = afile[:-1] + 'b'

target_temp_rtofs0, target_z_rtofs0, time_rtofs0 = get_latitudinal_transect_rtofs_ab_file_desn_layers(nz,oklon,oklats,afile,bfile)

for afile in artofs_files:
    bfile = afile[:-1] + 'b'

    target_temp_rtofs, target_z_rtofs, time_rtofs = get_latitudinal_transect_rtofs_ab_file_desn_layers(nz,oklon,oklats,afile,bfile)

    ###################################################################    #%% Plot transect

    fhour = bfile.split('/')[-1].split('.')[2]
    target_lat_rtofs = np.tile(lat_hycom[oklats,oklon],(nz,1))

    kw = dict(levels=np.arange(12,31.1,0.5))
    fig,ax = plt.subplots(figsize=(8,4))
    ctr = ax.contourf(target_lat_rtofs,-target_z_rtofs,target_temp_rtofs,cmap='Spectral_r',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(target_lat_rtofs,-target_z_rtofs,target_temp_rtofs,[26],colors='k')
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.set_ylim([-300,0])
    ax.set_title('RTOFS Forecast for ' + ' Cycle: ' + cycle + ' ' + fhour + ' Valid time ' +  str(time_rtofs[0])[0:13] + '\n' + 'Temperature')
    fname = cycle + '.rtofs.temp_trans.' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()

##################################################################################
