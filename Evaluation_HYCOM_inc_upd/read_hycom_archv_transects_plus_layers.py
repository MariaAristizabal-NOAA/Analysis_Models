#%% User input

# RTOFS grid file name
hycom_depth = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/hafs_hycom_hat10.basin.regional.depth_no_pacif'
hycom_grid = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/regional.grid'
archv_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2_HYCOM_increm_update/natl00l.2020072812.hafs_hycom_hat10.2020_210_18'
#afile = 'archv_rtofs'
#afile = 'archv_inc_mom6'
#afile = 'incupd.2020_210_12'
#afile = 'diff_temp_mom6_minus_rtofs'
#afile = 'diff_saln_mom6_minus_rtofs'

var_name = 'temp'
#colormap='nipy_spectral'
#colormap='seismic'
colormap='Spectral_r'
min_val = 0
max_val = 33
delta_val = 1  # delta in colorbar
delta_contour = 2 # delta in contour plot
units = '$^oC$'
lon_trans = -60
lat_trans = 30

'''
var_name = 'salin'
klayer = '1'
colormap='GnBu_r'
#colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = 30
max_val = 40
delta_val = 0.5
delta_contour = 2
units = ' '
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]
'''

#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime
import os
import os.path
import glob
import numpy.ma as ma
import struct

import sys
sys.path.append('/home/Maria.Aristizabal/Repos/NCEP_scripts/')
from utils4HYCOM import readBinz, readgrids, parse_z
from utils4HYCOM import readdepth, readVar

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

##################################################################################
def get_profile_rtofs_ab_file_desn_layers(nz,oklon,oklat,archv_file):

    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz))
    target_temp_rtofs[:] = np.nan
    target_z_rtofs = np.empty((nz))
    target_z_rtofs[:] = np.nan

    lines = [line.rstrip() for line in open(archv_file+'.b')]

    ztmp = readVar(archv_file,'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklat,oklon]
    for lyr in tuple(layers):
        #print(lyr)
        temp_RTOFS = readVar(archv_file,'archive','temp',[lyr+1])
        target_temp_rtofs[lyr] = temp_RTOFS[oklat,oklon]
        dp = readVar(archv_file,'archive','thknss',[lyr+1])/2/9806
        target_ztmp = np.append(target_ztmp,dp[oklat,oklon])

    target_z3d = np.cumsum(target_ztmp)              # [idm,jdm,kdm+1]
    target_z3d = np.squeeze(target_z3d[1:])             # [idm,jdm,kdm]
    target_z3d = np.asarray(target_z3d)
    target_z3d[target_z3d > 10**8] = np.nan
    target_z_rtofs = target_z3d

    return target_temp_rtofs, target_z_rtofs

##################################################################################

#%% Reading hycom grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(hycom_grid,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid,'plat:',[0]))
depth_hycom = np.asarray(readdepth(hycom_depth,'depth'))

#%% Getting grid points for transects
oklon = int(np.interp(lon_trans +360,lon_hycom[0,:],np.arange(len(lon_hycom[0,:]))))
oklat = int(np.interp(lat_trans,lat_hycom[:,0],np.arange(len(lat_hycom[:,0]))))

####################################################################
lines = [line.rstrip() for line in open(archv_file+'.b')]
#ijdm = idm*jdm
#npad = 4096-(ijdm%4096)
#fld = ma.array([],fill_value=1.2676506002282294e+30)

for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == 'field':
            nheading = n + 1
            print(line.split()[0])

nvar = []
dens_layer = []
for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == var_name:
            nvar.append(n - nheading)
            dens_layer.append(line.split()[5])
nk = len(nvar)

trans_zonal = np.empty([nk,int(idm)])
trans_zonal[:] = np.nan
trans_zonal_depth = np.empty([nk,int(idm)])
trans_zonal_depth[:] = np.nan
for x in np.arange(idm):
    print(x)
    trans_zonal[:,x], trans_zonal_depth[:,x] = get_profile_rtofs_ab_file_desn_layers(nk,x,oklat,archv_file)


trans_merid = np.empty([nk,int(jdm)])
trans_merid[:] = np.nan
trans_merid_depth = np.empty([nk,int(idm)])
trans_merid_depth[:] = np.nan
for y in np.arange(jdm):
    print(y)
    trans_zonal[:,y], trans_zonal_depth[:,y] = get_profile_rtofs_ab_file_desn_layers(nk,oklon,y,archv_file)

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax1 = plt.subplots(figsize = (11,6))
ax1.set_facecolor("bisque")
plt.contour(lon_hycom[0,:],-trans_zonal_depths,dens_layers,np.arange(min_val,max_val,delta_contour),colors='grey')
#plt.contourf(lon_hycom[0,:],-depths,np.flip(trans_zonal,axis=0),cmap=colormap ) #,**kw)
plt.contourf(lon_hycom[0,:],-trans_zonal_depth,trans_zonal,cmap=colormap,**kw)
#plt.axis('scaled')
cbar = plt.colorbar(fraction=0.025, pad=0.04)
cbar.set_label(units,fontsize = 14)
plt.title(var_name + ' Zonal Transect ',fontsize=18)
#plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
#plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
#plt.xlim(lon_lim)
#plt.ylim(lat_lim)

####################################################################

