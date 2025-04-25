#%% User input

# folder ab files HYCOM
folder_hycom1 = '/mnt/lfs4/HFIP/hwrfv3/Maria.Aristizabal/hafstmp/hafsv0.3_phase2_hycom_mvnest/2021082400/16W/ocn_prep/hycominit1/'
folder_hycom2 = '/mnt/lfs4/HFIP/hwrfv3/Maria.Aristizabal/hafsv0.3_phase2_hycom/fix/fix_hycom/'


# RTOFS grid file name
hycom_grid_global = folder_hycom1 + 'regional.grid'
hycom_depth_global = folder_hycom1 + 'regional.depth'
hycom_restart_global = folder_hycom1 + 'rtofs_glo.t00z.n00.restart'

hycom_grid_reg = folder_hycom2 + 'hafs_hycom_jtnh.basin.regional.grid'
hycom_depth_reg = folder_hycom2 + 'hafs_hycom_jtnh.basin.regional.depth'
hycom_restart_reg = folder_hycom1 + 'restart_out'

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
from utils4HYCOM import readdepth

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Reading hycom grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid_global+'.b')]
idm_global = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm_global = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom_global = np.array(readgrids(hycom_grid_global,'plon:',[0]))
lat_hycom_global = np.array(readgrids(hycom_grid_global,'plat:',[0]))
depth_hycom_glob = np.asarray(readdepth(hycom_depth_global,'depth'))

lines_grid = [line.rstrip() for line in open(hycom_grid_reg+'.b')]
idm_reg = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm_reg = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom_reg = np.array(readgrids(hycom_grid_reg,'plon:',[0]))
lat_hycom_reg = np.array(readgrids(hycom_grid_reg,'plat:',[0]))
depth_hycom_reg = np.asarray(readdepth(hycom_depth_reg,'depth'))

####################################################################
lines_restart = [line.rstrip() for line in open(hycom_restart_global+'.b')]
ijdm = idm_global*jdm_global
npad = 4096-(ijdm%4096)
fld2 = ma.array([],fill_value=1.2676506002282294e+30)
#fld2 = ma.array([],fill_value=1e30)

inFile = hycom_restart_global + '.a'

nvar = 248
fid = open(inFile,'rb')
fid.seek((nvar-1)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld_global = ma.reshape(fld,(jdm_global,idm_global))

print(np.max(fld_global))
print(np.min(fld_global))

kw = dict(levels=np.arange(-7,40))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom_global,lat_hycom_global,fld_global,cmap='Spectral_r',**kw)
plt.axis('scaled')
plt.colorbar()

####################################################################
lines_restart = [line.rstrip() for line in open(hycom_restart_reg+'.b')]
ijdm = idm_reg*jdm_reg
npad = 4096-(ijdm%4096)
fld2 = ma.array([],fill_value=1.2676506002282294e+30)
#fld2 = ma.array([],fill_value=1e30)

inFile = hycom_restart_reg + '.a'

nvar = 248
fid = open(inFile,'rb')
fid.seek((nvar-1)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld_reg = ma.reshape(fld,(jdm_reg,idm_reg))

print(np.max(fld_reg))
print(np.min(fld_reg))

kw = dict(levels=np.arange(-7,40))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom_reg,lat_hycom_reg,fld_reg,cmap='Spectral_r',**kw)
plt.axis('scaled')
plt.colorbar()

kw = dict(levels=np.arange(0,8000,500))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom_reg,lat_hycom_reg,depth_hycom_reg,cmap='Spectral_r',**kw)
plt.axis('scaled')
plt.colorbar()

