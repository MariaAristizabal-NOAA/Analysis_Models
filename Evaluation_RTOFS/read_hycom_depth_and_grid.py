#%% User input

#depthfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_nhc.basin.regional.depth'
#gridfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_nhc.basin.regional.grid'
#depthfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_jtnh.basin.regional.depth'
#gridfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_jtnh.basin.regional.grid'
depthfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_jtsh.basin.regional.depth'
gridfile_hycom = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_hycom/hafs_hycom_jtsh.basin.regional.grid'

#%%
import matplotlib.pyplot as plt
import numpy as np
import sys

import sys
sys.path.append('/home/maristiz/Repos/NCEP_scripts/')
from utils4HYCOM import readBinz, readgrids, parse_z
from utils4HYCOM import readdepth

#%% Reading hycom grid
lines_grid = [line.rstrip() for line in open(gridfile_hycom+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(gridfile_hycom,'plon:',[0]))
lat_hycom = np.array(readgrids(gridfile_hycom,'plat:',[0]))

#%% Reading depths
depth_hycom = readdepth(depthfile_hycom,'depth',pntidx=None)

#%% Extracting the longitudinal and latitudinal array size
afile = depthfile_hycom + '.a'
lines=[line.rstrip() for line in open(afile[:-2]+'.b')]
idm = int([line.split() for line in lines if 'i/jdm' in line][0][2])
jdm = int([line.split() for line in lines if 'i/jdm' in line][0][3].split(';')[0])

#%% All domain
fig, ax = plt.subplots(figsize=(12, 6))
plt.contourf(np.arange(1,idm+1),np.arange(1,jdm+1),depth_hycom,cmap=plt.cm.Spectral_r)
plt.colorbar()
plt.plot(np.arange(1,idm+1),np.tile(1,idm),'.k')
plt.plot(np.arange(1,idm+1),np.tile(jdm,idm),'.k')
plt.plot(np.tile(1,jdm),np.arange(1,jdm+1),'.k')
plt.plot(np.tile(idm,jdm),np.arange(1,jdm+1),'.k')
plt.title(depthfile_hycom)

fig, ax = plt.subplots(figsize=(12, 6))
plt.contourf(lon_hycom,lat_hycom,depth_hycom,cmap=plt.cm.Spectral_r)
plt.colorbar()
plt.title(gridfile_hycom)


