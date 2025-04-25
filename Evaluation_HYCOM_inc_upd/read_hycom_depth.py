#%% User input

depthfile_hycom = '/work/noaa/hwrf/save/maristiz/hafs_develop_202112_mom6_incr_update/fix/fix_hycom/hafs_hycom_hat10.basin.regional.depth'

#%%
import matplotlib.pyplot as plt
import numpy as np
import sys

import sys
sys.path.append('/home/maristiz/Utils/HYCOM_utils/')
from utils4HYCOM import readBinz, readgrids, parse_z
from utils4HYCOM import readdepth

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

