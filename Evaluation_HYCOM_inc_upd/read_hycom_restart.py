#%% User input

# RTOFS grid file name
hycom_grid = 'hafs_hycom_hat10.basin.regional.grid'
hycom_depth = 'hafs_hycom_hat10.basin.regional.depth'
#hycom_restart = '/work/noaa/hwrf/noscrub/maristiz/MOM6_analysis_files/hat10_mom6nc2archv_restart_2020_210_12_April19'
hycom_restart = '/work/noaa/hwrf/noscrub/maristiz/RTOFS/arch_restart/2020072812/restart_out'

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
sys.path.append('/home/maristiz/Utils/HYCOM_utils/')
from utils4HYCOM import readBinz, readgrids, parse_z
from utils4HYCOM import readdepth

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Reading hycom grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(hycom_grid,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid,'plat:',[0]))
depth_hycom = np.asarray(readdepth(hycom_depth,'depth'))

####################################################################
lines = [line.rstrip() for line in open(hycom_restart+'.b')]

field = 'temp'
nline_temp = [n for n,line in enumerate(open(hycom_restart+'.b')) if line.split(':')[0].split()[0] == field]

field = 'saln'
nline_saln = [n for n,line in enumerate(open(hycom_restart+'.b')) if line.split(':')[0].split()[0] == field]

ijdm = idm*jdm
npad = 4096-(ijdm%4096)
fld = ma.array([],fill_value=1.2676506002282294e+30)
#fld2 = ma.array([],fill_value=1e30)

inFile = hycom_restart + '.a'

nvar = 1
fid = open(inFile,'rb')
fid.seek((nvar-1)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld = ma.reshape(fld,(jdm,idm))

print(np.max(fld))
print(np.min(fld))
kw = dict(levels=np.arange(-6,6.1,0.5))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom,lat_hycom,fld,cmap='Spectral_r') #,**kw)
plt.axis('scaled')
plt.colorbar()
#plt.title('MOM6 temp - RTOFS temp')

####################################################################
'''
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
'''

