#%% User input

# RTOFS grid file name
hycom_grid = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_IC/2020082512_JTNC/regional.grid'
hycom_depth = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_IC/2020082512_JTNC/regional.depth'
#afile = 'archv_rtofs'
#afile = 'archv_inc_mom6'
afile = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_IC/2020082512_JTNC/archv_in'
#afile = 'diff_temp_mom6_minus_rtofs'
#afile = 'diff_saln_mom6_minus_rtofs'

var_name = 'temp'
klayer = '1'
#klayer = '30'

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
lines = [line.rstrip() for line in open(afile+'.b')]
ijdm = idm*jdm
npad = 4096-(ijdm%4096)
fld = ma.array([],fill_value=1.2676506002282294e+30)
#fld2 = ma.array([],fill_value=1e30)

inFile = afile + '.a'

for n,line in enumerate(lines):
    if line.split()[0] == 'field':
        nheading = n + 1
        print(line.split()[0])

for n,line in enumerate(lines):
    if line.split()[0] == var_name and line.split()[4] == klayer:
        nvar = n - nheading
        print(nvar)
        print(n)

fid = open(inFile,'rb')
#fid.seek((nvar-1)*4*(npad+ijdm),0)
fid.seek((nvar)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld = ma.reshape(fld,(jdm,idm))

fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan 
print(np.max(fld2))
print(np.min(fld2))

for n,line in enumerate(lines):
    if line.split()[0] == 'field':
        nheading = n + 1
        print(line.split()[0])
    
for n,line in enumerate(lines):
    if line.split()[0] == var_name and line.split()[4] == klayer:
        nvar = n - nheading    
        print(nvar)

kw = dict(levels=np.arange(-6,6.1,0.5))
#kw = dict(levels=np.arange(-2,2.1,0.5))
fig,ax1 = plt.subplots(figsize = (7,5))
#plt.contourf(lon_hycom,lat_hycom,fld,cmap='Spectral_r') #,**kw)
plt.contourf(lon_hycom,lat_hycom,fld,cmap='Spectral_r') #,**kw)
plt.axis('scaled')
plt.colorbar()
plt.title(var_name + ' layer ' + klayer,fontsize=14)
plt.text(260,-10,'max val = ' + str(np.max(fld2)),fontsize=14)
plt.text(260,-15,'min val = ' + str(np.min(fld2)),fontsize=14)

kw = dict(levels=np.arange(0,36,1))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom,lat_hycom,fld2,cmap='Spectral_r',**kw)
plt.axis('scaled')
plt.colorbar()
plt.title(var_name + ' layer ' + klayer,fontsize=14)
plt.xlim([122,124])
plt.ylim([23,25])
#plt.text(260,-10,'max val = ' + str(np.max(fld2)),fontsize=14)
#plt.text(260,-15,'min val = ' + str(np.min(fld2)),fontsize=14)

depth_hycom2 = np.copy(depth_hycom)
mask = depth_hycom2 > 10**5
depth_hycom2[mask] = np.nan 
# plot depth
kw = dict(levels=np.arange(0,8000,100))
fig,ax1 = plt.subplots(figsize = (7,5))
plt.contourf(lon_hycom,lat_hycom,depth_hycom,cmap='Spectral_r',**kw)
plt.axis('scaled')
plt.colorbar()
plt.title('Depth from regional.depth.a/b files ' ,fontsize=14)
plt.xlim([122,124])
plt.ylim([23,25])


####################################################################

