#%% User input

# RTOFS files
#hycom_depth = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/hafs_hycom_hat10.basin.regional.depth_no_pacif'
#hycom_grid = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/regional.grid'
hycom_depth = '/work/noaa/hwrf/noscrub/maristiz/RTOFS/fix/hafs_hycom_hat10.basin.regional.depth_99'
hycom_grid = '/work/noaa/hwrf/save/maristiz/hafs_develop_202112/fix/fix_hycom/hafs_hycom_hat10.basin.regional.grid'
#afile = '/work/noaa/hwrf/scrub/maristiz/direct_increm_update_OBC3DF_mom6_minus_rtofs/com/2020072812/00L/natl00l.2020072812.hafs_hycom_hat10.2020_210_18'
#afile = '/work/noaa/hwrf/noscrub/maristiz/MOM6_analysis_files/OBC3DF_conf/mom6.20200728/mom6_hat10.t00z.f12.archv'
#afile ='/work/noaa/hwrf/scrub/maristiz/direct_increm_update_OBC3DF_mom6_minus_rtofs/2020072812/00L/ocn_prep/hycominit1/archv2r_in'
#afile='/work/noaa/hwrf/scrub/maristiz/IC_from_OBC3DF_mom6_analysis/2020072812/00L/forecast/archv.2020_210_18'
afile='/work/noaa/hwrf/scrub/maristiz/hafs_restart_RTOFS/com/2020072812/00L/natl00l.2020072812.hafs_hycom_hat10.2020_210_18'


#afile = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/archv_mom6'
#afile = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/archv_rtofs'

#var_name = 'v-vel.'
klayer = '1'
#colormap='nipy_spectral'
colormap='OrRd'
#colormap='seismic'
#colormap='Spectral_r'
min_val = 0
max_val = 2
delta_val = 0.2  # delta in colorbar
#delta_contour = 2 # delta in contour plot
units = '$m/s$'
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]

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
#sys.path.append('/home/Maria.Aristizabal/Repos/NCEP_scripts/')
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
    if len(line) > 0:
        if line.split()[0] == 'field':
            nheading = n + 1
            print(line.split()[0])

for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == 'u-vel.' and line.split()[4] == klayer:
            nvar_u = n - nheading
            print(nvar_u)
            print(n)

for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == 'v-vel.' and line.split()[4] == klayer:
            nvar_v = n - nheading
            print(nvar_v)
            print(n)

fid = open(inFile,'rb')
#fid.seek((nvar-1)*4*(npad+ijdm),0)
fid.seek((nvar_u)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld_u = ma.reshape(fld,(jdm,idm))

fid = open(inFile,'rb')
#fid.seek((nvar-1)*4*(npad+ijdm),0)
fid.seek((nvar_v)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld_v = ma.reshape(fld,(jdm,idm))

fld = np.sqrt(fld_u**2 + fld_v**2)

fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan
oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
maxval = np.nanmax(fld2[:,oklon][oklat,:])
minval = np.nanmin(fld2[:,oklon][oklat,:])
meanval = np.nanmean(fld2[:,oklon][oklat,:])
print(maxval)
print(minval)
print(meanval)

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax1 = plt.subplots(figsize = (11,6))
ax1.set_facecolor("bisque")
#plt.contour(lon_hycom-360,lat_hycom,fld,np.arange(min_val,max_val,delta_contour),colors='grey')
#plt.contour(lon_hycom-360,lat_hycom,fld,[0],colors='k')
plt.contourf(lon_hycom-360,lat_hycom,fld,cmap=colormap,**kw)
plt.axis('scaled')
cbar = plt.colorbar(fraction=0.025, pad=0.04)
cbar.set_label(units,fontsize = 14)
plt.title('Currents' + ' layer ' + klayer,fontsize=18)
plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
plt.xlim(lon_lim)
plt.ylim(lat_lim)

####################################################################

