
#%% User input

ncfiles_dirs = ['/work/noaa/hwrf/scrub/maristiz/hafsv1_merge_mom6_20230703_hfsa_mom6_24h/2020082512/13L/ocn_prep/mom6_forcings/','/work/noaa/hwrf/scrub/maristiz/hafsv1_merge_mom6_20230703_hfsa_mom6_24h/com/2020082512/13L/']

cartopyDataDir = '/work/noaa/hwrf/local/share/cartopy/'

#######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates 
import xarray as xr
import glob
import os
import netCDF4
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter

import pyproj
import cartopy
import cartopy.crs as ccrs

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#plt.switch_backend('agg')
######################################################################################
#%% Read nc files

files0 = sorted(glob.glob(os.path.join(ncfiles_dirs[0],'*DSWRF.nc')))
ncdata = xr.open_dataset(files0[0],decode_times=False)
yh0 = np.asarray(ncdata['latitude'][:])
xh0 = np.asarray(ncdata['longitude'][:])
tt = np.asarray(ncdata['time'][:])
time0 = np.asarray([datetime(1970,1,1) + timedelta(seconds=t) for t in tt])
dswrf0 = np.asarray(ncdata['DSWRF_surface'][:])

files1 = sorted(glob.glob(os.path.join(ncfiles_dirs[1],'*DSWRF.nc')))

dswrf1 = np.empty((len(files1),1361,1681))
dswrf1[:] = np.nan
time1 = []
for n,file in enumerate(files1):
    ncdata = xr.open_dataset(file,decode_times=False)
    yh1 = np.asarray(ncdata['latitude'][:])
    xh1 = np.asarray(ncdata['longitude'][:])
    tt = np.asarray(ncdata['time'][:])
    time1.append(datetime(1970,1,1) + timedelta(seconds=tt[0]))
    dswrf1[n,:,:] = np.asarray(ncdata['DSWRF_surface'][:])

# time series

#cf = ax.contourf(xh00[okx00],yh0,dswrf0[nt,:,okx00].T,cmap='Spectral_r',levels=levels,extend='max')
xh00 = np.asarray([x-360 if x>=180 else x for x in xh0])
okx00 = np.argsort(xh00)
xh0_sort = xh00[okx00]
dswrf0_sort = dswrf0[:,:,okx00]

date_form = DateFormatter("%d-%H")
okyh0 = np.where(yh0 >= 10)[0][0] 
okxh0 = np.where(xh0_sort >= -128)[0][0] 
okyh1 = np.where(yh1 >= 10)[0][0] 
okxh1 = np.where(xh1 >= -128)[0][0] 
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(time0,dswrf0_sort[:,okyh0,okxh0],'.-',label='GFS forcing')
ax.plot(time1,dswrf1[:,okyh1,okxh1],'.-',label='FV3 HAFS-MOM6')
ax.xaxis.set_major_formatter(date_form)
plt.legend()

okyh0 = np.where(yh0 >= 10)[0][0] 
okxh0 = np.where(xh0_sort >= -55)[0][0] 
okyh1 = np.where(yh1 >= 10)[0][0] 
okxh1 = np.where(xh1 >= -55)[0][0] 
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(time0,dswrf0_sort[:,okyh0,okxh0],'.-',label='GFS forcing')
ax.plot(time1,dswrf1[:,okyh1,okxh1],'.-',label='FV3 HAFS-MOM6')
ax.xaxis.set_major_formatter(date_form)
plt.legend()

okyh0 = np.where(yh0 >= 23)[0][0] 
okxh0 = np.where(xh0_sort >= -86)[0][0] 
okyh1 = np.where(yh1 >= 23)[0][0] 
okxh1 = np.where(xh1 >= -86)[0][0] 
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(time0,dswrf0_sort[:,okyh0,okxh0],'.-',label='GFS forcing')
ax.plot(time1,dswrf1[:,okyh1,okxh1],'.-',label='FV3 HAFS-MOM6')
ax.xaxis.set_major_formatter(date_form)
plt.legend()

