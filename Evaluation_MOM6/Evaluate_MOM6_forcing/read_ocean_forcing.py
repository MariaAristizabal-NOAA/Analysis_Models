
#%% User input
ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-Aug9_2024/ocean_only/RTOFS_2DVAR_SST_baseline_tests/INPUT/ocean_forcings.nc']
#ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-Aug9_2024/ocean_only/RTOFS_2DVAR_SST_baseline_tests/INPUT/ocean_forcings.nc_orig']
#ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-Aug9_2024/ocean_only/hafs_mom6_hat10_ctrl/INPUT/ocean_forcings.nc']

cartopyDataDir = '/work/noaa/hwrf/local/share/cartopy/'

#######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates 
import xarray as xr
import glob
import netCDF4 as nc
from datetime import datetime, timedelta

import pyproj
import cartopy
import cartopy.crs as ccrs

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#plt.switch_backend('agg')
######################################################################################
#%% Read nc files

#ncdata = xr.open_dataset(ncfiles[0],decode_times=False)
ncdata = nc.Dataset(ncfiles[0],decode_times=False)
yh = np.asarray(ncdata['latitude'][:])
xhh = np.asarray(ncdata['longitude'][:])
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(1970,1,1) + timedelta(seconds=t) for t in tt])
netsw = np.asarray(ncdata['NETSW_surface'][:])
netlw = np.asarray(ncdata['NETLW_surface'][:])

nc_file = nc.Dataset(ncfiles[0],'a')
#netSW = np.copy(netsw)
#netSW[:,:,:] = 100
#nc_file['NETSW_surface'][:,:,:] = netSW
TT = np.copy(tt) 
#tt_ok = np.array([1.7052984e+09, 1.7053092e+09, 1.7053200e+09, 1.7053308e+09,1.7053416e+09, 1.7053524e+09])
#TT = tt_ok[0:5]
tt_ok = np.array([1598356800, 1598367600, 1598378400, 1598389200, 1598400000, 1598410800])
TT = tt_ok
nc_file['time'][:] = TT
nc_file.close()

######################################################################################

xh = np.asarray([x-360 if x>=180 else x for x in xhh])
okx = np.argsort(xh)
# Short wave
nt = 1
levels = np.arange(-0.1,1001,100)
ticks = np.arange(-0.1,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh[okx],yh,netsw[nt,:,okx].T,cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW GFS Input '+ str(time[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
#ax.set_xlim(-178,15)
#ax.set_ylim(-23,47)

# Long wave
nt = 1
levels = np.arange(-100,101,10)
ticks = np.arange(-100,101,10)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh[okx],yh,netlw[nt,:,okx].T,cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('LW GFS Input '+ str(time[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
#ax.set_xlim(-178,15)
#ax.set_ylim(-23,47)


'''
okyh0 = np.where(yh0 >= 23)[0][0] 

plt.figure()
plt.plot(xh00[okx00],dswrf0[1,okyh0,okx00],label='GFS')
plt.plot(xh3,dswrf3[2,okyh3,:],label='MOM6 ocean_only')
plt.plot(xh2,dswrf2[0,okyh2,:],label='atm HAFS-MOM6')
plt.plot(xh1,dswrf1[0,okyh1,:],label='ocean HAFS-MOM6')
plt.legend()
plt.xlim(-178,15)
'''
