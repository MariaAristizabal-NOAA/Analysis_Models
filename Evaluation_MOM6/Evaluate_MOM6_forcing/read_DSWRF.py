
#%% User input

ncfiles = ['/work/noaa/hwrf/scrub/maristiz/hafsv1_merge_mom6_20230703_hfsa_mom6/2020082512/13L/ocn_prep/mom6_forcings/gfs_global_2020082512_DSWRF.nc','/work/noaa/hwrf/scrub/maristiz/hafsv1_merge_mom6_20230703_hfsa_mom6_24h/com/2020082512/13L/13l.2020082512.hfsa.mom6.f003.nc','/work/noaa/hwrf/scrub/maristiz/hafsv1_merge_mom6_20230703_hfsa_mom6/com/2020082512/13L/13l.2020082512.hfsa.parent.atm.f003.DSWRF.nc','/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/hafs_mom6_nhc_test_obc/out1_yes_UV_ic_yes_OBCs_1day_run_corr_west_OK/ocean_hourly_2d.nc','/work/noaa/hwrf/noscrub/maristiz/HFSA_oper/2020082512/13l.2020082512.hfsa.hycom.2d.f003.nc','/work/noaa/hwrf/noscrub/maristiz/HFSA_oper/2020082512/13l.2020082512.hfsa.parent.atm.DSWRF.f003.nc']

cartopyDataDir = '/work/noaa/hwrf/local/share/cartopy/'

#######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates 
import xarray as xr
import glob
import netCDF4
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

ncdata = xr.open_dataset(ncfiles[0],decode_times=False)
yh0 = np.asarray(ncdata['latitude'][:])
xh0 = np.asarray(ncdata['longitude'][:])
tt = np.asarray(ncdata['time'][:])
time0 = np.asarray([datetime(1970,1,1) + timedelta(seconds=t) for t in tt])
dswrf0 = np.asarray(ncdata['DSWRF_surface'][:])

ncdata = xr.open_dataset(ncfiles[1],decode_times=False)
yh1 = np.asarray(ncdata['yh'][:])
xh1 = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time1 = np.asarray([datetime(2020,8,25,12) + timedelta(hours=t) for t in tt])
dswrf1 = np.asarray(ncdata['SW'][:])

ncdata = xr.open_dataset(ncfiles[2],decode_times=False)
yh2 = np.asarray(ncdata['latitude'][:])
xh2 = np.asarray(ncdata['longitude'][:])
tt = np.asarray(ncdata['time'][:])
time2 = np.asarray([datetime(1970,1,1) + timedelta(seconds=t) for t in tt])
dswrf2 = np.asarray(ncdata['DSWRF_surface'][:])

ncdata = xr.open_dataset(ncfiles[3],decode_times=False)
yh3 = np.asarray(ncdata['yh'][:])
xh3 = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time3 = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
dswrf3 = np.asarray(ncdata['SW'][:])

ncdata = xr.open_dataset(ncfiles[4],decode_times=False)
yh4 = np.asarray(ncdata['Latitude'][:])
xh4 = np.asarray(ncdata['Longitude'][:])
tt = np.asarray(ncdata['MT'][:])
time4 = np.asarray([datetime(1900,12,31) + timedelta(days=t) for t in tt])
dswrf4 = np.asarray(ncdata['surface_heat_flux'][:])

ncdata = xr.open_dataset(ncfiles[5],decode_times=False)
yh5 = np.asarray(ncdata['latitude'][:])
xh5 = np.asarray(ncdata['longitude'][:])
tt = np.asarray(ncdata['time'][:])
time5 = np.asarray([datetime(1970,1,1) + timedelta(seconds=t) for t in tt])
dswrf5 = np.asarray(ncdata['DSWRF_surface'][:])

######################################################################################

xh00 = np.asarray([x-360 if x>=180 else x for x in xh0])
okx00 = np.argsort(xh00)
# Short wave
nt = 1
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
#cf = ax.contourf(xh0,yh0,dswrf0[nt,:,:],cmap='Spectral_r',levels=levels,extend='max')
cf = ax.contourf(xh00[okx00],yh0,dswrf0[nt,:,okx00].T,cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW GFS Input '+ str(time0[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)

# Short wave
nt = 0
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh1,yh1,dswrf1[nt,:,:],cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW HAFS-MOM6 '+ str(time1[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)
            
# Short wave
nt = 0
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh2,yh2,dswrf2[nt,:,:],cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW FV3 from HAFS-MOM6 '+ str(time2[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)

# Short wave
nt = 2
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh3,yh3,dswrf3[nt,:,:],cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW ocean_only '+ str(time3[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)

# Short wave
nt = 0
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh4,yh4,dswrf4[nt,:,:],cmap='Spectral_r') #,levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True) #, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('Surface heat flux HYCOM '+ str(time4[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)

# Short wave
nt = 0
levels = np.arange(0,1001,100)
ticks = np.arange(0,1001,100)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh5,yh5,dswrf5[nt,:,:],cmap='Spectral_r',levels=levels,extend='max')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SW FV3 form HAFS-HYCOM '+ str(time5[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim(-178,15)
ax.set_ylim(-23,47)

okyh0 = np.where(yh0 >= 23)[0][0] 
okyh1 = np.where(yh1 >= 23)[0][0] 
okyh2 = np.where(yh2 >= 23)[0][0] 
okyh3 = np.where(yh3 >= 23)[0][0] 

plt.figure()
plt.plot(xh00[okx00],dswrf0[1,okyh0,okx00],label='GFS')
plt.plot(xh3,dswrf3[2,okyh3,:],label='MOM6 ocean_only')
plt.plot(xh2,dswrf2[0,okyh2,:],label='atm HAFS-MOM6')
plt.plot(xh1,dswrf1[0,okyh1,:],label='ocean HAFS-MOM6')
plt.legend()
plt.xlim(-178,15)

okyh2 = np.where(yh2 >= 23)[0][0] 
okyh5 = np.where(yh2 >= 23)[0][0] 

cf = ax.contourf(xh00[okx00],yh0,dswrf0[nt,:,okx00].T,cmap='Spectral_r',levels=levels,extend='max')
plt.figure()
plt.plot(xh2,dswrf2[0,okyh2,:],label='atm HAFS-MOM6')
plt.plot(xh5,dswrf5[0,okyh5,:],label='atm HAFS-HYCOM')
plt.legend()

