
#%% User input
#ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_uv_ic_file_5days/ocean_hourly_2d.nc'

ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_12_degree_RTOFS_TS_UV_SSH/ocean_hourly_2d.nc','/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_12_degree_RTOFS_TS_UV_SSH_test_OBCs/OBC_persistent_1day/ocean_hourly_2d.nc','/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_12_degree_RTOFS_TS_UV_SSH_test_OBCs/OBC_persistent_6days/ocean_hourly_2d.nc']

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
yh0 = np.asarray(ncdata['yh'][:])
xh0 = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time0 = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
ssh0 = np.asarray(ncdata['ssh'][:])
speed0 = np.asarray(ncdata['speed'][:])
sst0 = np.asarray(ncdata['sst'][:])
sss0 = np.asarray(ncdata['sss'][:])

ncdata = xr.open_dataset(ncfiles[1],decode_times=False)
yh1 = np.asarray(ncdata['yh'][:])
xh1 = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time1 = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
ssh1 = np.asarray(ncdata['ssh'][:])
speed1 = np.asarray(ncdata['speed'][:])
sst1 = np.asarray(ncdata['sst'][:])
sss1 = np.asarray(ncdata['sss'][:])

ncdata = xr.open_dataset(ncfiles[2],decode_times=False)
yh2 = np.asarray(ncdata['yh'][:])
xh2 = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time2 = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
ssh2 = np.asarray(ncdata['ssh'][:])
speed2 = np.asarray(ncdata['speed'][:])
sst2 = np.asarray(ncdata['sst'][:])
sss2 = np.asarray(ncdata['sss'][:])

######################################################################################

plt.figure(figsize=(8,5))
plt.plot(time0,ssh0[:,0,1000],label='Daily OBCs')
plt.plot(time1,ssh1[:,0,1000],label='Persitent OBCs one time level')
plt.plot(time2,ssh2[:,0,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSH South Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,ssh0[:,-1,1000],label='Daily OBCs')
plt.plot(time1,ssh1[:,-1,1000],label='Persitent OBCs one time level')
plt.plot(time2,ssh2[:,-1,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSH North Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,ssh0[:,10,-1],label='Daily OBCs')
plt.plot(time1,ssh1[:,10,-1],label='Persitent OBCs one time level')
plt.plot(time2,ssh2[:,10,-1],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSH East Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sst0[:,0,1000],label='Daily OBCs')
plt.plot(time1,sst1[:,0,1000],label='Persitent OBCs one time level')
plt.plot(time2,sst2[:,0,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SST South Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sst0[:,-1,1000],label='Daily OBCs')
plt.plot(time1,sst1[:,-1,1000],label='Persitent OBCs one time level')
plt.plot(time2,sst2[:,-1,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SST North Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sst0[:,10,-1],label='Daily OBCs')
plt.plot(time1,sst1[:,10,-1],label='Persitent OBCs one time level')
plt.plot(time2,sst2[:,10,-1],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SST East Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sss0[:,0,1000],label='Daily OBCs')
plt.plot(time1,sss1[:,0,1000],label='Persitent OBCs one time level')
plt.plot(time2,sss2[:,0,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSS South Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sss0[:,-1,1000],label='Daily OBCs')
plt.plot(time1,sss1[:,-1,1000],label='Persitent OBCs one time level')
plt.plot(time2,sss2[:,-1,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSS North Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,sss0[:,10,-1],label='Daily OBCs')
plt.plot(time1,sss1[:,10,-1],label='Persitent OBCs one time level')
plt.plot(time2,sss2[:,10,-1],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSS East Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,speed0[:,0,1000],label='Daily OBCs')
plt.plot(time1,speed1[:,0,1000],label='Persitent OBCs one time level')
plt.plot(time2,speed2[:,0,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('Speed South Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,speed0[:,-1,1000],label='Daily OBCs')
plt.plot(time1,speed1[:,-1,1000],label='Persitent OBCs one time level')
plt.plot(time2,speed2[:,-1,1000],label='Persitent OBCs six time levels')
plt.legend()
plt.title('Speed North Boundary')

plt.figure(figsize=(8,5))
plt.plot(time0,speed0[:,10,-1],label='Daily OBCs')
plt.plot(time1,speed1[:,10,-1],label='Persitent OBCs one time level')
plt.plot(time2,speed2[:,10,-1],label='Persitent OBCs six time levels')
plt.legend()
plt.title('Speed East Boundary')

# Gulf of Mexico
plt.figure(figsize=(8,5))
plt.plot(time0,ssh0[:,300,300],label='Daily OBCs')
plt.plot(time1,ssh1[:,300,300],label='Persitent OBCs one time level')
plt.plot(time2,ssh2[:,300,300],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSH Gulf Mexico')

plt.figure(figsize=(8,5))
plt.plot(time0,sst0[:,300,300],label='Daily OBCs')
plt.plot(time1,sst1[:,300,300],label='Persitent OBCs one time level')
plt.plot(time2,sst2[:,300,300],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SST Gulf Mexico')

plt.figure(figsize=(8,5))
plt.plot(time0,sss0[:,300,300],label='Daily OBCs')
plt.plot(time1,sss1[:,300,300],label='Persitent OBCs one time level')
plt.plot(time2,sss2[:,300,300],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSS Gulf Mexico')

plt.figure(figsize=(8,5))
plt.plot(time0,speed0[:,300,300],label='Daily OBCs')
plt.plot(time1,speed1[:,300,300],label='Persitent OBCs one time level')
plt.plot(time2,speed2[:,300,300],label='Persitent OBCs six time levels')
plt.legend()
plt.title('Speed Gulf Mexico')

# North Atlantic
plt.figure(figsize=(8,5))
plt.plot(time0,ssh0[:,300,700],label='Daily OBCs')
plt.plot(time1,ssh1[:,300,700],label='Persitent OBCs one time level')
plt.plot(time2,ssh2[:,300,700],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSH North Atlantic')

plt.figure(figsize=(8,5))
plt.plot(time0,sst0[:,300,700],label='Daily OBCs')
plt.plot(time1,sst1[:,300,700],label='Persitent OBCs one time level')
plt.plot(time2,sst2[:,300,700],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SST North Atlantic')

plt.figure(figsize=(8,5))
plt.plot(time0,sss0[:,300,700],label='Daily OBCs')
plt.plot(time1,sss1[:,300,700],label='Persitent OBCs one time level')
plt.plot(time2,sss2[:,300,700],label='Persitent OBCs six time levels')
plt.legend()
plt.title('SSS North Atlantic')

plt.figure(figsize=(8,5))
plt.plot(time0,speed0[:,300,700],label='Daily OBCs')
plt.plot(time1,speed1[:,300,700],label='Persitent OBCs one time level')
plt.plot(time2,speed2[:,300,700],label='Persitent OBCs six time levels')
plt.legend()
plt.title('Speed North Atlantic')

'''
for nt,tt in enumerate(time):
    print(nt)
    # Temperature surface
    levels = np.arange(15,34,1)
    ticks = np.arange(15,34,2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,sst[nt,:,:],cmap='Spectral_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$^oC$',fontsize=12)
    cs = ax.contour(xh,yh,sst[nt,:,:],[26],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SST '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False

    # Salinity surface
    levels = np.arange(32,38,0.5)
    ticks = np.arange(32,38,0.5)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,sss[nt,:,:],cmap='GnBu_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label(' ',fontsize=12)
    cs = ax.contour(xh,yh,sss[nt,:,:],[35,37],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSS '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False

    # SSH
    levels = np.arange(-100,110,10)
    ticks = np.arange(-100,110,20)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,ssh[nt,:,:]*100,cmap='nipy_spectral',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('cm',fontsize=12)
    #cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSH '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_ssh_f' + str(nt+1))
    plt.close()

    # Surface speed
    levels = np.arange(0,2.1,0.1)
    ticks = np.arange(0,2.1,0.2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,speed[nt,:,:],cmap='YlOrBr',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('m/s',fontsize=12)
    #cs = ax.contour(xh,yh,speed[nt,:,:],[0.3],colors='grey',alpha=0.2)
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Surface speed '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_ssv_f' + str(nt+1))
    plt.close()

'''
