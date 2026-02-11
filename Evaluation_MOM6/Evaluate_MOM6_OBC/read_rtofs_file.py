
#%% User input
ncfile_ts = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512/INPUT/rtofs_ts_ic.nc'
ncfile_ssh = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512/INPUT/rtofs_ssh_ic.nc'
ncfile_uv = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512/INPUT/rtofs_uv_ic.nc'

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

######################################################################################
#%% Read nc files
ncdata = xr.open_dataset(ncfile_ts,decode_times=False)
yh = np.asarray(ncdata['lath'][:])
xh = np.asarray(ncdata['lonh'][:])
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(1900,12,31,00) + timedelta(days=t) for t in tt])
sst = np.asarray(ncdata['Temp'][:])[0,0,:,:]
sss = np.asarray(ncdata['Salt'][:])[0,0,:,:]

ncdata = xr.open_dataset(ncfile_ssh,decode_times=False)
ssh = np.asarray(ncdata['ave_ssh'][:])[0,:,:]

ncdata = xr.open_dataset(ncfile_uv,decode_times=False)
u = np.asarray(ncdata['u'][:])
v = np.asarray(ncdata['v'][:])
speed = np.sqrt(u[0,0,:,:]**2 + v[0,0,:,:]**2)

######################################################################################
for nt,tt in enumerate(time[0:1]):

    # Temperature surface
    levels = np.arange(15,34,1)
    ticks = np.arange(15,34,2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,sst,cmap='Spectral_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$^oC$',fontsize=12)
    cs = ax.contour(xh,yh,sst,[26],colors='k')
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
    cf = ax.contourf(xh,yh,sss,cmap='GnBu_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label(' ',fontsize=12)
    cs = ax.contour(xh,yh,sss,[35,37],colors='k')
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
    cf = ax.contourf(xh,yh,ssh*100,cmap='nipy_spectral',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('cm',fontsize=12)
    cs = ax.contour(xh,yh,ssh*100,[0],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSH '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False

    # Surface speed
    levels = np.arange(0,2.1,0.1)
    ticks = np.arange(0,2.1,0.2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,speed,cmap='YlOrBr',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('m/s',fontsize=12)
    cs = ax.contour(xh,yh,speed,[35,37],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Surface speed '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False


'''
fig,ax = plt.subplots()
ax.axis('equal')
cf = ax.contourf(lon,lat,field[nt,:,:]) # ,**kw)
plt.colorbar(cf)
ax.set_title('Forcing file ' + str(time[nt]))
'''

'''
dates = matplotlib.dates.date2num(time)
labels = np.asarray([t.hour for t in time]) 
fig,ax = plt.subplots(figsize=(10,5))
fhour_even = [label%2==0 for label in labels] 
fhour_odd = [label%2!=0 for label in labels] 
ax.plot(dates,field_target,'*-k')
ax.plot(dates[fhour_even],field_target[fhour_even],'*g',label='ave from previous 6 hours')
ax.plot(dates[fhour_odd],field_target[fhour_odd],'*b',label='ave from previous 3 hours')
ax.set_xticks(dates[::2])
ax.set_xticklabels(labels[::2])
ax.set_xlabel('UTC start date '+str(time[0]),fontsize=14)
ax.grid()
ax.legend()
'''

