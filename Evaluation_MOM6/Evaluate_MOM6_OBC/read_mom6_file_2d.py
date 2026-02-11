
#%% User input
#ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_uv_ic_file_5days/ocean_hourly_2d.nc'

ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_4_degree_MOM6_TS_UV_SSH_SAVE_INITIAL_noshift_newflooding/ocean_hourly_2d.nc'

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

plt.switch_backend('agg')
######################################################################################
#%% Read nc file
ncdata = xr.open_dataset(ncfile,decode_times=False)
yh = np.asarray(ncdata['yh'][:])
xh = np.asarray(ncdata['xh'][:])
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
ssh = np.asarray(ncdata['ssh'][:])
speed = np.asarray(ncdata['speed'][:])
sst = np.asarray(ncdata['sst'][:])
sss = np.asarray(ncdata['sss'][:])

######################################################################################
#for nt,tt in enumerate(time[0:1]):
for nt,tt in enumerate(time):
    print(nt)
    '''
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
    '''

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

