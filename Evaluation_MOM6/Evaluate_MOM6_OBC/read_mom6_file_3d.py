
#%% User input
ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_quater_degree_interpolating/ocean_hourly_3d.nc'

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
xq = np.asarray(ncdata['xq'][:])
yh = np.asarray(ncdata['yh'][:])
z_l = np.asarray(ncdata['z_l'][:])
z_i = np.asarray(ncdata['z_i'][:])
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
xh = np.asarray(ncdata['xh'][:])
yq = np.asarray(ncdata['yq'][:])
u = np.asarray(ncdata['u'][:])
v = np.asarray(ncdata['v'][:])
temp = np.asarray(ncdata['temp'][:])
salt = np.asarray(ncdata['salt'][:])

######################################################################################
#for nt,tt in enumerate(time[0:3]):
for nt,tt in enumerate(time):
    print(nt)

    # Temperature
    levels = np.arange(26,30.2,0.2)
    #ticks = np.arange(15,34,2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,temp[nt,0,:,:],cmap='Spectral_r',levels=levels,extend='both')
    #cf = ax.contourf(xh,yh,temp[nt,0,:,:],cmap='prism',levels=levels,extend='both')
    plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True) #, ticks=ticks)
    cs = ax.contour(xh,yh,temp[nt,0,:,:],[26],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SST '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sst_f' + str(nt+1))
    plt.close()

    # Salinity
    levels = np.arange(32,38,0.5)
    ticks = np.arange(32,38,0.5)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,salt[nt,0,:,:],cmap='GnBu_r',levels=levels,extend='both')
    plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cs = ax.contour(xh,yh,salt[nt,0,:,:],[35,37],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSS '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sss_f' + str(nt+1))
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

