
#%% User input
ncfile = '/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/hafs_mom6_nhc_test_obc/ocn_2020_08_25_15.nc'

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
#%% Read nc file
ncdata = xr.open_dataset(ncfile,decode_times=False)
xq = np.asarray(ncdata['xq'][:])
yh = np.asarray(ncdata['yh'][:])
z_l = np.asarray(ncdata['z_l'][:])
z_i = np.asarray(ncdata['z_i'][:])
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(2020,8,25,12) + timedelta(days=t) for t in tt])
xh = np.asarray(ncdata['xh'][:])
yq = np.asarray(ncdata['yq'][:])
rho2_l = np.asarray(ncdata['rho2_l'][:])
rho2_i = np.asarray(ncdata['rho2_i'][:])
#u = np.asarray(ncdata['u'][:])
#v = np.asarray(ncdata['v'][:])
uh_rho_z = np.asarray(ncdata['uh_rho_z'][:])
vh_rho_z = np.asarray(ncdata['vh_rho_z'][:])
uh_rho_rho2 = np.asarray(ncdata['uh_rho_rho2'][:])
vh_rho_rho2 = np.asarray(ncdata['vh_rho_rho2'][:])

######################################################################################
#for nt,tt in enumerate(time[0:3]):
for nt,tt in enumerate(time):
    print(nt)

    # Temperature
    levels = np.arange(26,30.2,0.2)
    #ticks = np.arange(15,34,2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xq,yh,uh_rho_z[0,0,:,:],cmap='Spectral_r') #,levels=levels,extend='both')
    #cf = ax.contourf(xh,yh,temp[nt,0,:,:],cmap='prism',levels=levels,extend='both')
    plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True) #, ticks=ticks)
    #cs = ax.contour(xh,yh,temp[nt,0,:,:],[26],colors='k')
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    #ax.set_title('SST '+ str(time[nt]),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    #plt.savefig('mom6_sst_f' + str(nt+1))
    #plt.close()

div_x = np.diff(uh_rho_z[0,0,:,:],axis=1) 
div_y = np.diff(vh_rho_z[0,0,:,:],axis=0) 
div_z = -div_x - div_y
dx = np.diff(xq)*110*10**3 # meters
dy = np.diff(yq)*110*10**3 # meters
dA = dy.reshape(964,1)*dx.reshape(1,2413) # m^2
w2 = div_z/dA * 3600*24 # m/day

div_x = np.diff(uh_rho_rho2[0,0,:,:],axis=1)
div_y = np.diff(vh_rho_rho2[0,0,:,:],axis=0)
div_z = -div_x - div_y
dx = np.diff(xq)*110*10**3 # meters
dy = np.diff(yq)*110*10**3 # meters
dA = dy.reshape(964,1)*dx.reshape(1,2413) # m^2
w2 = div_z/dA * 3600*24 # m/day


fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,div_x,cmap='Spectral_r') #,levels=levels,extend='both')
plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True) #, ticks=ticks)
ax.coastlines(resolution='50m')


fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,w2,cmap='Spectral_r',levels=np.arange(-60,61,10),extend='both')
plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True) #, ticks=ticks)
ax.coastlines(resolution='50m')


