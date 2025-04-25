
#%% User input
ncfiles = ['/work/noaa/hwrf/scrub/maristiz/hafsv2_develop_July9_CTL_RTOFS_v2.5.test01_hfsa_dev/com/2024081200/05L/05l.2024081200.hfsa.mom6.f126.nc','/work/noaa/hwrf/scrub/maristiz/hafsv2_develop_July9_SST_SSS_2DVAR_RTOFS_v2.5.test01_hfsa_dev/com/2024081200/05L/05l.2024081200.hfsa.mom6.f126.nc']

#ncfiles = ['/work/noaa/hwrf/noscrub/maristiz/hafsv2_hnod_hnod/00l.2024041506.hfsa.mom6.f048.nc','/work/noaa/hwrf/noscrub/maristiz/hafsv2_SST_2DVAR_hnod/SST_SSS_adjust_prof_HYCOM_mld/00l.2024041506.hfsa.mom6.f048.nc']

#ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/RTOFS_2DVAR_SST_baseline/Exp_2024011506/ocn_2024_01_15_09.nc','/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/RTOFS_2DVAR_SST_parallel_test3/ocn_2024_01_15_09.nc']
#ncfiles = ['/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/RTOFS_2DVAR_SST_baseline/Exp_2024011506/ocn_2024_01_15_09.nc','/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/RTOFS_2DVAR_SST_parallel/Exp_2DVAR_SST_2024011506/ocn_2024_01_15_09.nc']

exp_labels = ['Original','Modified']

exp_colors = ['red','orange']

cartopyDataDir = '/work/noaa/hwrf/local/share/cartopy/'

#######################################################################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates 
import xarray as xr
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
# point of interest for vertical profiles
#x = 1093
#y = 34
#x = 1100
#y = 350
#x = 1490
#y = 830
#x = 200
#y = 350
#x = 667
#y = 352
#x = 1349
#y = 905
x = 1760
y = 891

temp_prof = np.empty((2,55))
temp_prof[:] = np.nan
salt_prof = np.empty((2,55))
salt_prof[:] = np.nan
sst = np.empty((2,964,2413))
sst[:] = np.nan
sss = np.empty((2,964,2413))
sss[:] = np.nan
for n,ncfile in enumerate(ncfiles):
    ncdata = xr.open_dataset(ncfile,decode_times=False)
    xh = np.asarray(ncdata['xh'][:])
    yh = np.asarray(ncdata['yh'][:])
    z_l = np.asarray(ncdata['z_l'][:])

    tt = ncdata['time']
    ref_year = int(tt.units.split(' ')[2].split('-')[0])
    ref_month = int(tt.units.split(' ')[2].split('-')[1])
    ref_day = int(tt.units.split(' ')[2].split('-')[2])
    ref_hour = int(tt.units.split(' ')[3].split('-')[0].split(':')[0])
    ftime = np.asarray(tt)[0]
    time = datetime(ref_year,ref_month,ref_day,ref_hour) + timedelta(hours=ftime)

    temp_prof[n,:] = np.asarray(ncdata['temp'][0,:,y,x])
    salt_prof[n,:] = np.asarray(ncdata['so'][0,:,y,x])
    SST = np.asarray(ncdata['temp'][0,0,:,:])
    SSS = np.asarray(ncdata['so'][0,0,:,:])
    sst[n,:,:] = SST
    sss[n,:,:] = SSS

    levels = np.arange(5,34.2,1)
    ticks = np.arange(15,34,2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,SST,cmap='Spectral_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=1.0, extendrect=True) #, ticks=ticks)
    cb.set_label('$^oC$',fontsize=12)
    cs = ax.contour(xh,yh,SST,[26],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.plot(xh[x],yh[y],'*',color='black',markersize=10)
    ax.set_title('SST '+ str(time)+' '+exp_labels[n],fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('SST_'+exp_labels[n]+'.png')

    levels = np.arange(32,38.1,0.2)
    ticks = np.arange(32,38.1,0.2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,SSS,cmap='ocean',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=1.0, extendrect=True) #, ticks=ticks)
    ax.plot(xh[x],yh[y],'*',color='black',markersize=10)
    cb.set_label(' ',fontsize=12)
    cs = ax.contour(xh,yh,SSS,[35,37],colors='k')
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSS '+ str(time)+' '+exp_labels[n],fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('SSS_'+exp_labels[n]+'.png')

sst_diff = sst[0,:,:]-sst[1,:,:]
np.where(sst_diff == np.nanmax(sst_diff))
levels = np.arange(-2,2.1,0.1)
ticks = np.arange(-2,2.1,0.1)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,sst_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=1.0, extendrect=True) #, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SST diff '+ str(time)+' '+exp_labels[n],fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.plot(xh[x],yh[y],'*',color='black',markersize=10)
plt.savefig('SST_diff.png')

sss_diff = sss[0,:,:]-sss[1,:,:]
np.where(sss_diff == np.nanmin(sss_diff))
levels = np.arange(-1,1.1,0.1)
ticks = np.arange(-1,1.11,0.1)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,sss_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=1.0, extendrect=True) #, ticks=ticks)
cb.set_label(' ',fontsize=12)
ax.coastlines(resolution='50m')
ax.set_title('SSS diff '+ str(time)+' '+exp_labels[n],fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.plot(xh[x],yh[y],'*',color='black',markersize=10)
plt.savefig('SSS_diff.png')

fig = plt.figure(figsize=(5,8))
plt.plot(temp_prof[0,:],-z_l,'o-',label=exp_labels[0])
plt.plot(temp_prof[1,:],-z_l,'o-',label=exp_labels[1])
plt.legend()
plt.ylim([-200,0])
plt.title('Temperature '+ str(time),fontsize=16)
plt.savefig('temp_prof.png')

fig = plt.figure(figsize=(5,8))
plt.plot(salt_prof[0,:],-z_l,'o-',label=exp_labels[0])
plt.plot(salt_prof[1,:],-z_l,'o-',label=exp_labels[1])
plt.legend()
plt.ylim([-200,0])
plt.title('Salinity '+ str(time),fontsize=16)
plt.savefig('salt_prof.png')

##############################################################

