
#%% User input
ncfiles = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv2p0p1a_2024rt/2024092512/09l/09l.2024092512.hfsa.mom6.f000.nc','/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv2p0p1a_2024rt_BBL_EFFIC/2024092512/09l/09l.2024092512.hfsa.mom6.f000.nc']

#ncfiles = ['/work/noaa/hwrf/noscrub/maristiz/hafsv2_hnod_hnod/00l.2024041506.hfsa.mom6.f126.nc','/work/noaa/hwrf/noscrub/maristiz/hafsv2_SST_2DVAR_hnod/SST_SSS_adjust_prof_HYCOM_mld/00l.2024041506.hfsa.mom6.f126.nc']

#ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_uv_ic_file_5days/ocean_hourly_2d.nc'

#ncfile = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_4_degree_MOM6_TS_UV_SSH_SAVE_INITIAL_noshift_newflooding/ocean_hourly_2d.nc'

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
#####################################################################
ncfile = ncfiles[0]
ncdata = xr.open_dataset(ncfile,decode_times=False)
yh = np.asarray(ncdata['yh'][:])
yq = np.asarray(ncdata['yq'][:])
xh = np.asarray(ncdata['xh'][:])
xq = np.asarray(ncdata['xq'][:])

SST = np.empty((len(ncfiles),len(yh),len(xh)))
SST[:] = np.nan
SSS = np.empty((len(ncfiles),len(yh),len(xh)))
SSS[:] = np.nan
SSH = np.empty((len(ncfiles),len(yh),len(xh)))
SSH[:] = np.nan
SSU = np.empty((len(ncfiles),len(yh),len(xq)))
SSU[:] = np.nan
SSV = np.empty((len(ncfiles),len(yq),len(xh)))
SSV[:] = np.nan
LW = np.empty((len(ncfiles),len(yh),len(xh)))
LW[:] = np.nan
SW = np.empty((len(ncfiles),len(yh),len(xh)))
SW[:] = np.nan
SEN = np.empty((len(ncfiles),len(yh),len(xh)))
SEN[:] = np.nan
LAT = np.empty((len(ncfiles),len(yh),len(xh)))
LAT[:] = np.nan

for exp,ncfile in enumerate(ncfiles):
    #%% Read nc file
    ncdata = xr.open_dataset(ncfile,decode_times=False)
    yh = np.asarray(ncdata['yh'][:])
    yq = np.asarray(ncdata['yq'][:])
    xh = np.asarray(ncdata['xh'][:])
    xq = np.asarray(ncdata['xq'][:])

    tt = ncdata['time']
    ref_year = int(tt.units.split(' ')[2].split('-')[0])
    ref_month = int(tt.units.split(' ')[2].split('-')[1])
    ref_day = int(tt.units.split(' ')[2].split('-')[2])
    ref_hour = int(tt.units.split(' ')[3].split('-')[0].split(':')[0])
    ftime = np.asarray(tt)[0]
    time = datetime(ref_year,ref_month,ref_day,ref_hour) + timedelta(hours=ftime)
    
    sst = np.asarray(ncdata['SST'][0,:,:])
    sss = np.asarray(ncdata['SSS'][0,:,:])
    ssh = np.asarray(ncdata['SSH'][0,:,:])
    ssu = np.asarray(ncdata['SSU'][0,:,:])
    ssv = np.asarray(ncdata['SSV'][0,:,:])
    sw = np.asarray(ncdata['SW'][0,:,:])
    lw = np.asarray(ncdata['LW'][0,:,:])
    sen = np.asarray(ncdata['sensible'][0,:,:])
    lat = np.asarray(ncdata['latent'][0,:,:])
    
    SST[exp,:,:] = sst
    SSS[exp,:,:] = sss
    SSH[exp,:,:] = ssh
    SSU[exp,:,:] = ssu
    SSV[exp,:,:] = ssv
    SW[exp,:,:] = sw
    LW[exp,:,:] = lw
    SEN[exp,:,:] = sen
    LAT[exp,:,:] = lat
    ######################################################################################
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
    ax.set_title('SST '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sst_exp_' + str(exp))
    #plt.close()
    
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
    ax.set_title('SSS '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sss_exp_' + str(exp))
    #plt.close()
    
    # SSH
    levels = np.arange(-100,110,10)
    ticks = np.arange(-100,110,20)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,ssh*100,cmap='nipy_spectral',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('cm',fontsize=12)
    #cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('SSH '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_ssh_exp_' + str(exp))
    #plt.close()
    
    # Surface u vel
    levels = np.arange(-1,1.1,0.1)
    ticks = np.arange(-1,1.1,0.2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xq,yh,ssu,cmap='YlOrBr',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('m/s',fontsize=12)
    #cs = ax.contour(xh,yh,speed[nt,:,:],[0.3],colors='grey',alpha=0.2)
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Surface u vel. '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_ssu_exp_' + str(exp))
    #plt.close()
    
    # Surface v vel
    levels = np.arange(-1,1.1,0.1)
    ticks = np.arange(-1,1.1,0.2)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yq,ssv,cmap='YlOrBr',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('m/s',fontsize=12)
    #cs = ax.contour(xh,yh,speed[nt,:,:],[0.3],colors='grey',alpha=0.2)
    #ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Surface v vel. '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_ssv_exp_' + str(exp))
    #plt.close()

    # SW
    levels = np.arange(0,401,50)
    ticks = np.arange(0,401,50)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,sw,cmap='YlOrRd',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$W/m^2$',fontsize=12)
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Short Wave '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sw_exp_' + str(exp))
    #plt.close()

    # LW
    levels = np.arange(-100,101,20)
    ticks = np.arange(-100,101,20)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,lw,cmap='PuOr_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$W/m^2$',fontsize=12)
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Longwave '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_lw_exp_' + str(exp))
    #plt.close()
    
    # latent
    levels = np.arange(-300,301,50)
    ticks = np.arange(-300,301,50)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,lat,cmap='PuOr_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$W/m^2$',fontsize=12)
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Latent '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_lat_exp_' + str(exp))
    #plt.close()
    
    # Sensible
    levels = np.arange(-100,101,10)
    ticks = np.arange(-100,101,20)
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    cf = ax.contourf(xh,yh,sen,cmap='PuOr_r',levels=levels,extend='both')
    cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
    cb.set_label('$W/m^2$',fontsize=12)
    ax.clabel(cs,inline=True)
    ax.coastlines(resolution='50m')
    ax.set_title('Sensible '+ str(time),fontsize=16)
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig('mom6_sen_exp_' + str(exp))
    #plt.close()
    
# Change figures
# SST
SST_diff = SST[0,:,:] - SST[1,:,:]
levels = np.arange(-2,2.1,0.2)
ticks = np.arange(-2,2.1,0.4)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,SST_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$^oC$',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SST diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_sst_diff')
#plt.close()

# SSS
SSS_diff = SSS[0,:,:] - SSS[1,:,:]
levels = np.arange(-1,1.1,0.1)
ticks = np.arange(-1,1.1,0.2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,SSS_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label(' ',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SSS diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_sss_diff')
#plt.close()

# SSH
SSH_diff = SSH[0,:,:] - SSH[1,:,:]
levels = np.arange(-10,11,1)
ticks = np.arange(-10,11,2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,SSH_diff*100,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('cm',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SSH diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_ssh_diff')
#plt.close()

# SSU
SSU_diff = SSU[0,:,:] - SSU[1,:,:]
levels = np.arange(-1,1,0.1)
ticks = np.arange(-1,1,0.2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xq,yh,SSU_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('m/s',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SSU diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_ssu_diff')
#plt.close()

# SSV
SSV_diff = SSV[0,:,:] - SSV[1,:,:]
levels = np.arange(-1,1,0.1)
ticks = np.arange(-1,1,0.2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yq,SSV_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('m/s',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SSV diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_ssv_diff')
#plt.close()

# SW
SW_diff = SW[0,:,:] - SW[1,:,:]
levels = np.arange(-2,2.1,0.2)
ticks = np.arange(-2,2.1,0.4)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,SW_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$W/cm^2$',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SW diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_sw_diff')
#plt.close()

# LW
LW_diff = LW[0,:,:] - LW[1,:,:]
levels = np.arange(-2,2.1,0.2)
ticks = np.arange(-2,2.1,0.4)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,LW_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$W/cm^2$',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('LW diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_lw_diff')
#plt.close()

# Latent
LAT_diff = LAT[0,:,:] - LAT[1,:,:]
levels = np.arange(-400,400.1,40)
ticks = np.arange(-400,400.1,80)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,LAT_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$W/cm^2$',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('LAT diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_lat_diff')
#plt.close()

# Sensible
SEN_diff = SEN[0,:,:] - SEN[1,:,:]
levels = np.arange(-200,200.1,20)
ticks = np.arange(-200,200.1,40)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,SEN_diff,cmap='bwr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=ticks)
cb.set_label('$W/cm^2$',fontsize=12)
#cs = ax.contour(xh,yh,ssh[nt,:,:]*100,[0],colors='k')
#ax.clabel(cs,inline=True)
ax.coastlines(resolution='50m')
ax.set_title('SEN diff '+ str(time),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
plt.savefig('mom6_sen_diff')
#plt.close()
