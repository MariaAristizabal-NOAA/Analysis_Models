
#%% User input

'''
dirout = '/work/noaa/hwrf/scrub/yongzuo/HAFS_MOM6_coupling_2022092312/2022092312/00L/forecast/'
prefix = 'ocn_2022_09_'
dd = dirout + prefix
ncfiles_3d = [dd+'23_15.nc',dd+'23_18.nc',dd+'23_21.nc',dd+'24_00.nc',dd+'24_03.nc',dd+'24_06.nc',dd+'24_09.nc']
#ncfiles_3d = [dd+'25_15.nc',dd+'25_18.nc',dd+'25_21.nc',dd+'26_00.nc',dd+'26_03.nc',dd+'26_06.nc',dd+'26_09.nc',dd+'26_12.nc',dd+'25_18.nc',dd+'26_21.nc',dd+'27_00.nc',dd+'27_03.nc',dd+'27_06.nc',dd+'27_09.nc',dd+'27_12.nc',dd+'27_15.nc',dd+'27_18.nc',dd+'27_21.nc',dd+'28_00.nc',dd+'28_03.nc',dd+'28_06.nc',dd+'28_09.nc',dd+'28_12.nc',dd+'28_15.nc',dd+'28_18.nc',dd+'28_21.nc',dd+'29_00.nc',dd+'29_03.nc',dd+'29_06.nc',dd+'29_09.nc',dd+'29_12.nc',dd+'29_15.nc',dd+'29_18.nc',dd+'29_21.nc',dd+'30_00.nc',dd+'30_03.nc',dd+'30_06.nc',dd+'30_09.nc',dd+'30_12.nc',dd+'30_15.nc',dd+'30_18.nc']
'''

ncfile_2d = '/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_4_degree_MOM6_TS_UV_SSH/ocean_hourly_2d.nc'

ncfiles_3d = ['/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_12_degree_RTOFS_TS_UV_SSH/ocean_hourly_3d.nc','/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_IC_from_1_4_degree_MOM6_TS_UV_SSH_SAVE_INITIAL_noshift_newflooding/ocean_hourly_3d.nc']

#ncfiles_3d = ['/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_ST_ic_only_5days_good_forcing/ocean_hourly_3d.nc','/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_ST_UV_ic_no_SSH_5days_good_forcing/ocean_hourly_3d.nc','/work/noaa/hwrf/save/maristiz/MOM6_NOAA_EMC/z_data_table_HAT10_2020082512_uv_ic_file_5days_good_forcing/ocean_hourly_3d.nc']

#labels_exp = ['TS','TS_UV','TS_UV_SSH']
#colors_exp = ['indianred','green','blueviolet']

labels_exp = ['IC 1/12 RTOFS','IC 1/4 MOM6']
colors_exp = ['indianred','green']

FC_file = '/work/noaa/hwrf/save/maristiz/scripts_to_evaluate_MOM6_IC/FC_cable_transport_2020.dat'

cartopyDataDir = '/work/noaa/hwrf/local/share/cartopy/'

#######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates 
import xarray as xr
import glob
import netCDF4
from datetime import datetime, timedelta
from pyproj import Geod
import cartopy.crs as ccrs
import math

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#-----------------------------------------------
# coordinate rotation
#-----------------------------------------------
def rotate(u,v,phi):
    # phi is in radians
    u2 =  u*math.cos(phi) + v*math.sin(phi)
    v2 = -u*math.sin(phi) + v*math.cos(phi)
    return u2,v2

################################################################
#%% Read 3d nc file

#%% Read MOM6 grid
ncdata = xr.open_dataset(ncfiles_3d[0],decode_times=False)
xq = np.asarray(ncdata['xq'][:])
yh = np.asarray(ncdata['yh'][:])
z_l = np.asarray(ncdata['z_l'][:])
z_i = np.asarray(ncdata['z_i'][:])
xh = np.asarray(ncdata['xh'][:])
yq = np.asarray(ncdata['yq'][:])

################################################################
#%% Read 2d nc file
ncdata = xr.open_dataset(ncfile_2d,decode_times=False)
tt = np.asarray(ncdata['time'][:])
time = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
speed = np.asarray(ncdata['speed'][:])

################################################################
# Define cross section Florida Cable
x1 = -80.0533746
y1 = 26.7153425
#x1 = -80.14
#y1 = 26.74
x2 = -78.7833
y2 = 26.5167
# Slope
m = (y1-y2)/(x1-x2)
# Intercept
b = y1 - m*x1

X = np.arange(x1,x2,0.05)
Y = b + m*X

dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
g = Geod(ellps='WGS84')

################################################################
# Find u and v points along the Florida cable
okxh = np.round(np.interp(X,xh,np.arange(len(xh)))).astype(int)
okxq = np.round(np.interp(X,xq,np.arange(len(xq)))).astype(int)
okyh = np.round(np.interp(Y,yh,np.arange(len(yh)))).astype(int)
okyq = np.round(np.interp(Y,yq,np.arange(len(yq)))).astype(int)

# delta in the z direction
delta_zi = np.diff(z_i) # in meters

# Distance in the along cable direction
_,_,dist_cable=g.inv(xh[okxh[0:-1]],yq[okyq[0:-1]],xh[okxh[1:]],yq[okyq[1:]])

# Make distance and depth matrixes
delta_x,delta_z = np.meshgrid(dist_cable,delta_zi)

# Calculate area matrix
area_u = delta_z * delta_x # in meters square

################################################################
# Read 3d nc file

time_vec = []
#trans_total = []
time_total = np.empty((len(ncfiles_3d),len(tt)))
time_total[:] = np.nan
trans_total = np.empty((len(ncfiles_3d),len(tt)))
trans_total[:] = np.nan
for nn,ncfile_3d in enumerate(ncfiles_3d):
    print(ncfile_3d)
    ncdata = xr.open_dataset(ncfile_3d,decode_times=False)
    tt = np.asarray(ncdata['time'][:])
    time = np.asarray([datetime(2020,1,1,12) + timedelta(days=t) for t in tt])
    #u = np.asarray(ncdata['u'][:])
    #v = np.asarray(ncdata['v'][:])
    # Find u and v component of velocity along the cable
    ucable = np.empty((len(time),area_u.shape[0],len(okyh)))
    ucable[:] = np.nan
    vcable = np.empty((len(time),area_u.shape[0],len(okyh)))
    vcable[:] = np.nan
    for i in np.arange(len(okyh)):
        print(i)
        ucable[:,:,i] = ncdata['u'][:,:,okyh[i],okxq[i]].values
        vcable[:,:,i] = ncdata['v'][:,:,okyq[i],okxh[i]].values
        #ucable[:,:,i] = ncdata['u'][:,:,okyh[i],okxq[i]]
        #vcable[:,:,i] = ncdata['v'][:,:,okyq[i],okxh[i]]
    # Rotate velocities according with the cable's angle
    cable_angle = math.atan(m)
    ucable_rot,vcable_rot = rotate(ucable,vcable,cable_angle)
    for t,tt in enumerate(time):
        print(t)
        time_vec.append(tt)
        time_total[nn,t] = matplotlib.dates.date2num(tt)
        trans_matrix = vcable_rot[t,:,1:] * area_u # cubic meter per second
        #trans_total.append(np.nansum(trans_matrix) / 10**6) # in Sverdrups
        trans_total[nn,t] = (np.nansum(trans_matrix) / 10**6) # in Sverdrups

################################################################
# Read Florida cable trasnport
ff = open(FC_file,'r')
f = ff.readlines()

fc_time = []
fc_transp = []
for l in f:
    if l[0] != '%':
        year = int(l.split()[0])
        month = int(l.split()[1])
        day = int(l.split()[2])
        fc_time.append(datetime(year,month,day))
        fc_transp.append(float(l.split()[3]))

fc_time = np.asarray(fc_time)
fc_transp = np.asarray(fc_transp)

################################################################
# Surface speed
nt = 0
levels = np.arange(0,2.1,0.1)
ticks = np.arange(0,2.1,0.2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,speed[nt,:,:],cmap='YlOrBr',levels=levels,extend='both')
cb = plt.colorbar(cf,orientation='vertical', pad=0.02, aspect=30, shrink=1, extendrect=True, ticks=ticks)
cb.set_label('m/s',fontsize=12)
ax.plot(X,Y,'-g')
#ax.plot(xh[okxh],yh[okyh],'*r')
ax.coastlines(resolution='10m')
ax.set_title('Surface speed '+ str(time[nt]),fontsize=16)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim([-84,-74])
ax.set_ylim([23,30])

# u and v points along the Florida cable
levels = np.arange(0,2.1,0.1)
ticks = np.arange(0,2.1,0.2)
fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
cf = ax.contourf(xh,yh,speed[nt,:,:],cmap='YlOrBr',levels=levels,extend='both')
ax.plot(X,Y,'-g')
#ax.plot(xq[okxq],yh[okyh],'*r',label='U points')
ax.plot(xh[okxh],yq[okyq],'*b',label='V points')
#ax.plot(xq[okxq],yq[okyq],'*k',label='Q points')
ax.legend()
ax.coastlines(resolution='10m')
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,alpha=0.1)
gl.top_labels = False
gl.right_labels = False
ax.set_xlim([-80.5,-78.5])
ax.set_ylim([26,27])

# Transport
dates = matplotlib.dates.date2num(time_vec)
labels = np.asarray([str(t)[8:10] for t in time_vec]) 
#xticks = np.asarray([dates[n] if t.hour==0 or t.hour==12 else '' for n,t in enumerate(time)]) 
#labels = np.asarray([str(t)[6:13]+'h' if t.hour==0 or t.hour==12 else '' for t in time]) 
fig,ax = plt.subplots(figsize=(10,5))
#ax.plot(dates,trans_total,'*g',label='MOM6')
for nn in np.arange(len(ncfiles_3d)):
    ax.plot(time_total[nn,:],trans_total[nn,:],'-*',color=colors_exp[nn],label=labels_exp[nn])
ax.plot(fc_time,fc_transp,'o-k',markersize=6,linewidth=3,markeredgecolor='gold',label='Flor. Cable')
ax.fill_between(fc_time,fc_transp+fc_transp*0.05,fc_transp-fc_transp*0.05,color='grey',alpha=0.3)
ax.set_xticks(dates[11::24])
ax.set_xticklabels(labels[11::24])
ax.set_xlabel('UTC start date '+str(time_vec[0]),fontsize=14)
ax.set_ylabel('Transport (Sv)',fontsize=14)
ax.grid()
ax.set_title('Florida Cable Transport',fontsize=16)
ax.set_xlim(dates[0],dates[-1])
ax.legend()
#ax.legend(loc='lower left')

