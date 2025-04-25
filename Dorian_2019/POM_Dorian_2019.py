
#%% User input
# Huyn-Sook folder for Dorian
folder = '/scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/'

# POM grid file name
grid_file = folder + 'dorian05l.2019083018.pom.grid.nc'
 
# POM files
prefix = 'dorian05l.2019083018.pom.'

# Track file
prefix_track = 'dorian05l.2019083018.trak.hwrf.atcfunix'

# Name of 3D variable
var_name = 't'

# point of interest
target_lon = -77.4
target_lat = 27.0

#%% 
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import netCDF4 
from datetime import datetime,timedelta
from matplotlib.dates import date2num, num2date
import matplotlib.dates as mdates
import os
import os.path
import glob
import cmocean
from mpl_toolkits.basemap import Basemap

# GOFS 3.1 surface temp followind Dorian
surf_temp31_track = np.array([29.164001, 28.683   , 28.373001, 28.467   , 28.898   , 28.804   ,
       28.848999, 29.217   ,       np.nan, 27.981   ,       np.nan, 26.542   ,
       26.835001, 26.464   , 26.176   , 26.154   , 27.42    , 26.211   ,
       27.49    , 27.586   , 27.825   ])

#%% Reading POM grid files
pom_grid = xr.open_dataset( grid_file)
lonc = np.asarray(pom_grid['east_e'][:])
latc = np.asarray( pom_grid['north_e'][:])
zlevc = np.asarray(pom_grid['zz'][:])
topoz = np.asarray(pom_grid['h'][:])

#%% Read track file operational
file_track_oper = folder + prefix_track
ff_oper = open(file_track_oper,'r')
f_oper = ff_oper.readlines()

latt_oper = []
lonn_oper = []
lead_time_oper = []
for l in f_oper:
    lat = float(l.split(',')[6][0:4])/10
    if l.split(',')[6][4] == 'N':
        lat = lat
    else:
        lat = -lat
    lon = float(l.split(',')[7][0:5])/10
    if l.split(',')[7][4] == 'E':
        lon = lon
    else:
        lon = -lon
    latt_oper.append(lat)
    lonn_oper.append(lon)
    lead_time_oper.append(int(l.split(',')[5][1:4]))

latt_oper = np.asarray(latt_oper)
lonn_oper = np.asarray(lonn_oper)
lead_time_track_oper, ind = np.unique(lead_time_oper,return_index=True)
lat_track_oper = latt_oper[ind]
lon_track_oper = lonn_oper[ind]

# Getting time of track
t0 = datetime.strptime(prefix.split('.')[1],'%Y%m%d%H')  
time_track = [t0 + timedelta(hours=int(hrs)) for hrs in lead_time_track_oper] 

#%% Getting list of POM files
ncfiles = sorted(glob.glob(os.path.join(folder,prefix+'*0*.nc')))

# Reading POM time
time_pom = []
for i,file in enumerate(ncfiles):
    print(i)
    pom = xr.open_dataset(file)
    tpom = pom['time'][:]
    timestamp_pom = date2num(tpom)[0]
    time_pom.append(num2date(timestamp_pom))

time_pom = np.asarray(time_pom)

'''
# Reading POM temperature
target_temp_pom = np.empty((len(ncfiles),len(zlevc),))
target_temp_pom[:] = np.nan
target_topoz_pom = np.empty((len(ncfiles),))
target_topoz_pom[:] = np.nan
time_pom = []
for x,file in enumerate(ncfiles):
    print(x)
    pom = xr.open_dataset(file)

    tpom = pom['time'][:]
    timestamp_pom = date2num(tpom)[0]
    time_pom.append(num2date(timestamp_pom))

    # Interpolating latg and longlider into RTOFS grid
    #sublonpom = np.interp(timestamp_pom,timestampg,long)
    #sublatpom = np.interp(timestamp_pom,timestampg,latg)
    oklonpom = np.int(np.round(np.interp(target_lon,lonc[0,:],np.arange(len(lonc[0,:])))))
    oklatpom = np.int(np.round(np.interp(target_lat,latc[:,0],np.arange(len(latc[:,0])))))

    target_temp_pom[x,:] = np.asarray(pom['t'][0,:,oklatpom,oklonpom])
    target_topoz_pom[x] = np.asarray(topoz[oklatpom,oklonpom])

timestamp_pom = date2num(time_pom)

z_matrix_pom = np.dot(target_topoz_pom.reshape(-1,1),zlevc.reshape(1,-1))

#%% Figure Temp time series profile
time_matrixpom = np.tile(date2num(time_pom),(z_matrix_pom.shape[1],1)).T
kw = dict(levels = np.linspace(17,31,15))

fig, ax = plt.subplots(figsize=(11, 3))
plt.ion()
plt.contour(time_matrixpom,z_matrix_pom,target_temp_pom,colors = 'lightgrey',**kw)
plt.contour(time_matrixpom,z_matrix_pom,target_temp_pom,[26],colors = 'k')
plt.contourf(time_matrixpom,z_matrix_pom,target_temp_pom,cmap=cmocean.cm.thermal,**kw)
#ax.set_ylim(36,22.5)
#ax.set_xlim(datetime(2018,10,8),datetime(2018,10,13))
ax.set_ylim(-300,0)
yl = ax.set_ylabel('Depth (m)',fontsize=16) #,labelpad=20)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature ($^\circ$C)',fontsize=16)
ax.set_title('HWRF-POM',fontsize=20)
xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
ax.xaxis.set_major_formatter(xfmt)

file = "/home/Maria.Aristizabal/Dorian_2019/Figures/HWRF_pom_temp_Dorian.png"
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Surface temperature time series
fig,ax1 = plt.subplots(figsize=(11, 3))
plt.plot(timestamp_pom,target_temp_pom[:,0],'o-',linewidth=2)

t0 = datetime(2019,8,25)
deltat= timedelta(1)
xticks = [t0+nday*deltat for nday in np.arange(15)]
xticks = np.asarray(xticks)
plt.xticks(xticks)
xfmt = mdates.DateFormatter('%d \n %b')
ax1.xaxis.set_major_formatter(xfmt)
#plt.xlim([time31[okt,time31[oktime31[-1]]])
plt.ylabel('$^oC$',fontsize = 14)
tDorian = np.tile(datetime(2019,9,2,0),len(np.arange(22,30,0.1)))
plt.plot(tDorian,np.arange(22,30,0.1),'--k')
plt.title('Surface Temperature POM',fontsize=16)
plt.grid(True)
#plt.legend(loc='upper left',bbox_to_anchor=(0.6,1.3))

file = "/home/Maria.Aristizabal/Dorian_2019/Figures/surf_temp_POM_Dorian.png"
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
''' 

#%% Define categoty vector frem HWRF output based on max wind vel at 10 m
cat_hwrf = [None]*len(time_track)
cat_hwrf[0:3] = np.tile('cat3',3)
cat_hwrf[3:24] = np.tile('cat4',21)
cat_hwrf[24:27] = np.tile('cat3',3)
cat_hwrf[27:30] = np.tile('cat4',3)
cat_hwrf[30:36] = np.tile('cat3',6)
cat_hwrf[36] = 'cat2'
cat_hwrf[37:40] = np.tile('cat1',3)
cat_hwrf[40:43] = np.tile('ts',3)

time_pom = []
# Heat capacity in J/(kg K)
cp = 3985 
#%% Reading pom surface temperature
for i,file in enumerate(ncfiles):
    print(i)
    pom = xr.open_dataset(file)

    tpom = pom['time'][:]
    timestamp_pom = date2num(tpom)[0]
    time_pom.append(num2date(timestamp_pom))

    surf_temp_pom = np.asarray(pom['t'][0,0,:,:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])

    okt = np.where(date2num(time_pom[i]) == date2num(time_track))[0][0] 
    
    # OHC
    OHC_pom = np.empty((len(np.arange(botm,top)),len(np.arange(left,right))))
    OHC_pom[:] = np.nan
    for j, index in enumerate(np.arange(left,right)):
        for i, index in enumerate(np.arange(botm,top)):
            print(i,' ' ,j)
            ok26 = T31[:,i,j] >= 26
            rho0 = np.nanmean(D31[ok26,i,j])
            OHC[i,j] = cp * rho0 * np.trapz(T31[ok26,i,j]-26,df.depth[ok26])
    '''
    #%% map SST 
    m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
    x, y = m(*np.meshgrid(lonc[0,:],latc[:,0]))
    plt.figure()
    plt.ion()
    m.drawcoastlines()
    m.fillcontinents()
    m.drawmapboundary()
    kw = dict(levels=np.linspace(23,32,10))
    plt.contourf(x,y,surf_temp_pom,cmap=cmocean.cm.thermal,**kw)
    c = plt.colorbar()
    #q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
    #xq,yq = m(-78,12.5)
    #plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    xc, yc = m(-77.4,27.0)
    plt.plot(xc,yc,'*g',markersize=5)
    xt, yt = m(lon_track_oper,lat_track_oper)
    plt.plot(xt[okt],yt[okt],'or',markersize=5,label='Dorian Forecast '+cat_hwrf[okt])
    plt.plot(xt[0:okt],yt[0:okt],'o',color='grey',markersize=2)
    plt.title('POM SST on '+str(time_pom[i])[2:18],fontsize=16)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=16)
    c.ax.tick_params(labelsize=14)
    plt.legend()
    '''	

'''
target_topoz_pom = np.empty((len(ncfiles),))
target_topoz_pom[:] = np.nan
surf_temp_track = []
temp_prof_track = np.empty([len(zlevc),len(time_pom)])
temp_prof_track[:] = np.nan
#%% Reading pom temperature following Dorian
for i,file in enumerate(ncfiles):
    print(i)
    pom = xr.open_dataset(file)

    okt = np.where(date2num(time_pom[i]) == date2num(time_track))[0][0]
    print(okt)
    oklat = int(np.round(np.interp(lat_track_oper[okt],latc[:,0],np.arange(len(latc[:,0])))))
    oklon = int(np.round(np.interp(lon_track_oper[okt],lonc[0,:],np.arange(len(lonc[0,:])))))
    
    temp_prof_track[:,i] = pom['t'][0,:,oklat,oklon]
    surf_temp_track.append(np.asarray(pom['t'][0,0,oklat,oklon]))
    target_topoz_pom[i] = np.asarray(topoz[oklat,oklon])

z_matrix_pom = np.dot(target_topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
surf_temp_track = np.asarray(surf_temp_track)    
surf_temp_track[surf_temp_track == 0] = np.nan
temp_prof_track[temp_prof_track == 0] = np.nan

#%% Surface temperature following the storm
fig,ax1 = plt.subplots(figsize=(11, 3))
plt.plot(date2num(time_pom),surf_temp_track,'o-',color='mediumorchid',linewidth=2,label='HWRF-POM')
plt.plot(date2num(time_pom),surf_temp31_track,'o-',color='darkgreen',linewidth=2,label='GOFS 3.1')
xfmt = mdates.DateFormatter('%d \n %b')
ax1.xaxis.set_major_formatter(xfmt)
plt.ylabel('$^oC$',fontsize = 14)
plt.title('Surface Temperature POM following Dorian',fontsize=16)
plt.grid(True)
plt.legend(fontsize=14)
#plt.legend(loc='upper left',bbox_to_anchor=(0.6,1.3))

#%% Figure Temp time series profile
time_matrixpom = np.tile(date2num(time_pom),(z_matrix_pom.shape[0],1))
kw = dict(levels = np.linspace(17,31,15))

fig, ax = plt.subplots(figsize=(11, 3))
plt.ion()
plt.contour(time_matrixpom,z_matrix_pom,temp_prof_track,colors = 'lightgrey',**kw)
plt.contour(time_matrixpom,z_matrix_pom,temp_prof_track,[26],colors = 'k')
plt.contourf(time_matrixpom,z_matrix_pom,temp_prof_track,cmap=cmocean.cm.thermal,**kw)
#ax.set_ylim(36,22.5)
#ax.set_xlim(datetime(2018,10,8),datetime(2018,10,13))
ax.set_ylim(-300,0)
yl = ax.set_ylabel('Depth (m)',fontsize=16) #,labelpad=20)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Temperature ($^\circ$C)',fontsize=16)
ax.set_title('HWRF-POM following Dorian',fontsize=20)
xfmt = mdates.DateFormatter('%H:%Mh\n%d-%b')
ax.xaxis.set_major_formatter(xfmt)
'''
