#%% User input
cycle = '2020042500'

# folder ab files HYCOM
folder_RTOFS_DA = '/scratch2/NOS/nosofs/Maria.Aristizabal/RTOFS-DA/data_' + cycle
prefix_RTOFS_DA = 'archv'

# RTOFS grid file name
folder_RTOFS_DA_grid_depth = '/scratch2/NOS/nosofs/Maria.Aristizabal/RTOFS-DA/GRID_DEPTH/'
RTOFS_DA_grid = folder_RTOFS_DA_grid_depth + 'regional.grid'
RTOFS_DA_depth = folder_RTOFS_DA_grid_depth + 'regional.depth'

#%% 
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime
import os
import os.path
import glob

import sys
sys.path.append('/home/aristizabal/NCEP_scripts/')
from utils4HYCOM import readBinz, readgrids
from utils4HYCOM import readdepth

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)


#%% Reading RTOFS grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(RTOFS_DA_grid+'.b')]
lon_RTOFS_DA = np.array(readgrids(RTOFS_DA_grid,'plon:',[0]))
lat_RTOFS_DA = np.array(readgrids(RTOFS_DA_grid,'plat:',[0]))

depth_RTOFS_DA = np.asarray(readdepth(RTOFS_DA_depth,'depth'))

#%%%%%%%%%%%
#%% Reading RTOFS ab files
    
nz = 41
layers = np.arange(0,nz)
target_temp_RTOFS = np.empty(((tmax-tmin).days+1,nz,))
target_temp_RTOFS[:] = np.nan
target_zRTOFS = np.empty(((tmax-tmin).days+1,nz,))
target_zRTOFS[:] = np.nan
timeRTOFS = []
for tt in np.arange(0,(tmax-tmin).days+1):
    print(tt)
    t = tmin+timedelta(np.int(tt))
    if t.day < 9: 
        afile = sorted(glob.glob(os.path.join(Dir_rtofs+str(tmin.year)+str(tmin.month)+'0'+str(tmin.day+1+tt)+'/',prefix_ab+'*.a')))
    else:
        afile = sorted(glob.glob(os.path.join(Dir_rtofs+str(tmin.year)+str(tmin.month)+str(tmin.day+1+tt)+'/',prefix_ab+'*.a')))
      
    lines = [line.rstrip() for line in open(afile[0][:-2]+'.b')]    
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero=datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days))
    timeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT) 

    #oklonRTOFS = 2508
    #oklatRTOFS = 1858
    
    # Interpolating latg and long into RTOFS grid
    sublonRTOFS = np.interp(timestampRTOFS,timestampg,target_lon)
    sublatRTOFS = np.interp(timestampRTOFS,timestampg,target_lat)
    oklonRTOFS = np.int(np.round(np.interp(sublonRTOFS,hlon[0,:],np.arange(len(hlon[0,:])))))
    oklatRTOFS = np.int(np.round(np.interp(sublatRTOFS,hlat[:,0],np.arange(len(hlat[:,0])))))
    
    ztmp=readVar(afile[0][:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklatRTOFS,oklonRTOFS]
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(afile[0][:-2],'archive',var_name,[lyr+1])
        target_temp_RTOFS[tt,lyr] = temp_RTOFS[oklatRTOFS,oklonRTOFS]
        
        dp=readVar(afile[0][:-2],'archive','thknss',[lyr+1])/2/9806
        
        target_ztmp = np.append(target_ztmp,dp[oklatRTOFS,oklonRTOFS])
        
    target_z3d=np.cumsum(target_ztmp)              # [idm,jdm,kdm+1]
    target_z3d=np.squeeze(target_z3d[1:])             # [idm,jdm,kdm]
    target_z3d=np.array(target_z3d)
    target_z3d[target_z3d > 10**8] = np.nan
    target_zRTOFS[tt,:] = target_z3d
    
timeRTOFS = np.asarray(timeRTOFS)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%





'''
# Extracting the longitudinal and latitudinal size array
idm=int([line.split() for line in lines_grid if 'longitudinal' in line][0][0])
jdm=int([line.split() for line in lines_grid if 'latitudinal' in line][0][0])

afiles = sorted(glob.glob(os.path.join(folder_RTOFS_DA,prefix_RTOFS_DA+'*.a')))

# Reading depths
lines=[line.rstrip() for line in open(afiles[0][:-2]+'.b')]
z = []
for line in lines[10:]:
    if line.split()[0]=='temp':
        print(line)
        z.append(float(line.split()[1]))
depth_RTOFS-DA = np.asarray(z) 

nz = len(depth_RTOFS_DA) 
'''
N = 0
var = 'temp'

afiles = sorted(glob.glob(os.path.join(folder_RTOFS_DA,prefix_RTOFS_DA+'*.a')))    
file = afiles[N]

# Reading 3D variable from binary file 
var_rtofs = readBinz(file[:-2],'3z',var)
temp_RTOFS_DA = var_rtofs[:,:,0]

#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_RTOFS_DA = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_RTOFS_DA = mdates.num2date(timestamp_RTOFS_DA)

#%%
kw = dict(levels = np.linspace(28,30,41))
plt.figure()
plt.contourf(lon_RTOFS_DA[0,:]-360,lat_RTOFS_DA[:,0],temp_RTOFS_DA,cmap=cmocean.cm.thermal,**kw)
plt.colorbar()
plt.title('RTOFS-DA SST \n on '+str(time_RTOFS_DA)[0:13],fontsize=14)
plt.axis('scaled')
