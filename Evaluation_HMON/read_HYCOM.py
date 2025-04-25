#%% User input
cycle = '2021082418'

# folder ab files HYCOM
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
hycom_folder = scratch_folder + 'HMOM2021_HYCOM/' + cycle + '/'

# RTOFS grid file name
hycom_grid_depth_folder = scratch_folder + 'HMOM2021_HYCOM/GRID_DEPTH/'
hycom_grid_file = hycom_grid_depth_folder + 'hmon_rtofs_hat10.basin.regional.grid'
hycom_depth_file = hycom_grid_depth_folder + 'hmon_rtofs_hat10.basin.regional.depth'

#%% 
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime
import os
import os.path
import glob
import cmocean

import sys
sys.path.append('/home/Maria.Aristizabal/Repos/NCEP_scripts/')
from utils4HYCOM2 import readBinz, readgrids, readdepth, readVar

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)


#%% Reading HYCOM grid
print('Retrieving coordinates from HYCOM')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid_file+'.b')]
lon_hycom = np.array(readgrids(hycom_grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid_file,'plat:',[0]))

depth_hycom = np.asarray(readdepth(hycom_depth_file,'depth'))

#%% finding list of hycom files
ahycom_files = sorted(glob.glob(os.path.join(hycom_folder,'*.a')))
bhycom_files = sorted(glob.glob(os.path.join(hycom_folder,'*.b')))

#%%%%%%%%%%%
#%% Reading hycom ab files
# Reading 3D variable from binary file
for n,file in enumerate(ahycom_files[::4]):
    temp_hycom = readBinz(file[:-2],'3z','temp')
    sst_hycom = temp_hycom[:,:,0]
    salt_hycom = readBinz(ahycom_files[n][:-2],'3z','salinity')
    sss_hycom = salt_hycom[:,:,0]

    #Reading time stamp
    year = int(file.split('/')[-1].split('.')[1][0:4])
    month = int(file.split('/')[-1].split('.')[1][4:6])
    day = int(file.split('/')[-1].split('.')[1][6:8])
    hour = int(file.split('/')[-1].split('.')[1][8:10])
    dt = int(file.split('/')[-1].split('.')[3][1:])
    timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
    time_hycom = mdates.num2date(timestamp_hycom)

    #%%
    kw = dict(levels = np.arange(15,33,1))
    fig, ax = plt.subplots()
    plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],sst_hycom,cmap='Spectral_r',**kw)
    plt.colorbar()
    cs = plt.contour(lon_hycom[0,:]-360,lat_hycom[:,0],sst_hycom,[26],colors='grey')
    #ax.clabel(cs, cs.levels, inline=True,fontsize=10)
    plt.title('SST HYCOM on '+str(time_hycom)[0:13],fontsize=14)
    plt.axis('scaled')

    kw = dict(levels = np.arange(32,38,0.3))
    plt.figure()
    plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],sss_hycom,cmap=cmocean.cm.haline,**kw)
    plt.colorbar()
    lb = plt.contour(lon_hycom[0,:]-360,lat_hycom[:,0],sss_hycom,[35,37],colors='grey')
    plt.title('SSS HYCOM on '+str(time_hycom)[0:13],fontsize=14)
    plt.axis('scaled')

##################################################################################
'''
nz = 41
layers = np.arange(0,nz)
#target_temp_hycom = np.empty(((tmax-tmin).days+1,nz,))
#target_temp_hycom[:] = np.nan
#target_zRTOFS = np.empty(((tmax-tmin).days+1,nz,))
#target_zRTOFS[:] = np.nan
#timeRTOFS = []
for n in np.arange(len(ahycom_files)):
    print(n)

    lines = [line.rstrip() for line in open(ahycom_files[n][:-2]+'.b')]    
    #time_stamp = lines[-1].split()[2]
    #hycom_days = lines[-1].split()[3]
    #tzero = datetime(1901,1,1,0,0)
    #timeRT = tzero+timedelta(float(hycom_days))
    #timeRTOFS.append(timeRT)
    #timestampRTOFS = mdates.date2num(timeRT) 
    time_str = lines[-1].split('date:')[1].split('[')[0]
    
    ztmp = readVar(ahycom_files[n][:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    #target_ztmp = ztmp[oklatRTOFS,oklonRTOFS]
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(afile[0][:-2],'archive',var_name,[lyr+1])
        #target_temp_RTOFS[tt,lyr] = temp_RTOFS[oklatRTOFS,oklonRTOFS]
        
        dp=readVar(afile[0][:-2],'archive','thknss',[lyr+1])/2/9806
        
        target_ztmp = np.append(target_ztmp,dp[oklatRTOFS,oklonRTOFS])
        
    target_z3d=np.cumsum(target_ztmp)              # [idm,jdm,kdm+1]
    target_z3d=np.squeeze(target_z3d[1:])             # [idm,jdm,kdm]
    target_z3d=np.array(target_z3d)
    target_z3d[target_z3d > 10**8] = np.nan
    target_zRTOFS[tt,:] = target_z3d
    
timeRTOFS = np.asarray(timeRTOFS)
'''
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
