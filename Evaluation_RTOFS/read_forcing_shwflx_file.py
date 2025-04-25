#%% User input

# folder ab files HYCOM
save_folder = '/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
#hycom_fix_folder = save_folder + 'hafs_develop_202105/fix/fix_hycom/'
#hycom_fix_folder = '/scratch1/NCEPDEV/hwrf/noscrub/hafs-fix-files/hafs-20210520-fix/fix/fix_hycom/'
#hycom_fix_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/SubDomain_v11/NHC/for_fix/'
#hycom_fix_folder = '/scratch2/NCEPDEV/marine/Hyun.Sook.Kim/SubDomain_v11/NHC/for_fix/'
#hycom_fix_folder = '/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/hafs_develop_202105/fix/fix_hycom/'
hycom_fix_folder = '/scratch1/NCEPDEV/hwrf/scrub/Maria.Aristizabal/hafsv1_final_test_emc_graphics_hfsa/2020082512/13L/ocn_prep/hycominit2/'


basin = 'nhc'
#basin = 'hat10'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
coastline_file = scratch_folder +'Data/coastlines_180x180.dat'

##################################################################################
def coast180(costline_file):
   inf=open(coastline_file,'r')
   c1 = []
   c2 = []
   hsk=np.genfromtxt(inf)
   c1=hsk[:,0]
   c2=hsk[:,1]
   return (c1,c2)

##################################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
def find_grid_position_hycom(lon,lat,target_lon,target_lat):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),\
                                             left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklon,oklat

##################################################################################
def get_profile_rtofs_ab_file_desn_layers(nz,oklon,oklat,afile,bfile):

    #nz = 41
    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz))
    target_temp_rtofs[:] = np.nan
    target_z_rtofs = np.empty((nz))
    target_z_rtofs[:] = np.nan
    timeRTOFS = []

    lines = [line.rstrip() for line in open(bfile)]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days))
    timeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT)

    ztmp = readVar(afile[:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklat,oklon]
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(afile[:-2],'archive','temp',[lyr+1])
        target_temp_rtofs[lyr] = temp_RTOFS[oklat,oklon]
        dp = readVar(afile[:-2],'archive','thknss',[lyr+1])/2/9806
        target_ztmp = np.append(target_ztmp,dp[oklat,oklon])

    target_z3d = np.cumsum(target_ztmp)              # [idm,jdm,kdm+1]
    target_z3d = np.squeeze(target_z3d[1:])             # [idm,jdm,kdm]
    target_z3d = np.asarray(target_z3d)
    target_z3d[target_z3d > 10**8] = np.nan
    target_z_rtofs = target_z3d

    time_rtofs = np.asarray(timeRTOFS)

    return target_temp_rtofs, target_z_rtofs, time_rtofs

##################################################################################
#%% 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
#import xarray as xr
import os
import glob
import sys
#import cmocean

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar, parse_b

sys.path.append(folder_myutils)
#from my_models_utils import geo_coord_to_HYCOM_coord

# Increase fontsize of labels globally
#plt.rc('xtick',labelsize=14)
#plt.rc('ytick',labelsize=14)
#plt.rc('legend',fontsize=14)

##################################################################################
#%% hycom files
#hycom_depth_file = hycom_fix_folder + 'hafs_hycom_'+basin+'.basin.regional.depth'
hycom_depth_file = hycom_fix_folder + 'regional.depth'

hycom_grid_file = hycom_fix_folder + 'regional.grid'

#hycom_forc_chl_file = hycom_fix_folder + 'hafs_hycom_nhc.basin.forcing.chl'

hycom_forc_chl_file = hycom_fix_folder + 'forcing.shwflx'

#################################################################################
#%% Reading bathymetry data
#ncbath = xr.open_dataset(bath_file)
#bath_lat = ncbath.variables['lat'][:]
#bath_lon = ncbath.variables['lon'][:]
#bath_elev = ncbath.variables['elevation'][:]

'''
oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]
'''

#################################################################################
#%% Reading coastline data

cx,cy = coast180(coastline_file)
#################################################################################
#%% Reading HYCOM grid  and depth from ab files
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid_file+'.b')]
lon_hycom = np.array(readgrids(hycom_grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid_file,'plat:',[0]))

depth_hycom = np.asarray(readdepth(hycom_depth_file,'depth'))

'''
kw = dict(levels=np.arange(-7000,1,100)) 
plt.figure()
plt.plot(cx,cy,'-k')
#plt.plot(lon_hycom[::25,::25]-360,lat_hycom[::25,::25],'.k',markersize=1)
plt.contourf(lon_hycom,lat_hycom,-depth_hycom,cmap='YlGnBu_r',**kw,vmax=0)
cl = plt.colorbar(shrink = 0.7)
cl.set_label('Meters',fontsize=14)
plt.axis('scaled')
plt.xlim([lon_hycom.min(),lon_hycom.max()])
plt.ylim([lat_hycom.min(),lat_hycom.max()])
plt.title('HYCOM Bathymetry',fontsize=14)
'''

##################################################################################
#%% Reading HYCOM forcing ab files
import numpy as np
import numpy.ma as ma
import struct

lines = [line.rstrip() for line in open(hycom_forc_chl_file+'.b')]

idm = int([line.split() for line in lines if 'i/jdm' in line][0][2])
jdm = int([line.split() for line in lines if 'i/jdm' in line][0][3])
ijdm = idm*jdm
npad = 4096 - (ijdm%4096)
fld2 = ma.array([],fill_value=1e30)

inFile = hycom_forc_chl_file+'.a'

fid = open(inFile,'rb')
# read first field 
#nvar = 0
nvar = 3
fid.seek((nvar)*4*(npad+ijdm),0)
fld=fid.read(ijdm*4)
fld=struct.unpack('>'+str(ijdm)+'f',fld)
fld=np.array(fld)
fld=ma.reshape(fld,(jdm,idm))

print(np.max(fld))
print(np.min(fld))

end_heading = 'i/jdm'
for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == end_heading:
            nheading = n + 1
            print(line.split()[0])

# Read field time series
var_name = 'shwflx:'
oklon = np.where(lon_hycom[0,:] - 360 >= 0.0)[0][0]
oklat = np.where(lat_hycom[:,0] >= 0.0)[0][0]
fld_time_series = []
time_stamp = []
for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == var_name:
            nvar = n - nheading
            print(nvar)
            fid.seek((nvar)*4*(npad+ijdm),0)
            fld=fid.read(ijdm*4)
            fld=struct.unpack('>'+str(ijdm)+'f',fld)
            fld=np.array(fld)
            fld=ma.reshape(fld,(jdm,idm))
            print(np.max(fld))
            print(np.min(fld))
            fld_time_series.append(fld[oklat,oklon])
            time_stamp.append(float(lines[nvar+nheading].split('=')[1].split()[0]))

fld_time_series = np.asarray(fld_time_series)
time_stamp = np.asarray(time_stamp)
#time = np.asarray([datetime(1900,1,1) - timedelta(days = 365 + 30 + 2) + timedelta(days = t) for t in time_stamp])
time = np.asarray([datetime(1901,1,1) - timedelta(days = 1) + timedelta(days = t) for t in time_stamp])

# Read field for only one time
nvar = 3
fid.seek((nvar)*4*(npad+ijdm),0)
fld=fid.read(ijdm*4)
fld=struct.unpack('>'+str(ijdm)+'f',fld)
fld=np.array(fld)
fld=ma.reshape(fld,(jdm,idm))
print(np.max(fld))
print(np.min(fld))

# Plot field all damain
kw = dict(levels=np.arange(0,1051,150)) 
fig,ax1 = plt.subplots(figsize = (10,7))
ax1.plot(cx,cy,'-k')
plt.contourf(lon_hycom-360,lat_hycom,fld,**kw) #,cmap='Spectral_r') #,**kw)
plt.axis('scaled') #, adjustable='box', anchor='C'). 
plt.colorbar()
plt.xlim([-180,15])
plt.ylim([-25,46])
plt.title(str(time[nvar]))
#ax1.set_xlim([-100,-5])
#ax1.set_ylim([0,48])
#ax1.plot(np.arange(-10,-7.4),np.tile(2.4,len(np.arange(-10,-7.4))),'-k')
#ax1.plot(np.tile(-10,len(np.arange(1,2.4))),np.arange(1,2.4),'-k')
#ax1.text(-100,50,'(a)',fontsize=14) 

dates = mdates.date2num(time)
labels = np.asarray([t.hour for t in time])
oklab = labels == 12
fig,ax = plt.subplots(figsize=(10,5))
ax.plot(time,fld_time_series,'.-k')
ax.set_xticks(dates[oklab])
ax.set_xticklabels(labels[oklab])
ax.set_xlabel('UTC start date '+str(time[0]),fontsize=14)
ax.grid()

##################################################################################
