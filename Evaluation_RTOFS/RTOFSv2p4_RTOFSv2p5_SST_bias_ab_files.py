#%% User input

# forecasting cycle to be used
#cycles = ['20240901','20240902','20240903','20240904','20240905','20240905','20240906','20240907','20240908','20240909','20240910','20240911','20240912','20240913','20240914']
cycles = ['20240901','20240902','20240903','20240904','20240905','20240905','20240906','20240907','20240908','20240909','20240910','20240911','20240912','20240913','20240914','20240915','20240916','20240917','20240918','20240919','20240920','20240921','20240922','20240923','20240923','20240924','20240925','20240926','20240927','20240928','20240929','20240930','20241001','20241001','20241002','20241003','20241004','20241005','20241006','20241007','20241008','20241009','20241010','20241011','20241012','20241013','20241014','20241015']

scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
folder_rtofs1 = scratch_folder + 'RTOFS/'
folder_rtofs2 = scratch_folder + 'RTOFS_v2.5.test01/'

forec_hours = ['n00','f06','f12','f18']
prefix = 'rtofs_glo.t00z.'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

var_name = 'temp'
klayer = '1'
colormap='turbo'
#colormap='Spectral_r'
#colormap='jet'
min_val = 25
max_val = 32.4
delta_val = 0.4  # delta in colorbar
delta_contour = 1 # delta in contour plot
units = '$^oC$'
#lon_lim = [-180,180]
#lat_lim = [-90,90]

# North Atlantic
lon_lim = [-100,-60]
lat_lim = [0,40]

# West Pacific
#lon_lim = [-180,-90]
#lat_lim = [-20,45]

##################################################################################
def read_field_klayer(rtofs_file,var_name,klayer):

    lines = [line.rstrip() for line in open(rtofs_file+'.b')]
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)
    fld = ma.array([],fill_value=1.2676506002282294e+30)

    inFile = rtofs_file + '.a'

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == 'field':
                nheading = n + 1
                print(line.split()[0])

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == var_name and line.split()[4] == klayer:
                nvar = n - nheading
                print(nvar)
                print(n)

    fid = open(inFile,'rb')
    #fid.seek((nvar-1)*4*(npad+ijdm),0)
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    fld = np.array(fld)
    fld = ma.reshape(fld,(jdm,idm))

    return fld

##################################################################################
#%% 
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import glob
import sys
import numpy.ma as ma
import struct

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

sys.path.append(folder_myutils)
from my_models_utils import geo_coord_to_HYCOM_coord

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#################################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]


oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

##################################################################################
#%% Reading HYCOM grid from ab files
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(rtofs_grid_file+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
lat = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth = np.asarray(readdepth(rtofs_depth_file,'depth'))

lon_min_hycom, _  = geo_coord_to_HYCOM_coord(lon_lim[0],lat_lim[0])
lon_max_hycom, _  = geo_coord_to_HYCOM_coord(lon_lim[1],lat_lim[1])

oklon = np.logical_and(lon[0,:] >= lon_min_hycom,lon[0,:] <= lon_max_hycom)
oklat = np.logical_and(lat[:,0] >= lat_lim[0],lat[:,0] <= lat_lim[1])

lonn = lon[oklat,:][:,oklon]
latt = lat[oklat,:][:,oklon]

##################################################################################
print('Reading ab files')

time = []
timestamp = []
mean_sst1 = []
std_sst1 = []
mean_sst2 = []
std_sst2 = []
mean_bias = []
std_bias = []
bias_SST = np.empty((len(cycles),len(forec_hours),latt.shape[0],latt.shape[1]))
bias_SST[:] = np.nan
for c,cycle in enumerate(cycles):
    print(cycle)
    for f,fh in enumerate(forec_hours):
        print(fh)
        file1 = glob.glob(folder_rtofs1 + 'rtofs.' + cycle + '/' + prefix + fh + '.archv.a')[0]
        file2 = glob.glob(folder_rtofs2 + 'rtofs.' + cycle + '/' + prefix + fh + '.archv.a')[0]

        lines = [line.rstrip() for line in open(file1[:-1]+'b')]
        time_stamp = lines[-1].split()[2]
        hycom_days = lines[-1].split()[3]
        tzero = datetime(1901,1,1,0,0)
        timeRT = tzero+timedelta(float(hycom_days)-1)
        time.append(timeRT)
        timestamp.append(mdates.date2num(timeRT))
    
        fld1 = read_field_klayer(file1[:-2],var_name,klayer)
        #fld1 = np.copy(fldd1)
        mask = fld1 > 10**5
        fld1[mask] = np.nan
        SST1 = fld1[:,oklon][oklat,:]
        mean_sst1.append(np.nanmean(SST1))
        std_sst1.append(np.nanstd(SST1))
    
        fld2 = read_field_klayer(file2[:-2],var_name,klayer)
        #fld1 = np.copy(fldd1)
        mask = fld2 > 10**5
        fld2[mask] = np.nan
        SST2 = fld2[:,oklon][oklat,:]
        mean_sst2.append(np.nanmean(SST2))
        std_sst2.append(np.nanstd(SST2))
    
        bias_sst = SST2 - SST1
        bias_SST[c,f,:,:] = bias_sst
        mean_bias.append(np.nanmean(bias_sst))
        std_bias.append(np.nanstd(bias_sst))
    
        '''
        cflevels = np.arange(20,32,1)
        plt.figure()
        #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn-360,latt,SST1,levels=cflevels,cmap='RdBu_r',extend='both')
        plt.colorbar(extendrect=True)
        plt.axis('scaled')
        plt.ylim(lat_lim)
        plt.xlim(lon_lim)
        plt.title('SST Diff RTOFSv2.5 - RTOFSv2.4  '+ str(timeRT)[0:13])

        cflevels = np.arange(-2,2.1,0.2)
        plt.figure()
        #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
        #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
        plt.contourf(lonn-360,latt,bias_sst,levels=cflevels,cmap='RdBu_r',extend='both')
        plt.colorbar(extendrect=True)
        plt.axis('scaled')
        plt.ylim(lat_lim)
        plt.xlim(lon_lim)
        plt.title('SST Diff RTOFSv2.5 - RTOFSv2.4  '+ str(str(timeRT)[0:13])[0:13])
    
        plt.figure()
        plt.plot(SST1,SST2,'.',color='green')
        plt.plot(np.arange(-5,36),np.arange(-5,36),color='grey',linewidth=3)
        plt.xlabel('RTOFSv2.4')
        plt.ylabel('RTOFSv2.5')
        plt.title(str(timeRT)[0:13])
        plt.xlim([18,35])
        plt.ylim([18,35])
        '''

#Bias_SST = np.reshape(bias_SST,(48*4,1710,742))
Bias_SST = np.reshape(bias_SST,(len(cycles)*len(forec_hours),latt.shape[0],latt.shape[1]))
time_mean_Bias_SST = np.nanmean(Bias_SST,axis=0)
Bias_SST_std = np.nanstd(Bias_SST,axis=0)

cflevels = np.arange(-1,1.1,0.2)
plt.figure()
#plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
#plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lonn-360,latt,time_mean_Bias_SST,levels=cflevels,cmap='RdBu_r',extend='both')
plt.colorbar(extendrect=True)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('SST Bias RTOFSv2.5 - RTOFSv2.4')

cflevels = np.arange(0,1.1,0.1)
plt.figure()
#plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
#plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.contourf(lonn-360,latt,Bias_SST_std,levels=cflevels,cmap='RdBu_r',extend='max')
plt.colorbar(extendrect=True)
plt.axis('scaled')
plt.ylim(lat_lim)
plt.xlim(lon_lim)
plt.title('STD: SST Bias RTOFSv2.5 - RTOFSv2.4')

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(time,mean_sst1,'.-',color='magenta',label='RTOFSv2.4')
plt.plot(time,mean_sst2,'.-',color='salmon',label='RTOFSv2.5')
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.legend()
plt.grid(True)
plt.title('Domain Mean SST')
plt.ylabel('$^oC$',fontsize=14)

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(time,np.zeros(len(time)),'-k')
plt.plot(time,mean_bias,'.-',color='cyan')
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.ylim([-0.2,0.2])
plt.grid(True)
plt.title('SST Bias RTOFSv2.5 - RTOFSv2.4')
plt.ylabel('$^oC$',fontsize=14)


