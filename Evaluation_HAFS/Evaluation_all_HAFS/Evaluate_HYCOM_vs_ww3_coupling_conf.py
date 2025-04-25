#%% User input

# forecasting cycle to be used
cycle = '2020082506'

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

list_exp = ['a2o2a_a2w2a','a2o2a_a2w','a2w2a','a2w','abw']

folder_hafs_phase3_final = scratch_folder + 'HAFSv0p2a_phase3_final/' + cycle + '/00L/'
folder_hafs_phase3_wave = scratch_folder + 'HAFSv0p2a_phase3_wave/' + cycle + '/00L/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
coast_line_file = scratch_folder + 'Data/coastlines_180x180.dat'

best_track_file = scratch_folder + 'bdeck/bal132020.dat'
GFS_track_file = scratch_folder + 'adeck/aal132020.dat'

#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = home_folder + 'Analysis/Evaluation_HAFS/Evaluation_HAFSv0p2a_phase3_wave/Laura_2020/Figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob
import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_HYCOM,\
                            figure_transect_temp, glider_data_vector_to_array,\
                            grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
def coast180(file):
   inf=open(file,'r')
   c1 = []
   c2 = []
   hsk=np.genfromtxt(inf)
   c1=hsk[:,0]
   c2=hsk[:,1]
   return (c1,c2)

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

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

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read GFS track
lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)
 
#################################################################################
#%% Reading HYCOM grid
# Reading lat and lon
i=0
hycom_grid = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'forecast_'+list_exp[i]+'/','regional.grid.a')))[0][0:-2]

lines_grid = [line.rstrip() for line in open(hycom_grid+'.b')]
lon_hycom = np.array(readgrids(hycom_grid,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid,'plat:',[0]))

# Extracting the longitudinal and latitudinal size array
idm=int([line.split() for line in lines_grid if 'longitudinal' in line][0][0])
jdm=int([line.split() for line in lines_grid if 'latitudinal' in line][0][0])
'''
afiles = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'forecast_'+list_exp[i]+'/','*archv*.a')))

# Reading depths
lines = [line.rstrip() for line in open(afiles[0][:-2]+'.b')]
z = []
for line in lines[6:]:
    if line.split()[2]=='temp':
        print(line.split()[1])
        z.append(float(line.split()[1]))
depth_HYCOM = np.asarray(z)

nz = len(depth_HYCOM_exp)
'''
#################################################################################
n = 1 # 12 hours forecast

files_hafs_phase3_hycom = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'forecast_'+list_exp[i]+'/','*archv*.a')))

file = files_hafs_phase3_hycom[n]

#Reading time stamp
year = file.split('/')[-1].split('.')[1].split('_')[0]
day_of_year = file.split('/')[-1].split('.')[1].split('_')[1]
hour = file.split('/')[-1].split('.')[1].split('_')[2]

tyme_hycom.append(datetime.strptime(year + '-' + day_of_year + '-' + hour, "%Y-%j-%H"))

sst_hycom = readBinz(file[:-2],'3z','temp')[:,:,0]
sss_hycom = readBinz(file[:-2],'3z','salin')[:,:,0]
uvel_hycom = readBinz(file[:-2],'3z','u-vel.')[:,:,0]
vvel_hycom = readBinz(file[:-2],'3z','v-vel.')[:,:,0]
 
lonh,lath = HYCOM_coord_to_geo_coord(lon_hycom[0,:],lat_hycom[:,0])
'''   
    kw = dict(levels=np.arange(20,33,1))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lonh,lath,sst_hycom,cmap=cmocean.cm.thermal,**kw)
    cb = plt.colorbar()
    #plt.plot(lon_best_track,lat_best_track,'.-',color='grey',label='Best Track')
    #plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.legend(loc='lower left')
    plt.title('SST '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('$(^oC)$',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

    kw = dict(levels=np.arange(20,33,1))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lonh,lath,sss_hycom,cmap='Spectral_r',**kw)
    cb = plt.colorbar()
    #plt.plot(lon_best_track,lat_best_track,'.-',color='grey',label='Best Track')
    #plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.legend(loc='lower left')
    plt.title('SST '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('$(^oC)$',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])
'''



'''
#******************************************************************
n = 5 # 12 hours forecast
uwind_wave = np.empty((len(list_exp),n,441,901)
vwind_wave[:] = np.nan
vwind_wave = np.empty((len(list_exp),n,441,901))
vwind_wave[:] = np.nan
hs_wave = np.empty((len(list_exp),n,441,901))
hs_wave[:] = np.nan
fp_wave = np.empty((len(list_exp),n,441,901))
fp_wave[:] = np.nan
dp_wave = np.empty((len(list_exp),n,441,901))
dp_wave[:] = np.nan

for i in np.arange(len(list_exp)):

    files_hafs_phase3_wave = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'forecast_'+list_exp[i]+'/','*ww3.field*.nc')))
    print(files_hafs_phase3_wave)
    WW3 = xr.open_dataset(files_hafs_phase3_wave[0],decode_times=True)

    # read WW3 time
    t = np.asarray(WW3.variables['time'][:])
    timestamp = mdates.date2num(t)
    time = mdates.num2date(timestamp)

    lat = np.asarray(WW3.variables['latitude'])
    lon = np.asarray(WW3.variables['longitude'])

    mapsta = np.asarray(WW3.variables['MAPSTA'])
    uwind = np.asarray(WW3.variables['uwnd'][:])
    vwind = np.asarray(WW3.variables['vwnd'][:])
    hs = np.asarray(WW3.variables['hs'][:])
    fp = np.asarray(WW3.variables['fp'][:])
    dp = np.asarray(WW3.variables['dp'][:])

    uwind_wave[i,:,:,:] = uwind[0:n,:,:]
    vwind_wave[i,:,:,:] = vwind[0:n,:,:]
    hs_wave[i,:,:,:] = hs[0:n,:,:]
    fp_wave[i,:,:,:] = fp[0:n,:,:]
    dp_wave[i,:,:,:] = dp[0:n,:,:]

    lonm , latm = np.meshgrid(lon,lat)
    #okt = np.where(mdates.date2num(time_best_track)>=mdates.date2num(time[n]))[0][0]

    kw = dict(levels=np.arange(-32,33,2))
    plt.figure()
    plt.plot(cx,cy,'-k')
    #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,uwind[n-1,:,:],cmap='seismic',**kw)
    cb = plt.colorbar()
    #plt.plot(lon_best_track,lat_best_track,'.-',color='grey',label='Best Track')
    #plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.quiver(lonm[::20,::20],latm[::20,::20],uwind[n,::20,::20],vwind[n,::20,::20],scale=200)
    #plt.legend(loc='lower left')
    plt.title('Eastward Wind '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('m/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

    kw = dict(levels=np.arange(-32,33,2))
    plt.figure()
    plt.plot(cx,cy,'-k')
    #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,vwind[n-1,:,:],cmap='seismic',**kw)
    cb = plt.colorbar()
    #plt.plot(lon_best_track, lat_best_track,'.-',color='k',label='Best Track')
    #plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.quiver(lonm[::20,::20],latm[::20,::20],uwind[n,::20,::20],vwind[n,::20,::20],scale=200)
    #plt.legend(loc='lower left')
    plt.title('Northward Wind '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('m/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

    kw = dict(levels=np.arange(0,12.3,0.2))
    plt.figure()
    plt.plot(cx,cy,'-k')
    #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,hs[n-1,:,:],cmap='Spectral_r',**kw)
    cb = plt.colorbar()
    plt.title('Sea Surface Wave Significant Height \n '+list_exp[i]+' '+str(time[n-1])[0:13])
    cb.set_label('m',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

    kw = dict(levels=np.arange(0,1.05,0.05))
    plt.figure()
    plt.plot(cx,cy,'-k')
    #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,fp[n-1,:,:],cmap='Spectral_r',**kw)
    cb = plt.colorbar()
    plt.title('Dominat Wave Frequency \n'+list_exp[i]+' '+str(time[n-1])[0:13])
    cb.set_label('1/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

    kw = dict(levels=np.arange(0,361,10))
    plt.figure()
    plt.plot(cx,cy,'-k')
    #plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    #plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,dp[n-1,:,:],cmap='twilight',**kw)
    cb = plt.colorbar()
    plt.title('Sea Surface Wave Peak Direction \n'+list_exp[i]+' '+str(time[n-1])[0:13])
    cb.set_label('1/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])

ok_lon = np.where(lon >= -87)[0][0]

plt.figure()
plt.plot(lat, uwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, uwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon],'.-',label=list_exp[2],color='olivedrab')
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon],'.-',label=list_exp[3],color='darkcyan')
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon],'.-',label=list_exp[4],color='darkviolet')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)

plt.figure()
plt.plot(lat, vwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, vwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, vwind_wave[2,n-1,:,ok_lon],'.-',label=list_exp[2],color='olivedrab')
plt.plot(lat, vwind_wave[3,n-1,:,ok_lon],'.-',label=list_exp[3],color='darkcyan')
plt.plot(lat, vwind_wave[4,n-1,:,ok_lon],'.-',label=list_exp[4],color='darkviolet')
plt.legend()
plt.title('Northward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)

plt.figure()
plt.plot(lat, hs_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, hs_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, hs_wave[2,n-1,:,ok_lon],'.-',label=list_exp[2],color='olivedrab')
plt.plot(lat, hs_wave[3,n-1,:,ok_lon],'.-',label=list_exp[3],color='darkcyan')
plt.plot(lat, hs_wave[4,n-1,:,ok_lon],'.-',label=list_exp[4],color='darkviolet')
plt.legend()
plt.title('Wave Significant Height Longitudinal Transect',fontsize=14)
plt.ylabel('(m)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)

plt.figure()
plt.plot(lat, fp_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, fp_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, fp_wave[2,n-1,:,ok_lon],'.-',label=list_exp[2],color='olivedrab')
plt.plot(lat, fp_wave[3,n-1,:,ok_lon],'.-',label=list_exp[3],color='darkcyan')
plt.plot(lat, fp_wave[4,n-1,:,ok_lon],'.-',label=list_exp[4],color='darkviolet')
plt.legend()
plt.title('Dominant Wave Frequency Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)

plt.figure()
plt.plot(lat, dp_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, dp_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, dp_wave[2,n-1,:,ok_lon],'.-',label=list_exp[2],color='olivedrab')
plt.plot(lat, dp_wave[3,n-1,:,ok_lon],'.-',label=list_exp[3],color='darkcyan')
plt.plot(lat, dp_wave[4,n-1,:,ok_lon],'.-',label=list_exp[4],color='darkviolet')
plt.legend()
plt.title('Wave Peak Direction Longitudinal Transect',fontsize=14)
plt.ylabel('Degrees',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)

'''
