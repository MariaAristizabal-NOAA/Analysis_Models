#%% User input

# forecasting cycle to be used
cycle = '2020082506'

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

list_exp = ['a2o2a_a2w2a','a2o2a_a2w','a2o_a2w2a','a2o_a2w','abo_abw','a2w2a','a2w','w2a','abw']

#folder_hafs_phase3_final = scratch_folder + 'HAFSv0p2a_phase3_final/' + cycle + '/00L/'
folder_hafs_phase3_wave = scratch_folder + 'HAFSv0p2a_phase3_wave/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
coast_line_file = scratch_folder + 'Data/coastlines_180x180.dat'

best_track_file = scratch_folder + 'bdeck/bal132020.dat'
GFS_track_file = scratch_folder + 'adeck/aal132020.dat'

#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = scratch_folder + 'figures/'

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
 
#******************************************************************
n = 4 # 12 hours forecast
uwind_wave = np.empty((len(list_exp),n,441,901))
uwind_wave[:] = np.nan
vwind_wave = np.empty((len(list_exp),n,441,901))
vwind_wave[:] = np.nan
hs_wave = np.empty((len(list_exp),n,441,901))
hs_wave[:] = np.nan
fp_wave = np.empty((len(list_exp),n,441,901))
fp_wave[:] = np.nan
dp_wave = np.empty((len(list_exp),n,441,901))
dp_wave[:] = np.nan

for i in np.arange(len(list_exp)):

    files_hafs_phase3_wave = sorted(glob.glob(os.path.join(folder_hafs_phase3_wave+'/'+cycle+'/00L/'+'hafs_202108_'+list_exp[i]+'/com/'+cycle+'/00L/','*ww3_ounf.nc')))
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

    uwind_wave[i,:,:,:] = uwind[0:n,:,:]
    vwind_wave[i,:,:,:] = vwind[0:n,:,:]
    hs_wave[i,:,:,:] = hs[0:n,:,:]
    fp_wave[i,:,:,:] = fp[0:n,:,:]

    lonm , latm = np.meshgrid(lon,lat)
    okt = np.where(mdates.date2num(time_best_track)>=mdates.date2num(time[n]))[0][0]

    '''
    kw = dict(levels=np.arange(-32,33,2))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,uwind[n-1,:,:],cmap='seismic',**kw)
    cb = plt.colorbar()
    plt.plot(lon_best_track,lat_best_track,'.-',color='grey',label='Best Track')
    plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.quiver(lonm[::20,::20],latm[::20,::20],uwind[n,::20,::20],vwind[n,::20,::20],scale=200)
    #plt.legend(loc='lower left')
    plt.title('Eastward Wind '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('m/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])
    file_name = folder_fig + 'uwind_' + list_exp[i]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    kw = dict(levels=np.arange(-32,33,2))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,vwind[n-1,:,:],cmap='seismic',**kw)
    cb = plt.colorbar()
    plt.plot(lon_best_track, lat_best_track,'.-',color='k',label='Best Track')
    plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    #plt.quiver(lonm[::20,::20],latm[::20,::20],uwind[n,::20,::20],vwind[n,::20,::20],scale=200)
    #plt.legend(loc='lower left')
    plt.title('Northward Wind '+list_exp[i]+' \n'+str(time[n-1])[0:13])
    cb.set_label('m/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])
    file_name = folder_fig + 'vwind_' + list_exp[i]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    kw = dict(levels=np.arange(0,12.3,0.2))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,hs[n-1,:,:],cmap='Spectral_r',**kw)
    cb = plt.colorbar()
    plt.plot(lon_best_track, lat_best_track,'.-',color='k',label='Best Track')
    plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    plt.title('Sea Surface Wave Significant Height \n '+list_exp[i]+' '+str(time[n-1])[0:13])
    cb.set_label('m',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])
    file_name = folder_fig + 'hs_' + list_exp[i]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    kw = dict(levels=np.arange(0,1.05,0.05))
    plt.figure()
    #plt.plot(cx,cy,'-k')
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contourf(lon,lat,fp[n-1,:,:],cmap='Spectral_r',**kw)
    cb = plt.colorbar()
    plt.plot(lon_best_track, lat_best_track,'.-',color='k',label='Best Track')
    plt.plot(lon_best_track[okt],lat_best_track[okt],'.-',color='lime')
    plt.title('Dominat Wave Frequency \n'+list_exp[i]+' '+str(time[n-1])[0:13])
    cb.set_label('1/s',fontsize=14)
    plt.axis('scaled')
    plt.ylim([0,45])
    plt.xlim([-100,-10])
    file_name = folder_fig + 'fp_' + list_exp[i]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)
    '''

ok_lon = np.where(lon >= -87)[0][0]
ok_lat = np.where(lat >= 24.2)[0][0]

plt.figure()
plt.plot(lat, uwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, uwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, uwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet' 
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[0,n-1,ok_lat,:],'.-',label=list_exp[0],color='firebrick')
plt.plot(lon, uwind_wave[1,n-1,ok_lat,:],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lon, uwind_wave[2,n-1,ok_lat,:],'o-',label=list_exp[2],color='chocolate')
plt.plot(lon, uwind_wave[3,n-1,ok_lat,:],'s-',label=list_exp[3],color='orangered')
plt.plot(lon, uwind_wave[4,n-1,ok_lat,:],'^-',label=list_exp[4],color='gold')
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lon, uwind_wave[6,n-1,ok_lat,:],'.-',label=list_exp[6],color='lime')
plt.plot(lon, uwind_wave[7,n-1,ok_lat,:],'.-',label=list_exp[7],color='c')
plt.plot(lon, uwind_wave[8,n-1,ok_lat,:],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

###################################################################################
#%% figures for different groups

#%% atm-ocean atm-wave
plt.figure()
plt.plot(lat, uwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, uwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet2'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[0,n-1,ok_lat,:],'.-',label=list_exp[0],color='firebrick')
plt.plot(lon, uwind_wave[1,n-1,ok_lat,:],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lon, uwind_wave[2,n-1,ok_lat,:],'o-',label=list_exp[2],color='chocolate')
plt.plot(lon, uwind_wave[3,n-1,ok_lat,:],'s-',label=list_exp[3],color='orangered')
plt.plot(lon, uwind_wave[4,n-1,ok_lat,:],'^-',label=list_exp[4],color='gold')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat2'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% atm-wave
plt.figure()
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, uwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet3'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lon, uwind_wave[6,n-1,ok_lat,:],'.-',label=list_exp[6],color='lime')
plt.plot(lon, uwind_wave[7,n-1,ok_lat,:],'.-',label=list_exp[7],color='c')
plt.plot(lon, uwind_wave[8,n-1,ok_lat,:],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat3'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% atm-wave
plt.figure()
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon]-uwind_wave[6,n-1,:,ok_lon],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet8'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lon, uwind_wave[6,n-1,ok_lat,:],'.-',label=list_exp[6],color='lime')
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:]-uwind_wave[6,n-1,ok_lat,:],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat8'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% atm-wave
plt.figure()
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, uwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon]-uwind_wave[8,n-1,:,ok_lon],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet9'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[7,n-1,ok_lat,:],'.-',label=list_exp[7],color='c')
plt.plot(lon, uwind_wave[8,n-1,ok_lat,:],'.-',label=list_exp[8],color='deepskyblue')
plt.plot(lon, uwind_wave[7,n-1,ok_lat,:]-uwind_wave[8,n-1,ok_lat,:],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat9'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% atm-wave
plt.figure()
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon]-uwind_wave[7,n-1,:,ok_lon],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet10'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lon, uwind_wave[7,n-1,ok_lat,:],'.-',label=list_exp[7],color='c')
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:]-uwind_wave[7,n-1,ok_lat,:],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat10'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% equivalent conf for the wave model
plt.figure()
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon]-uwind_wave[5,n-1,:,ok_lon],'o-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet4'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[2,n-1,ok_lat,:],'o-',label=list_exp[2],color='chocolate')
plt.plot(lon, uwind_wave[5,n-1,ok_lat,:],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lon, uwind_wave[2,n-1,ok_lat,:]-uwind_wave[5,n-1,ok_lat,:],'o-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat4'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% equivalent conf for the wave model
plt.figure()
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, uwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon]-uwind_wave[6,n-1,:,ok_lon],'s-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet6'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[3,n-1,ok_lat,:],'s-',label=list_exp[3],color='orangered')
plt.plot(lon, uwind_wave[6,n-1,ok_lat,:],'.-',label=list_exp[6],color='lime')
plt.plot(lon, uwind_wave[3,n-1,ok_lat,:]-uwind_wave[6,n-1,ok_lat,:],'.-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat6'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


#%% equivalent conf for the ocean model
plt.figure()
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, uwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon]-uwind_wave[8,n-1,:,ok_lon],'^-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet7'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lon, uwind_wave[4,n-1,ok_lat,:],'^-',label=list_exp[4],color='gold')
plt.plot(lon, uwind_wave[8,n-1,ok_lat,:],'.-',label=list_exp[8],color='deepskyblue')
plt.plot(lon, uwind_wave[4,n-1,ok_lat,:]-uwind_wave[8,n-1,ok_lat,:],'^-',label='Diff',color='k')
plt.legend()
plt.title('Eastward Wind Latitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'uwind_transet_lat7'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)


##########################################################################
plt.figure()
plt.plot(lat, uwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, uwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, uwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, uwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, uwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, uwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, uwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, uwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, uwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Eastward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
plt.xlim([21.5,24])
plt.ylim([0,20])
file_name = folder_fig + 'uwind_transet'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, vwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, vwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, vwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, vwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, vwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, vwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, vwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, vwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, vwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Northward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'vwind_transet' 
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, vwind_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, vwind_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, vwind_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, vwind_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, vwind_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, vwind_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, vwind_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, vwind_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, vwind_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Northward Wind Longitudinal Transect',fontsize=14)
plt.ylabel('(m/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
plt.xlim([21.5,24])
file_name = folder_fig + 'vwind_transet_zoom'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, hs_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, hs_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, hs_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, hs_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, hs_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, hs_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, hs_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, hs_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, hs_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Wave Significant Height Longitudinal Transect',fontsize=14)
plt.ylabel('(m)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'hs_transet' 
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, hs_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, hs_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, hs_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, hs_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, hs_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, hs_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, hs_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, hs_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, hs_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Wave Significant Height Longitudinal Transect',fontsize=14)
plt.ylabel('(m)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.xlim([21.5,24.5])
plt.grid(True)
file_name = folder_fig + 'hs_transet_zoom'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, fp_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, fp_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, fp_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, fp_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, fp_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, fp_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, fp_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, fp_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, fp_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Dominant Wave Frequency Longitudinal Transect',fontsize=14)
plt.ylabel('(1/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
file_name = folder_fig + 'fp_transet' 
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

plt.figure()
plt.plot(lat, fp_wave[0,n-1,:,ok_lon],'.-',label=list_exp[0],color='firebrick')
plt.plot(lat, fp_wave[1,n-1,:,ok_lon],'.-',label=list_exp[1],color='lightcoral')
plt.plot(lat, fp_wave[2,n-1,:,ok_lon],'o-',label=list_exp[2],color='chocolate')
plt.plot(lat, fp_wave[3,n-1,:,ok_lon],'s-',label=list_exp[3],color='orangered')
plt.plot(lat, fp_wave[4,n-1,:,ok_lon],'^-',label=list_exp[4],color='gold')
plt.plot(lat, fp_wave[5,n-1,:,ok_lon],'.-',label=list_exp[5],color='darkgreen')
plt.plot(lat, fp_wave[6,n-1,:,ok_lon],'.-',label=list_exp[6],color='lime')
plt.plot(lat, fp_wave[7,n-1,:,ok_lon],'.-',label=list_exp[7],color='c')
plt.plot(lat, fp_wave[8,n-1,:,ok_lon],'.-',label=list_exp[8],color='deepskyblue')
plt.legend()
plt.title('Dominant Wave Frequency Longitudinal Transect',fontsize=14)
plt.ylabel('(1/s)',fontsize=14)
plt.xlabel('Lat',fontsize=14)
plt.grid(True)
plt.xlim([16,20])
file_name = folder_fig + 'fp_transet_zoom'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

