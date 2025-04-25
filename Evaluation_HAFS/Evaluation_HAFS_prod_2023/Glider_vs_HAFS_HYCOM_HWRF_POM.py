#%% User input

# forecasting cycle to be used

# Idalia
cycle = '2023082718'
storm_num = '10'
basin = 'al'
url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2023/usf-jaialai-20230806T1200.nc' 
# time storm is closest to glider
ystorm = 2023
mstorm = 8
dstorm = 29 
hstorm = 18

'''
# Franklin
cycle = '2023082212'
storm_num = '08'
basin = 'al'
url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2023/SG683-20230713T1239.nc' 
# time storm is closest to glider
ystorm = 2023
mstorm = 8
dstorm = 24 
hstorm = 6 
'''

exp_names = ['HFSA_oper']
exp_labels = ['HFSA']
exp_colors = ['darkviolet']

lon_lim = [-80,-60.0]
lat_lim = [10.0,30.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder1 = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch_folder2 = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch_folder1 + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

bath_file = scratch_folder1 +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'
#folder_glider_comp = '/home/Maria.Aristizabal/Repos/glider_model_comparisons_Python/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
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
                            figure_transect_time_vs_depth, glider_data_vector_to_array,\
                            grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

################################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

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
#%% Read glider data
gdata = xr.open_dataset(url_glider)#,decode_times=False)

dataset_id = gdata.id.split('_')[0]
temperature = np.asarray(gdata.variables['temperature'][:])
salinity = np.asarray(gdata.variables['salinity'][:])
density = np.asarray(gdata.variables['density'][:])
latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
depth = np.asarray(gdata.depth)

time = np.asarray(gdata.time)
time = np.asarray(mdates.num2date(mdates.date2num(time)))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
temperat =  temperature[oktimeg].T
salinit =  salinity[oktimeg].T
densit =  density[oktimeg].T
latitud = latitude[oktimeg]
longitud = longitude[oktimeg]
depthh = depth[oktimeg].T
timee = time[oktimeg]

contour_plot = 'no' # default value is 'yes'
delta_z = 0.2     # default value is 0.3

depthg, timeg, tempg, latg, long = glider_data_vector_to_array(depthh,timee,temperat,latitud,longitud)
depthg, timeg, saltg, latg, long = glider_data_vector_to_array(depthh,timee,salinit,latitud,longitud)

ok = np.where(np.isfinite(timeg[0,:]))[0]
timegg = timeg[0,ok]
tempgg = tempg[:,ok]
saltgg = saltg[:,ok]
depthgg = depthg[:,ok]
longg = long[0,ok]
latgg = latg[0,ok]

tempg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(tempgg,timegg,depthgg,delta_z)
saltg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(saltgg,timegg,depthgg,delta_z)

#tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
tstamp_glider = timeg_gridded

# Conversion from glider longitude and latitude to HYCOM convention
target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
lon_glider = target_lonG
lat_glider = target_latG

#################################################################################
#%% Loop the experiments
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan
rmw_track = np.empty((len(folder_exps),43))
rmw_track[:] = np.nan
target_temp_10m = np.empty((len(folder_exps),22))
target_temp_10m[:] = np.nan
target_salt_10m = np.empty((len(folder_exps),22))
target_salt_10m[:] = np.nan
target_time = np.empty((len(folder_exps),22))
target_time[:] = np.nan
target_OHC = np.empty((len(folder_exps),22))
target_OHC[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    

    if exp_labels[i] == 'HFSA' or exp_labels[i] == 'HFSB':
        #%% Get list files
        files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))

        #%% Reading HAFS/HYCOM grid
        hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
        lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
        lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
        depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

        #%% Read HAFS/HYCOM time
        time_hycom = []
        timestamp_hycom = []
        for n,file in enumerate(files_hafs_hycom):
            print(file)
            HYCOM = xr.open_dataset(file)
            t = HYCOM['MT'][:]
            timestamp = mdates.date2num(t)[0]
            time_hycom.append(mdates.num2date(timestamp))
            timestamp_hycom.append(timestamp)

        time_hycom = np.asarray(time_hycom)
        timestamp_hycom = np.asarray(timestamp_hycom)

        #%% Retrieve glider transect from HAFS_HYCOM
        ncfiles = files_hafs_hycom
        lon = lon_hafs_hycom
        lat = lat_hafs_hycom
        depth = depth_hafs_hycom

        if np.min(lon) < 0:
            lon_glid = longg
        else: 
            lon_glid = lon_glider
        lat_glid = lat_glider

        target_tt, target_temp, target_salt = \
        get_glider_transect_from_HAFS_HYCOM(ncfiles,lon,lat,depth,tstamp_glider,lon_glid,lat_glid)
        timestamp = mdates.date2num(target_tt)
        target_time[i,:] = timestamp
        depth = -depth
        target_tt = np.asarray(target_tt)

        # OHC
        for t in np.arange(len(target_tt)):
            temp = target_temp[:,t]
            salt = target_salt[:,t]
            dens = sw.dens(salt,temp,depth) 
            target_OHC[i,t] = OHC_from_profile(depth,temp,dens)

        # Vertical profile
        okt_model = np.where(mdates.date2num(target_tt) == mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)))[0][0]  
        plt.figure(figsize=(5,7)) 
        plt.plot(target_temp[:,okt_model-2:okt_model],depth,'o-r',alpha=0.2) 
        plt.plot(target_temp[:,okt_model:okt_model+2],depth,'o-b',alpha=0.2) 
        plt.plot(target_temp[:,okt_model-2:okt_model-1],depth,'o-r',alpha=0.2,label='12h before storm passage') 
        plt.plot(target_temp[:,okt_model:okt_model+1],depth,'o-b',alpha=0.2,label='12h after storm passage') 
        plt.ylim([-200,0])
        plt.xlim([17,31])
        plt.legend(loc='lower left')
        plt.title(exp_labels[i])
        plt.ylabel('Depth (m)',fontsize=14)
        plt.xlabel('Temperature ($^oC$)',fontsize=14)
 
        #%% Get storm track from trak atcf files
        if exp_labels[i] == 'HFSA':
            file_track = glob.glob(os.path.join(folder,'*hfsa.trak.atcfunix'))[0]
        if exp_labels[i] == 'HFSB':
            file_track = glob.glob(os.path.join(folder,'*hfsb.trak.atcfunix'))[0]
    
    if exp_labels[i] == 'HWRF':
        #%% Get list files
        files_pom = sorted(glob.glob(os.path.join(folder,'*pom*00*.nc')))    

        #%% Reading POM grid
        print('Retrieving coordinates from POM')
        grid_file = glob.glob(os.path.join(folder,'*pom.grid.nc'))[0]
        pom_grid = xr.open_dataset(grid_file)
        lon_pom = np.asarray(pom_grid['east_e'][:])
        lat_pom = np.asarray( pom_grid['north_e'][:])
        zlevc = np.asarray(pom_grid['zz'][:])
        topoz = np.asarray(pom_grid['h'][:])    
  
        #%% Read POM time
        time_pom = []
        for n,file in enumerate(files_pom):
            print(file)
            pom = xr.open_dataset(file)
            #ocean = xr.open_dataset(files_ocean[0],decode_times=False)
            t = pom.variables['time'][:]
            timestamp = mdates.date2num(t)[0]
            time_pom.append(mdates.num2date(timestamp))

        time_pom = np.asarray(time_pom)

        #%% Retrieve glider transect from HWRF_POM
        ncfiles = files_pom
        lon = lon_pom
        lat = lat_pom

        if np.min(lon) < 0:
            lon_glid = longg
        else:
            lon_glid = lon_glider
        lat_glid = lat_glider

        #target_t, target_temp, target_salt = \
        #get_glider_transect_from_pom(ncfiles,lon,lat,zlevc,topoz,tstamp_glider,lon_glid,lat_glid)
        #timestamp = mdates.date2num(target_t)

        target_temp = np.empty((len(zlevc),len(ncfiles)))
        target_temp[:] = np.nan
        target_salt = np.empty((len(zlevc),len(ncfiles)))
        target_salt[:] = np.nan
        target_depth = np.empty((len(zlevc),len(ncfiles)))
        target_depth[:] = np.nan
        target_topoz = np.empty((len(ncfiles)))
        target_topoz[:] = np.nan
        target_timestamp = np.empty((len(ncfiles)))
        target_timestamp[:] = np.nan

        for x,file in enumerate(files_pom):
            print(x)
            pom = xr.open_dataset(file)
            t = pom['time'][:]
            target_timestamp[x] = mdates.date2num(t)[0]

            # Interpolating latg and longlider into HYCOM grid
            sublon = np.interp(timestamp,tstamp_glider,lon_glid)
            sublat = np.interp(timestamp,tstamp_glider,lat_glid)
            oklon = int(np.round(np.interp(sublon,lon_pom[0,:],np.arange(len(lon_pom[0,:])))))
            oklat = int(np.round(np.interp(sublat,lat_pom[:,0],np.arange(len(lat_pom[:,0])))))

            target_temp[:,x] = np.asarray(pom['t'][0,:,oklat,oklon])
            target_salt[:,x] = np.asarray(pom['s'][0,:,oklat,oklon])
            target_topoz[x] = np.asarray(topoz[oklat,oklon])

        z_matrix_pom = np.dot(target_topoz.reshape(-1,1),zlevc.reshape(1,-1)).T
        depth = z_matrix_pom

        target_temp[target_temp == 0.0] = np.nan
        target_salt[target_salt == 0.0] = np.nan

        target_time[i,0:len(target_timestamp)] = target_timestamp
        target_tt = np.tile(target_timestamp,(len(zlevc),1))

        # OHC
        for t in np.arange(target_tt.shape[1]):
            temp = target_temp[:,t]
            salt = target_salt[:,t]
            dept = depth[:,t]
            dens = sw.dens(salt,temp,dept) 
            target_OHC[i,t] = OHC_from_profile(dept,temp,dens)

        #%% Get storm track from trak atcf files
        file_track = glob.glob(os.path.join(folder,'*trak.hwrf.atcfunix'))[0]

        # Vertical profile
        okt_model = np.where(target_tt[0,:] >= mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)))[0][0]  
        plt.figure(figsize=(5,7)) 
        plt.plot(target_temp[:,okt_model-2:okt_model],depth[:,okt_model-2:okt_model],'o-r',alpha=0.2)
        plt.plot(target_temp[:,okt_model:okt_model+2],depth[:,okt_model:okt_model+2],'o-b',alpha=0.2) 
        plt.plot(target_temp[:,okt_model-2:okt_model-1],depth[:,okt_model-2:okt_model-1],'o-r',alpha=0.2,label='12h before storm passage') 
        plt.plot(target_temp[:,okt_model:okt_model+1],depth[:,okt_model:okt_model+1],'o-b',alpha=0.2,label='12h after storm passage') 
        plt.ylim([-200,0])
        plt.xlim([17,31])
        plt.legend(loc='lower left')
        plt.title(exp_labels[i])
        plt.ylabel('Depth (m)',fontsize=14)
        plt.xlabel('Temperature ($^oC$)',fontsize=14)
        
    max_depth = 200
    kw_temp = dict(levels = np.arange(11,33,1))
    figure_transect_time_vs_depth(target_tt,depth,target_temp,date_ini,date_end,max_depth,kw_temp,'Spectral_r','($^oC$)')
    plt.contour(target_tt,depth,target_temp,levels=[26],colors = 'k')
    plt.title(exp_labels[i],fontsize=16)
    fig_name = 'temp_transect_'+exp_labels[i]+'_'+cycle+'.png'
    plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    max_depth = 200
    kw_salt = dict(levels = np.arange(33,37.5,0.3))
    figure_transect_time_vs_depth(target_tt,depth,target_salt,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
    plt.title(exp_labels[i],fontsize=16)
    fig_name = 'salt_transect_'+exp_labels[i]+'_'+cycle+'.png'
    plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    # Temp at 10 meters depth
    if depth.ndim == 2:
        for t in np.arange(depth.shape[1]):
            if np.nanmean(depth) < 0:
                okd = np.where(depth[:,t] <= -10)[0][0]
            else:
                okd = np.where(depth[:,t] >= 10)[0][0]
            target_temp_10m[i,t] = target_temp[okd,t]
            target_salt_10m[i,t] = target_salt[okd,t]
            #target_time[i,:] = timestamp 
            
    if depth.ndim == 1:
        if np.nanmean(depth) < 0:
            okd = np.where(depth <= -10)[0][0]
        else:
            okd = np.where(depth >= 10)[0][0]
        target_temp_10m[i,:] = target_temp[okd,:]
        target_salt_10m[i,:] = target_salt[okd,:]
    
    # Vertical profile
    '''
    okt_model = np.where(mdates.date2num(target_tt) == mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)))[0][0]  
    plt.figure(figsize=(5,7)) 
    plt.plot(target_temp[:,okt_model-2:okt_model],depth,'o-r',alpha=0.2) 
    plt.plot(target_temp[:,okt_model:okt_model+2],depth,'o-b',alpha=0.2) 
    plt.plot(target_temp[:,okt_model-2:okt_model-1],depth,'o-r',alpha=0.2,label='12h before storm passage') 
    plt.plot(target_temp[:,okt_model:okt_model+1],depth,'o-b',alpha=0.2,label='12h after storm passage') 
    plt.ylim([-200,0])
    plt.xlim([17,31])
    plt.legend()
    plt.title(exp_labels[i])
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Temperature ($^oC$)',fontsize=14)
    '''

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], rmw_track[i,0:okn] = get_storm_track_and_int(file_track,storm_num)

##################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
fig_name = 'track_glider_'+exp_labels[i]+'_'+cycle+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(long)-2,np.nanmax(long)+2])
plt.ylim([np.nanmin(latg)-2,np.nanmax(latg)+2])
fig_name = 'track_glider2_'+exp_labels[i]+'_'+cycle+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

##################################################################
#%% Figure time series temp at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
tempg10 = tempg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,tempg10,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)): 
    plt.plot(target_time[i,:],target_temp_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Temperature ($^oC$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Temperature at 10 m depth',fontsize=16)
plt.grid(True)
fig_name = 'temp_10m_'+cycle+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% Figure time series salt at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
saltg10 = saltg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,saltg10,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)):
    plt.plot(target_time[i,:],target_salt_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Salt',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Salinity at 10 m depth',fontsize=16)
plt.grid(True)
fig_name = 'salt_10m_'+cycle+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% Figures transects

#%% Glider
max_depth = 200
kw_temp = dict(levels = np.arange(11,33,1))
figure_transect_time_vs_depth(timegg,-depthg_gridded,tempg_gridded,date_ini,date_end,max_depth,kw_temp,'Spectral_r','($^oC$)')
plt.contour(timegg,-depthg_gridded,tempg_gridded,levels=[26],colors = 'k')
plt.title(dataset_id,fontsize=16)
fig_name = 'temp_transect_glider_'+dataset_id+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

max_depth = 200
kw_salt = dict(levels = np.arange(33,37.5,0.3))
figure_transect_time_vs_depth(timegg,-depthg_gridded,saltg_gridded,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
plt.title(dataset_id,fontsize=16)
fig_name = 'salt_transect_glider_'+dataset_id+'.png'
plt.savefig(fig_name,format='png',bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#################################################################################
#%% OHC time series

OHCg = np.empty((len(timegg)))
OHCg[:] = np.nan
for t in np.arange(len(timeg_gridded)):
    temp = tempg_gridded[:,t]
    salt = saltg_gridded[:,t]
    dept = depthg_gridded
    dens = sw.dens(salt,temp,dept) 
    OHCg[t] = OHC_from_profile(dept,temp,dens)

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,OHCg,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)):
    plt.plot(target_time[i,:],target_OHC[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('$(kJ/cm^2)$',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Ocean Heat Content',fontsize=16)
plt.grid(True)


#################################################################################
#%% Figures profiles before and after storm passage
okt_glider_aft = np.logical_and(timeg_gridded >= mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)),\
                 timeg_gridded <= mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)+timedelta(0.5))) 

okt_glider_bef = np.logical_and(timeg_gridded >= mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm)-timedelta(0.5)),\
                 timeg_gridded <= mdates.date2num(datetime(ystorm,mstorm,dstorm,hstorm))) 

plt.figure(figsize=(5,7)) 
plt.plot(tempg_gridded[:,okt_glider_bef],-depthg_gridded,'-r',alpha=0.2)
plt.plot(tempg_gridded[:,okt_glider_aft],-depthg_gridded,'-b',alpha=0.2)
plt.plot(tempg_gridded[:,okt_glider_bef][:,0],-depthg_gridded,'-r',alpha=0.2,label='12h before storm passage')
plt.plot(tempg_gridded[:,okt_glider_aft][:,0],-depthg_gridded,'-b',alpha=0.2,label='12h after storm passage')
plt.ylim([-200,0])
plt.xlim([17,31])
plt.legend(loc='lower left')
plt.title(dataset_id.split('-')[0])
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('Temperature ($^oC$)',fontsize=14)
plt.ylabel('Depth (m)')

