#%% User input

# forecasting cycle to be used
# Polo 17e
cycle = '2026092100'
storm_num = '17'
basin = 'ep'
storm_id = '17e'
storm_name= 'polo'
models = ['hfsa']

exp_names = ['HFSA_oper']
exp_labels = ['HFSA']
exp_colors = ['darkviolet']

lon_lim = [-115,-95]
lat_lim = [5,25]

scratch_folder = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/'
abdeck_folder = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/abdeck/'
folder_myutils= '/ncrc/home1/Maria.Aristizabal/Maria_Utils/'

cartopyDataDir = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/local/share/cartopy/'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

#url_glider = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/Data/Gliders/2026/sg623-20260707T0000.nc'
url_glider = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/Data/Gliders/2026/sg625-20260630T0000.nc'

################################################################################
import sys
import glob
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_glider_transect_from_HAFS_OCEAN,\
                            figure_transect_time_vs_depth, glider_data_vector_to_array,\
                            grid_glider_data

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder + exp_names[i] + '/' + cycle + '/' + storm_id + '/')

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
#target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
#lon_glider = target_lonG
#lat_glider = target_latG
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
target_temp_10m = np.empty((len(folder_exps),43))
target_temp_10m[:] = np.nan
target_salt_10m = np.empty((len(folder_exps),43))
target_salt_10m[:] = np.nan
target_time = np.empty((len(folder_exps),43))
target_time[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    #%% Get list files
    files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*storm.atm*.grb2')))

    hafs_ocean = xr.open_dataset(files_hafs_ocean[0])
    lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
    lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
    depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])

    #%% Get storm track from trak atcf files
    file_track = folder+storm_id+'.'+cycle+'.'+models[i]+'.trak.atcfunix'

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], rmw_track[i,0:okn] = get_storm_track_and_int(file_track,storm_num)

    #%% Read HAFS/Fv3 time
    '''
    time_fv3 = []
    for n,file in enumerate(files_hafs_fv3):
        print(file)
        #FV3 = xr.open_dataset(file)
        #t = FV3.variables['time'][:]
        #timestamp = mdates.date2num(t)[0]
        #time_fv3.append(mdates.num2date(timestamp))

        FV3 = xr.open_dataset(file,engine="pynio")
        t0 = FV3.variables['TMP_P0_L1_GLL0'].attrs['initial_time']
        dt = FV3.variables['TMP_P0_L1_GLL0'].attrs['forecast_time'][0]
        time_fv3.append(datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt)))

    time_fv3 = np.asarray(time_fv3)
    '''

    #%% Read HAFS/MOM6 time
    time_ocean = []
    timestamp_ocean = []
    for n,file in enumerate(files_hafs_ocean):
        print(file)
        OCEAN = xr.open_dataset(file)
        t = OCEAN['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_ocean.append(mdates.num2date(timestamp))
        timestamp_ocean.append(timestamp)

    time_ocean = np.asarray(time_ocean)
    timestamp_ocean = np.asarray(timestamp_ocean)

    #%% Retrieve glider transect from HAFS_HYCOM
    ncfiles = files_hafs_ocean
    lon = lon_hafs_ocean
    lat = lat_hafs_ocean
    depth = depth_hafs_ocean
    time_name = 'time'
    temp_name = 'temp'
    salt_name = 'so'

    if np.min(lon) < 0:
        lon_glid = longg
    else: 
        lon_glid = lon_glider
    lat_glid = latgg

    target_t, target_temp_hafs_ocean, target_salt_hafs_ocean = \
    get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

    timestamp = mdates.date2num(target_t)
    
    max_depth = 200
    kw_temp = dict(levels = np.arange(11,33,1))
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_temp_hafs_ocean,date_ini,date_end,max_depth,kw_temp,'Spectral_r','Degrees')
    plt.title(exp_labels[i],fontsize=16)

    max_depth = 200
    kw_salt = dict(levels = np.arange(33,35.8,0.2))
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_salt_hafs_ocean,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
    plt.title(exp_labels[i],fontsize=16)

    # Temp at 10 meters depth
    okd = np.where(depth_hafs_ocean <= 10)[0]
    target_temp_10m[i,:] = target_temp_hafs_ocean[okd[-1],:]
    target_salt_10m[i,:] = target_salt_hafs_ocean[okd[-1],:]
    target_time[i,:] = timestamp 
  
##################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')

for i in np.arange(len(exp_names)): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')

plt.xlim(lon_lim)
plt.ylim(lat_lim)

# Add gridlines and labels
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='blue',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(long)-1,np.nanmax(long)+1])
plt.ylim([np.nanmin(latg)-1,np.nanmax(latg)+1])

# Add gridlines and labels
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-180., 180.+1, 2))
gl.ylocator = mticker.FixedLocator(np.arange(-90., 90.+1, 2))
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

##################################################################
#%% Figure time series temp at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
tempg10 = tempg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,tempg10,'o-',color='blue',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)): 
    plt.plot(target_time[i,:],target_temp_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Temperature ($^oC$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Temperature at 10 m depth',fontsize=16)
plt.grid(True)

#################################################################################
#%% Figure time series salt at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
saltg10 = saltg_gridded[okd[-1],:]

fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,saltg10,'o-',color='blue',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)):
    plt.plot(target_time[i,:],target_salt_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
plt.legend()
plt.ylabel('Salt',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Salinity at 10 m depth',fontsize=16)
plt.grid(True)

#################################################################################
#%% Figures transects

#%% Glider
max_depth = 200
kw_temp = dict(levels = np.arange(11,34,1))
figure_transect_time_vs_depth(timegg,-depthg_gridded,tempg_gridded,date_ini,date_end,max_depth,kw_temp,'Spectral_r','Degress')
plt.title(dataset_id,fontsize=16)

max_depth = 200
kw_salt = dict(levels = np.arange(33,35.8,0.2))
figure_transect_time_vs_depth(timegg,-depthg_gridded,saltg_gridded,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
plt.title(dataset_id,fontsize=16)

####################################################
# alternative way to produce contour figure
#Time window
'''
year_ini = int(date_ini.split('/')[0])
month_ini = int(date_ini.split('/')[1])
day_ini = int(date_ini.split('/')[2])
year_end = int(date_end.split('/')[0])
month_end = int(date_end.split('/')[1])
day_end = int(date_end.split('/')[2])
tini = datetime(year_ini, month_ini, day_ini)
tend = datetime(year_end, month_end, day_end)
'''

okt = np.isfinite(temperat)
tempera_colors = temperat[okt]
timee_colors = timee[okt]
depthh_colors = depthh[okt]

min_val = 11
max_val = 33
dt = 1
levels = np.arange(min_val,max_val+dt,dt)
lev_norm = (levels-levels[0])/(levels[-1]-levels[0])
tempera_colors_norm = (tempera_colors - min_val)/(max_val-min_val)

color_map = 'Spectral_r'
cmap_mod = plt.get_cmap(color_map,len(levels))
new_cmap = ListedColormap(cmap_mod(np.linspace(0, 1, len(levels))))
colors = new_cmap(np.arange(len(levels)))

colorss = np.empty((tempera_colors.shape[0],colors.shape[1]))
colorss[:] = [1,1,1,1]
for i,temp in enumerate(tempera_colors_norm):
    if temp < 0:
        colorss[i,:] = colors[0,:]
    else:
        okp = np.where(lev_norm <= temp)[0][-1]
        colorss[i,:] = colors[okp,:]

levels_cont = np.arange(min_val,max_val+1+dt,dt)
kw = dict(levels = levels_cont)
fig, ax = plt.subplots(figsize=(8, 4))
plt.scatter(timee_colors,-depthh_colors,marker='o',s=20,color=colorss)
cs = plt.contourf(timegg,-depthg_gridded,tempg_gridded,cmap=new_cmap,**kw)
cbar = plt.colorbar(cs)
ax.set_ylabel('Depth (m)',fontsize=14)
cbar.ax.set_ylabel('$C^o$',fontsize=14)
xvec = [tini + timedelta(int(dt)) for dt in np.arange((tend-tini).days+1)[::2]]
plt.xticks(xvec,fontsize=12)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.ylim(-np.abs(max_depth),0)
plt.xlim(tini,tend)
plt.title(dataset_id,fontsize=16)
