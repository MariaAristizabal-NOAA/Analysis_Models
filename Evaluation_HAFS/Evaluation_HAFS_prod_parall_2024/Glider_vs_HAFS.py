#%% User input

# Beryl
#cycle = '2024081312'

# Ernesto
#cycle = '2024081612'
#storm_num = '05'
#basin = 'al'
#storm_id = '05l'
#storm_name = 'Ernesto'

# Helena
cycle = '2024092500'
storm_num = '09'
basin = 'al'
storm_id = '09l'
storm_name = 'Helena'
# time eye passage closest to glider
teye = '2024092618' 

# Milton
'''
cycle = '2024100800'
storm_num = '14'
basin = 'al'
storm_id = '14l'
storm_name = 'Milton'
# time eye passage closest to glider
teye = '2024100918' 
'''

exp_names = ['HFSB_oper','hafs_20241220_v2p1b_baseline','hafs_20250306_v2p1b_hb43']
exp_labels = ['HFSB_oper','HAFSv2.1B baseline','HAFSv2.1B final']
exp_colors = ['lime','olive','darkorange']
hafs_ab = ['hfsb','hfsb','hfsb']
ocean = ['hycom','hycom','mom6']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']

'''
exp_names = ['HFSA_oper','hafs_20241220_v2p1a_baseline','hafs_20250210_v2p1a_ha30']
exp_labels = ['HFSA_oper','HAFSv2.1A baseline','HAFSv2.1A final']
exp_colors = ['purple','dodgerblue','#00c8c8']
hafs_ab = ['hfsa','hfsa','hfsa']
ocean = ['mom6','mom6','mom6']
scratch_folder = ['/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/','/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/']
'''

url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2024/usf-jaialai-20240920T0000.nc' 
#url_glider = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Gliders/2024/bios_jack-20240812T1810.nc' 

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_OCEAN,\
                            figure_transect_time_vs_depth,\
                            glider_data_vector_to_array,grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    folder_exps.append(scratch_folder[i] + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

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
delta_z = 1     # default value is 0.3

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

lon_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),len(time_fv3)))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),len(time_fv3)))
int_track[:] = np.nan
rmw_track = np.empty((len(folder_exps),len(time_fv3)))
rmw_track[:] = np.nan
target_temp_10m = np.empty((len(folder_exps),len(time_fv3)))
target_temp_10m[:] = np.nan
target_salt_10m = np.empty((len(folder_exps),len(time_fv3)))
target_salt_10m[:] = np.nan
target_time = np.empty((len(folder_exps),len(time_fv3)))
target_time[:] = np.nan
target_temp_hafs_ocean = np.empty((len(folder_exps),55,len(time_fv3)))
target_temp_hafs_ocean[:] = np.nan
target_salt_hafs_ocean = np.empty((len(folder_exps),55,len(time_fv3)))
target_salt_hafs_ocean[:] = np.nan
target_depth_hafs_ocean = np.empty((len(folder_exps),55))
target_depth_hafs_ocean[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)    
    #%% Get list files
    if ocean[i] == 'hycom':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon_hafs_ocean = np.asarray(hafs_ocean['Longitude'][:])
        lat_hafs_ocean = np.asarray(hafs_ocean['Latitude'][:])
        depth_hafs_ocean = np.asarray(hafs_ocean['Z'][:])
    if ocean[i] == 'mom6':
        files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
        hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
        lon_hafs_ocean = np.asarray(hafs_ocean['xh'][:])
        lat_hafs_ocean = np.asarray(hafs_ocean['yh'][:])
        depth_hafs_ocean = np.asarray(hafs_ocean['z_l'][:])

    #%% Get storm track from trak atcf files
    if hafs_ab[i] == 'hfsa':
        file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
    if hafs_ab[i] == 'hfsb':
        file_track = folder + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
    print(file_track)

    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], rmw_track[i,0:okn] = get_storm_track_and_int(file_track,storm_num)

    #%% Read time
    time_ocean = []
    timestamp_ocean = []
    for n,file in enumerate(files_hafs_ocean):
        print(file)
        hafs_ocean = xr.open_dataset(file)
        if ocean[i] == 'hycom':
            t = hafs_ocean['MT'][:]
        if ocean[i] == 'mom6':
            t = hafs_ocean['time'][:]
    timestamp = mdates.date2num(t)[0]
    time_ocean.append(mdates.num2date(timestamp))
    timestamp_ocean.append(timestamp)

    time_ocean = np.asarray(time_ocean)
    timestamp_ocean = np.asarray(timestamp_ocean)

    #%% Retrieve glider transect
    ncfiles = files_hafs_ocean
    lon = lon_hafs_ocean
    lat = lat_hafs_ocean
    depth = depth_hafs_ocean
    if ocean[i] == 'hycom':
        time_name = 'MT'
        temp_name = 'temperature'
        salt_name = 'salinity'
    if ocean[i] == 'mom6':
        time_name = 'time'
        temp_name = 'temp'
        salt_name = 'so'

    if np.min(lon) < 0:
        lon_glid = longg
    else: 
        lon_glid = lon_glider
    lat_glid = lat_glider

    target_t, target_temp_hafs_oc, target_salt_hafs_oc = \
    get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

    target_temp_hafs_ocean[i,0:target_temp_hafs_oc.shape[0],0:target_temp_hafs_oc.shape[1]] = target_temp_hafs_oc
    target_salt_hafs_ocean[i,0:target_temp_hafs_oc.shape[0],0:target_temp_hafs_oc.shape[1]] = target_salt_hafs_oc
    target_depth_hafs_ocean[i,0:len(depth_hafs_ocean)] = depth_hafs_ocean

    timestamp = mdates.date2num(target_t)
    
    max_depth = 200
    kw_temp = dict(levels = np.arange(18,33,1))
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_temp_hafs_oc,date_ini,date_end,max_depth,kw_temp,'Spectral_r','$C^o$')
    plt.title(exp_labels[i],fontsize=16)

    max_depth = 200
    levels = np.arange(35.6,37.1,0.1)
    kw_salt = dict(levels = levels)
    figure_transect_time_vs_depth(np.asarray(target_t),-depth,target_salt_hafs_oc,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
    plt.contour(np.asarray(target_t),-depth,target_salt_hafs_oc,levels)
    plt.title(exp_labels[i],fontsize=16)

    # Temp at 10 meters depth
    okd = np.where(depth_hafs_ocean <= 10)[0]
    target_temp_10m[i,0:len(target_t)] = target_temp_hafs_oc[okd[-1],:]
    target_salt_10m[i,0:len(target_t)] = target_salt_hafs_oc[okd[-1],:]
    target_time[i,0:len(target_t)] = timestamp 
  
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
plt.plot(long[0,:], latg[0,:],'.-',color='orange',label='Glider Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')

###################################################################%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.plot(long[0,:], latg[0,:],'.-',color='orange',label='Glider Track')
plt.legend(loc='lower right',bbox_to_anchor=[1.2,0])
#plt.legend(loc='lower right',bbox_to_anchor=[1.3,0.8])
#plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(long)-4,np.nanmax(long)+4])
plt.ylim([np.nanmin(latg)-4,np.nanmax(latg)+4])

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

#################################################################################
#%% Figure time series salt at 10 m depth

okd = np.where(depthg_gridded <= 10)[0]
saltg10 = saltg_gridded[okd[-1],:]
#,'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k'
fig, ax = plt.subplots(figsize=(8,3))
plt.plot(timegg,saltg10,'o-',color='k',label=dataset_id.split('-')[0],markeredgecolor='k')
for i in np.arange(len(exp_names)):
    plt.plot(target_time[i,:],target_salt_10m[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
#plt.legend()
plt.legend(loc='upper right',bbox_to_anchor=[1.1,1.3])
plt.ylabel('Salt',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.title('Salinity at 10 m depth',fontsize=16)
plt.grid(True)

#################################################################################
#%% Figures transects

#%% Glider
max_depth = 200
date_inikw_temp = dict(levels = np.arange(18,32,1))
figure_transect_time_vs_depth(timegg,-depthg_gridded,tempg_gridded,date_ini,date_end,max_depth,kw_temp,'Spectral_r','$C^o$')
plt.title(dataset_id,fontsize=16)

max_depth = 200
kw_salt = dict(levels = np.arange(35.6,37+0.1,0.1))
figure_transect_time_vs_depth(timegg,-depthg_gridded,saltg_gridded,date_ini,date_end,max_depth,kw_salt,'YlGnBu_r',' ')
plt.title(dataset_id,fontsize=16)

# alternative way to produce contour figure
#Time window
year_ini = int(date_ini.split('/')[0])
month_ini = int(date_ini.split('/')[1])
day_ini = int(date_ini.split('/')[2])
year_end = int(date_end.split('/')[0])
month_end = int(date_end.split('/')[1])
day_end = int(date_end.split('/')[2])
tini = datetime(year_ini, month_ini, day_ini)
tend = datetime(year_end, month_end, day_end)

okt = np.isfinite(temperat)
tempera_colors = temperat[okt]
timee_colors = timee[okt]
depthh_colors = depthh[okt]

min_val = 18
max_val = 31
dt = 1
levels = np.arange(min_val,max_val+dt,dt)
lev_norm = (levels-levels[0])/(levels[-1]-levels[0])
tempera_colors_norm = (tempera_colors - min_val)/(max_val-min_val)

color_map = 'Spectral_r'
cmap_mod = plt.get_cmap(color_map,len(levels))
new_cmap = ListedColormap(cmap_mod(np.linspace(0, 1, len(levels))))
colors = new_cmap(np.arange(len(levels)))

#plt.figure()
#for x in np.arange(len(levels)):
#    plt.plot(x,x,marker='o',color=colors[x,:])

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

# Salinity
oks = np.isfinite(salinit)
salinity_colors = salinit[oks]
timee_colors = timee[oks]
depthh_colors = depthh[oks]

#kw_salt = dict(levels = np.arange(34.9,36.7,0.1))
min_val = 35.6
max_val = 37
dt = 0.1
levels = np.arange(min_val,max_val,dt)
lev_norm = (levels-levels[0])/(levels[-1]-levels[0])
salinity_colors_norm = (salinity_colors - min_val)/(max_val-min_val)

color_map = 'YlGnBu_r'
cmap_mod = plt.get_cmap(color_map,len(levels))
new_cmap = ListedColormap(cmap_mod(np.linspace(0, 1, len(levels))))
colors = new_cmap(np.arange(len(levels)))

colorss = np.empty((salinity_colors.shape[0],colors.shape[1]))
colorss[:] = [1,1,1,1]
for i,temp in enumerate(salinity_colors_norm):
    if temp < 0:
        colorss[i,:] = colors[0,:]
    else:
        okp = np.where(lev_norm <= temp)[0][-1]
        colorss[i,:] = colors[okp,:]

levels_cont = np.arange(min_val,max_val+dt,dt)
kw = dict(levels = levels_cont)
fig, ax = plt.subplots(figsize=(8, 4))
plt.scatter(timee_colors,-depthh_colors,marker='o',s=20,color=colorss)
cs = plt.contourf(timegg,-depthg_gridded,saltg_gridded,cmap=new_cmap,**kw)
cbar = plt.colorbar(cs)
ax.set_ylabel('Depth (m)',fontsize=14)
cbar.ax.set_ylabel(' ',fontsize=14)
xvec = [tini + timedelta(int(dt)) for dt in np.arange((tend-tini).days+1)[::2]]
plt.xticks(xvec,fontsize=12)
xfmt = mdates.DateFormatter('%b-%d')
ax.xaxis.set_major_formatter(xfmt)
plt.ylim(-np.abs(max_depth),0)
plt.xlim(tini,tend)
plt.title(dataset_id,fontsize=16)

# Figures vertical profiles
y = int(teye[0:4])
m = int(teye[4:6])
d = int(teye[6:8])
h = int(teye[8:10])
time_eye_passage = mdates.date2num(datetime(y,m,d,h,0,0)) 
before_eye = timeg[0,:] < time_eye_passage 
after_eye = timeg[0,:] >= time_eye_passage 

fig, ax = plt.subplots(figsize=(8,11))
plt.plot(tempg[0,before_eye],-depthg[0,before_eye],'.',color='r',label='Before')
plt.plot(tempg[0,after_eye],-depthg[0,after_eye],'.',color='b',label='After')
plt.plot(tempg[:,before_eye],-depthg[:,before_eye],'.',color='r',alpha=0.5,markersize=1)
plt.plot(tempg[:,after_eye],-depthg[:,after_eye],'.',color='b',alpha=0.3,markersize=1)
plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Temperature ($^oC$)',fontsize=16)
plt.legend()

# Initial conditions profiles
okt_glider = np.logical_and(timeg[0,:] >= mdates.date2num(tini + timedelta(hours=0)),timeg[0,:]<mdates.date2num(tini + timedelta(hours=3)))
fig, ax = plt.subplots(figsize=(8,11))
plt.plot(tempg[0,okt_glider],-depthg[0,okt_glider],'-',color='black',label=dataset_id.split('-')[0])
plt.plot(tempg[:,okt_glider],-depthg[:,okt_glider],'-',color='black',linewidth=3)
for i,ff in enumerate(folder_exps): 
    okt_model = target_time[i,:] < mdates.date2num(tini + timedelta(hours=3)) 
    plt.plot(target_temp_hafs_ocean[i,:,okt_model][0,:],-target_depth_hafs_ocean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')

plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Temperature ($^oC$)',fontsize=16)
plt.title(str(tini + timedelta(hours=0))+ '-' + str(tini + timedelta(hours=3)),fontsize=16)
plt.legend()
plt.ylim([-160,0])
plt.xlim([16,31])

# Initial conditions profiles
okt_glider = np.logical_and(timeg[0,:] >= mdates.date2num(tini + timedelta(hours=0)),timeg[0,:]<mdates.date2num(tini + timedelta(hours=3)))
fig, ax = plt.subplots(figsize=(8,11))
plt.plot(saltg[0,okt_glider],-depthg[0,okt_glider],'-',color='black',label=dataset_id.split('-')[0])
plt.plot(saltg[:,okt_glider],-depthg[:,okt_glider],'-',color='black')
for i,ff in enumerate(folder_exps):
    okt_model = target_time[i,:] < mdates.date2num(tini + timedelta(hours=3))
    plt.plot(target_salt_hafs_ocean[i,:,okt_model][0,:],-target_depth_hafs_ocean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')

plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Salinity ',fontsize=16)
plt.title(str(tini + timedelta(hours=0))+ '-' + str(tini + timedelta(hours=3)),fontsize=16)
plt.legend()
plt.ylim([-140,0])
#plt.xlim([19,31])


# Subsequent time profiles
dt_ini = 42
dt_end = 45
okt_glider = np.logical_and(timeg[0,:] >= mdates.date2num(tini + timedelta(hours=dt_ini)),timeg[0,:]<mdates.date2num(tini + timedelta(hours=dt_end)))

fig, ax = plt.subplots(figsize=(8,11))
plt.plot(tempg[0,okt_glider],-depthg[0,okt_glider],'-',color='black',label=dataset_id.split('-')[0])
plt.plot(tempg[:,okt_glider],-depthg[:,okt_glider],'-',color='black',linewidth=3)
for i,ff in enumerate(folder_exps): 
    okt_model = np.logical_and(target_time[i,:] >= mdates.date2num(tini + timedelta(hours=dt_ini)),target_time[i,:]<mdates.date2num(tini + timedelta(hours=dt_end)))
    plt.plot(target_temp_hafs_ocean[i,:,okt_model][0,:],-target_depth_hafs_ocean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')

plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Temperature ($^oC$)',fontsize=16)
plt.title(str(tini + timedelta(hours=dt_ini))+ '-' + str(tini + timedelta(hours=dt_end)),fontsize=16)
plt.legend()
plt.ylim([-160,0])
plt.xlim([16,31])

fig, ax = plt.subplots(figsize=(8,11))
plt.plot(saltg[0,okt_glider],-depthg[0,okt_glider],'-',color='black',label=dataset_id.split('-')[0])
plt.plot(saltg[:,okt_glider],-depthg[:,okt_glider],'-',color='black')
for i,ff in enumerate(folder_exps):
    okt_model = np.logical_and(target_time[i,:] >= mdates.date2num(tini + timedelta(hours=dt_ini)),target_time[i,:]<mdates.date2num(tini + timedelta(hours=dt_end)))
    plt.plot(target_salt_hafs_ocean[i,:,okt_model][0,:],-target_depth_hafs_ocean[i,:],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')

plt.ylabel('Depth (m)',fontsize=16)
plt.xlabel('Salinity',fontsize=16)
plt.title(str(tini + timedelta(hours=dt_ini))+ '-' + str(tini + timedelta(hours=dt_end)),fontsize=16)
plt.legend()
plt.ylim([-140,0])
#plt.xlim([19,31])

