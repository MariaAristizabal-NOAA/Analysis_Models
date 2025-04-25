#%% User input
# forecasting cycle to be used

# User input
# Encompases Helene and Milton 2024
yyyymmddhh_ini = '2024092312'
# number of days for tha aanlysis from yyyymmddhh_ini
delta_days = 48 

file_drifter = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Scripts_lagrang_drifters/LDL_AtlanticHurricane2024.nc'
folder_SOHCS = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NESDIS_OHC/2024/'
folder_CoralTemp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Coral_Reef_Watch/2024/'
folder_METOP = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/METOP_SST/2024/'
folder_GeoPolarBlend = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Geo_Polar_Blend_SST/2024/'

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

############################################################
def get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output):
    subvar_obs = []
    subvar_model = []

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]
    var_obs = var_obss[oklo]

    okt_sort = np.argsort(timestamp_obs)
    timestamp_obs_sort = timestamp_obs[okt_sort]
    lon_obs_sort = lon_obs[okt_sort]
    lat_obs_sort = lat_obs[okt_sort]
    var_obs_sort = var_obs[okt_sort]

    target_time = []

    for x,file in enumerate(files_model):
        print(x,' ',file)
        model = xr.open_dataset(file)
        if file.split('.')[-1] == 'nc':
            t =  mdates.num2date(mdates.date2num(np.asarray(model[time_name][:])))[0]
            ti = t - timedelta(hours=int(freq_output/2))
            tii = mdates.date2num(ti)
            te = t + timedelta(hours=int(freq_output/2))
            tee = mdates.date2num(te)
            timestamp = mdates.date2num(t)
            target_time.append(mdates.num2date(timestamp))
            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])
        if file.split('.')[-1] == 'grb2':
            t0 = model[var_name].attrs['initial_time']
            dt = model[var_name].attrs['forecast_time'][0]
            t = datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt))
            timestamp = mdates.date2num(t)
            target_time.append(t)
            lon_model = np.asarray(model.lon_0)
            lat_model = np.asarray(model.lat_0)

        okt = np.logical_and(timestamp_obs_sort > tii,timestamp_obs_sort < tee) 
        sublon_obs = lon_obs_sort[okt]
        sublat_obs = lat_obs_sort[okt]

        if len(sublon_obs) != 0 or len(sublat_obs) != 0:
            okposs = []
            if lon_model.ndim == 1 and lat_model.ndim == 1:
                oklonn = np.round(np.interp(sublon_obs,lon_model,np.arange(len(lon_model)))).astype(int)
                oklatt = np.round(np.interp(sublat_obs,lat_model,np.arange(len(lat_model)))).astype(int)
                for i in np.arange(len(oklonn)):
                    okposs.append([oklatt[i],oklonn[i]])
            if lon_model.ndim == 2 and lat_model.ndim == 2:
                for i in np.arange(len(sublon_obs)):
                    oklatt, oklonn = find_grid_position_hycom(lat_model,lon_model,sublat_obs[i],sublon_obs[i])
                    okposs.append([oklatt,oklonn])

            okpos,okpos_ind = np.unique(okposs,return_index=True,axis=0)

            subvar_obs.append(var_obs_sort[okt][okpos_ind].tolist())

            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name])[0,depth_level,okpos[:,0],okpos[:,1]]
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name])[0,okpos[:,0],okpos[:,1]]
            var[var==0] = np.nan
            subvar_model.append(var.tolist())

    return subvar_obs,subvar_model

#####################################################################
def find_grid_position_hycom(lat,lon,target_lat,target_lon):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklat,oklon

#####################################################################
def taylor_template(angle_lim,std_lim):

    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure(figsize=(9,7))
    tr = PolarAxes.PolarTransform()

    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))

    STDgrid=np.arange(0,std_lim,.5)
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))

    ra0, ra1 =0, angle_lim
    cz0, cz1 = 0, std_lim
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)

    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["left"].label.set_text("Normalized Standard Deviation")
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")

    ax1.axis["bottom"].set_visible(False)
    ax1 = ax1.get_aux_axes(tr)

    plt.grid(linestyle=':',alpha=0.5)

    return fig,ax1

#####################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            get_var_from_model_following_trajectory


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#####################################################################
#%% Time window
date_ini = yyyymmddhh_ini[0:4]+'/'+yyyymmddhh_ini[4:6]+'/'+yyyymmddhh_ini[6:8]+'/'+yyyymmddhh_ini[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(days=delta_days)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#####################################################################
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

#####################################################################
#%% Read drifter data

url = file_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
platform_code = np.asarray(gdata.platform_code)
wind_speed = np.asarray(gdata.windspd)
sea_level_press = np.asarray(gdata.slp)
sea_surface_temp = np.asarray(gdata.sst)
sea_surface_salt = np.asarray(gdata.salinity)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(times)
times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(times) >= mdates.date2num(tini),\
                         mdates.date2num(times) <= mdates.date2num(tend))

# Fields within time window
timeDr = times[oktimeg]
timestampDr = timestamps[oktimeg]
latDr = latitude[oktimeg]
lonDr = longitude[oktimeg]
platform_codeDr = platform_code[oktimeg]
windspdDr = wind_speed[oktimeg]
slpDr = sea_level_press[oktimeg]
sstDr = sea_surface_temp[oktimeg]
sssDr = sea_surface_salt[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
platform_codeD = platform_codeDr[oklat][oklon]
windspdD = windspdDr[oklat][oklon]
slpD = slpDr[oklat][oklon]
sstD = sstDr[oklat][oklon]
sssD = sssDr[oklat][oklon]

codes = np.unique(platform_codeD)

number_obs = lonD.shape[0]

# Figure track
lev = np.arange(-9000,9100,100)
fig,ax = plt.subplots(figsize=(9,4))
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lonD, latD,'.',color='orangered',label='Drifters')
plt.legend(loc='upper right',bbox_to_anchor=[1.45,1.0])
#plt.legend(loc='upper right')
plt.title(date_ini[0:13]+'-'+date_end[0:13])
plt.axis('scaled')
plt.xlim([np.nanmin(lonD)-1,np.nanmax(lonD)+1])
plt.ylim([np.nanmin(latD)-1,np.nanmax(latD)+1])

#####################################################################
#%% Read NESDIS SST

files_sohcs = sorted(glob.glob(os.path.join(folder_SOHCS,'ohc_na*2024*.nc')))
SOHCS = xr.open_dataset(files_sohcs[0],decode_times=False)
lon = np.asarray(SOHCS['longitude'][:])

var_name = 'sst'
time_name = 'time'
lon_name = 'longitude'
lat_name = 'latitude'

files_model = files_sohcs
depth_level = 0
freq_output_hours = 24

timestamp_obss = timestampD
lon_obss = lonD
lat_obss = latD
var_obss = sstD

subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output_hours)

############################################################
subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])
ok1 = np.isfinite(subsst_obs_arr)
subsst_obs_arra = subsst_obs_arr[ok1]
subsst_model_arra = subsst_model_arr[ok1]
ok2 = np.isfinite(subsst_model_arra)
subsst_obs_array = subsst_obs_arra[ok2]
subsst_model_array = subsst_model_arra[ok2]
number_obs_sohcs = lonD.shape[0]
number_obs_used_sohcs = len(subsst_obs_array)
std_subsst_obs_sohcs = np.nanstd(subsst_obs_array)
std_subsst_model_sohcs = np.nanstd(subsst_model_array)
bias_sohcs = np.nanmean(subsst_obs_array) - np.nanmean(subsst_model_array)
corr_sohcs = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]

plt.figure()
plt.plot(subsst_obs_array,subsst_model_array,'.',color='darkorange',markersize=7,markeredgecolor='k')
plt.plot(np.arange(33),np.arange(33),'-',color='silver',linewidth=2)
plt.ylim([17,33])
plt.xlim([17,33])
plt.xlabel('Drifters Observations',fontsize=14)
plt.ylabel('SOHCS',fontsize=14)
plt.title('SST',fontsize=18)
plt.text(26,23,'Total Observations = ' + str(number_obs_sohcs))
plt.text(26,22,'Observations used = ' + str(number_obs_used_sohcs))
plt.text(26,21,'Bias = ' + str(np.round(bias_sohcs,2)))
plt.text(26,20,'STD obs = ' + str(np.round(std_subsst_obs_sohcs,2)))
plt.text(26,19,'STD model = ' + str(np.round(std_subsst_model_sohcs,2)))
plt.text(26,18,'Corr = ' + str(np.round(corr_sohcs,2)))

#####################################################################
#%% Read Coral Reef Watch SST

files_coraltemp = sorted(glob.glob(os.path.join(folder_CoralTemp,'coraltemp_v3.1_*.nc')))
CoralTemp = xr.open_dataset(files_coraltemp[0],decode_times=False)
lon = np.asarray(CoralTemp['lon'][:])

var_name = 'analysed_sst'
time_name = 'time'
lon_name = 'lon'
lat_name = 'lat'

files_model = files_coraltemp
depth_level = 0
freq_output_hours = 24

timestamp_obss = timestampD
lon_obss = lonD
lat_obss = latD
var_obss = sstD

subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output_hours)

############################################################
subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])
ok1 = np.isfinite(subsst_obs_arr)
subsst_obs_arra = subsst_obs_arr[ok1]
subsst_model_arra = subsst_model_arr[ok1]
ok2 = np.isfinite(subsst_model_arra)
subsst_obs_array = subsst_obs_arra[ok2]
subsst_model_array = subsst_model_arra[ok2]
number_obs_coraltemp = lonD.shape[0]
number_obs_used_coraltemp = len(subsst_obs_array)
std_subsst_obs_coraltemp = np.nanstd(subsst_obs_array)
std_subsst_model_coraltemp = np.nanstd(subsst_model_array)
bias_coraltemp = np.nanmean(subsst_obs_array) - np.nanmean(subsst_model_array)
corr_coraltemp = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]

plt.figure()
plt.plot(subsst_obs_array,subsst_model_array,'.',color='lime',markersize=7,markeredgecolor='k')
plt.plot(np.arange(33),np.arange(33),'-',color='silver',linewidth=2)
plt.ylim([17,33])
plt.xlim([17,33])
plt.xlabel('Drifters Observations',fontsize=14)
plt.ylabel('CoralTemp',fontsize=14)
plt.title('SST',fontsize=18)
plt.text(26,23,'Total Observations = ' + str(number_obs_coraltemp))
plt.text(26,22,'Observations used = ' + str(number_obs_used_coraltemp))
plt.text(26,21,'Bias = ' + str(np.round(bias_coraltemp,2)))
plt.text(26,20,'STD obs = ' + str(np.round(std_subsst_obs_coraltemp,2)))
plt.text(26,19,'STD model = ' + str(np.round(std_subsst_model_coraltemp,2)))
plt.text(26,18,'Corr = ' + str(np.round(corr_coraltemp,2)))

######################################################################
#%% Read METOP SST

files_metop = sorted(glob.glob(os.path.join(folder_METOP,'*OSISAF-L3C_GHRSST-SSTsubskin-AVHRR_SST_METOP_B_GLB-sstglb_metop01*.nc')))
METOP = xr.open_dataset(files_metop[0],decode_times=False)
lon = np.asarray(METOP['lon'][:])

var_name = 'sea_surface_temperature'
time_name = 'time'
lon_name = 'lon'
lat_name = 'lat'

files_model = files_metop
depth_level = 0
freq_output_hours = 24

timestamp_obss = timestampD
lon_obss = lonD
lat_obss = latD
var_obss = sstD

subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output_hours)

############################################################
subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])
ok1 = np.isfinite(subsst_obs_arr)
subsst_obs_arra = subsst_obs_arr[ok1]
subsst_model_arra = subsst_model_arr[ok1]
ok2 = np.isfinite(subsst_model_arra)
subsst_obs_array = subsst_obs_arra[ok2]
subsst_model_array = subsst_model_arra[ok2] - 272.15
number_obs_metop = lonD.shape[0]
number_obs_used_metop = len(subsst_obs_array)
std_subsst_obs_metop = np.nanstd(subsst_obs_array)
std_subsst_model_metop = np.nanstd(subsst_model_array)
bias_metop = np.nanmean(subsst_obs_array) - np.nanmean(subsst_model_array)
corr_metop = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]

plt.figure()
plt.plot(subsst_obs_array,subsst_model_array,'.',color='cyan',markersize=7,markeredgecolor='k')
plt.plot(np.arange(33),np.arange(33),'-',color='silver',linewidth=2)
plt.ylim([17,33])
plt.xlim([17,33])
plt.xlabel('Drifters Observations',fontsize=14)
plt.ylabel('METOP',fontsize=14)
plt.title('SST',fontsize=18)
plt.text(26,23,'Total Observations = ' + str(number_obs_metop))
plt.text(26,22,'Observations used = ' + str(number_obs_used_metop))
plt.text(26,21,'Bias = ' + str(np.round(bias_metop,2)))
plt.text(26,20,'STD obs = ' + str(np.round(std_subsst_obs_metop,2)))
plt.text(26,19,'STD model = ' + str(np.round(std_subsst_model_metop,2)))
plt.text(26,18,'Corr = ' + str(np.round(corr_metop,2)))

#####################################################################
#%% Read OSPO SST

files_geopolarblend = sorted(glob.glob(os.path.join(folder_GeoPolarBlend,'*OSPO-L4_GHRSST-SSTfnd-Geo_Polar_Blended-GLOB-v02.0-fv01.0*.nc')))
GeoPolarBlend = xr.open_dataset(files_geopolarblend[0],decode_times=False)
lon = np.asarray(GeoPolarBlend['lon'][:])

var_name = 'analysed_sst'
time_name = 'time'
lon_name = 'lon'
lat_name = 'lat'

files_model = files_geopolarblend
depth_level = 0
freq_output_hours = 24

timestamp_obss = timestampD
lon_obss = lonD
lat_obss = latD
var_obss = sstD

subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_name,lon_name,lat_name,var_name,depth_level,freq_output_hours)

############################################################
subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])
ok1 = np.isfinite(subsst_obs_arr)
subsst_obs_arra = subsst_obs_arr[ok1]
subsst_model_arra = subsst_model_arr[ok1]
ok2 = np.isfinite(subsst_model_arra)
subsst_obs_array = subsst_obs_arra[ok2]
subsst_model_array = subsst_model_arra[ok2] - 272.15
number_obs_geopolarblend = lonD.shape[0]
number_obs_used_geopolarblend = len(subsst_obs_array)
std_subsst_obs_geopolarblend = np.nanstd(subsst_obs_array)
std_subsst_model_geopolarblend = np.nanstd(subsst_model_array)
bias_geopolarblend = np.nanmean(subsst_obs_array) - np.nanmean(subsst_model_array)
corr_geopolarblend = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]

plt.figure()
plt.plot(subsst_obs_array,subsst_model_array,'.',color='magenta',markersize=7,markeredgecolor='k')
plt.plot(np.arange(33),np.arange(33),'-',color='silver',linewidth=2)
plt.ylim([17,33])
plt.xlim([17,33])
plt.xlabel('Drifters Observations',fontsize=14)
plt.ylabel('GeoPolarBlend',fontsize=14)
plt.title('SST',fontsize=18)
plt.text(26,23,'Total Observations = ' + str(number_obs_geopolarblend))
plt.text(26,22,'Observations used = ' + str(number_obs_used_geopolarblend))
plt.text(26,21,'Bias = ' + str(np.round(bias_geopolarblend,2)))
plt.text(26,20,'STD obs = ' + str(np.round(std_subsst_obs_geopolarblend,2)))
plt.text(26,19,'STD model = ' + str(np.round(std_subsst_model_geopolarblend,2)))
plt.text(26,18,'Corr = ' + str(np.round(corr_geopolarblend,2)))

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Drifters'+' =  '+str(number_obs),color='blue',markersize=10,markeredgecolor='k')

theta = np.arccos(corr_sohcs)
rr = std_subsst_model_sohcs/std_subsst_obs_sohcs
ax1.plot(theta,rr,'s',color = 'orange',markersize=8,markeredgecolor='k',label='SOHCS'+' = '+str(np.round(number_obs_used_sohcs,2)))

theta = np.arccos(corr_coraltemp)
rr = std_subsst_model_coraltemp/std_subsst_obs_coraltemp
ax1.plot(theta,rr,'s',color = 'lime',markersize=8,markeredgecolor='k',label='CoralTemp'+' = '+str(np.round(number_obs_used_coraltemp,2)))

theta = np.arccos(corr_geopolarblend)
rr = std_subsst_model_geopolarblend/std_subsst_obs_geopolarblend
ax1.plot(theta,rr,'s',color = 'magenta',markersize=8,markeredgecolor='k',label='GeoPolarBlend'+' = '+str(np.round(number_obs_used_geopolarblend,2)))

theta = np.arccos(corr_metop)
rr = std_subsst_model_metop/std_subsst_obs_metop
ax1.plot(theta,rr,'s',color = 'cyan',markersize=8,markeredgecolor='k',label='METOP'+' = '+str(np.round(number_obs_used_metop,2)))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('SST ',fontsize=18)

plt.text(1.4,1.2,'Bias',fontsize=14)
plt.text(1.4,1.1,'SOHCS = '+str(np.round(bias_sohcs,3)),fontsize=14)
plt.text(1.4,1.0,'CoralTemp = '+str(np.round(bias_coraltemp,3)),fontsize=14)
plt.text(1.4,0.9,'GeoPolarBlend = '+str(np.round(bias_geopolarblend,3)),fontsize=14)
plt.text(1.4,0.8,'METOP = '+str(np.round(bias_metop,3)),fontsize=14)

#plt.savefig('fig7',bbox_inches = 'tight',pad_inches = 0.1)

##################################################################
