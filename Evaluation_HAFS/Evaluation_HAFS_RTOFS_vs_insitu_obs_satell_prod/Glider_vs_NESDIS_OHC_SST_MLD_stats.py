# User input

# User input
# Encompases Helene and Milton 2024
yyyymmddhh_ini = '2024092312'
# number of days for tha aanlysis from yyyymmddhh_ini
delta_days = 48

folder_glider = '/scratch2/NCEPDEV/hurricane/noscrub/Maria.Aristizabal/Data/Gliders/2024/' 
folder_NESDIS = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/NESDIS_OHC/2024/'

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

artofs_folder = scratch_folder + 'RTOFS/'
# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

# folder utils for Hycom 
folder_myutils= '/home/Maria.Aristizabal/Utils/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'

#####################################################################
def glider_data_vector_to_array2(depth,time,var,lat,lon):

    tstamp = mdates.date2num(timee)
    tun = np.unique(tstamp)
    length = []
    for t,tt in enumerate(tun):
        okt = tstamp == tt
        length.append(tstamp[okt].shape)

    depthg = np.empty((np.max(length),len(tun)))
    depthg[:] = np.nan
    timeg = np.empty((np.max(length),len(tun)))
    timeg[:] = np.nan
    varg = np.empty((np.max(length),len(tun)))
    varg[:] = np.nan
    long = np.empty((np.max(length),len(tun)))
    long[:] = np.nan
    latg = np.empty((np.max(length),len(tun)))
    latg[:] = np.nan
    for t,tt in enumerate(tun):
        print(t, tt)
        okt = tstamp == tt
        timeg[:,t] = np.tile(tt,(1,np.max(length)))
        varg[0:len(var[okt]),t] = var[okt]
        depthg[0:len(depth[okt]),t] = depth[okt]
        long[0:len(lon[okt]),t] = lon[okt]
        latg[0:len(lat[okt]),t] = lat[okt]

    return depthg, timeg, varg, latg, long

######################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
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
def get_glider_transect_from_RTOFS_ncfiles(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,long,latg,tstamp_glider):

    target_temp = np.empty((len(depth),len(ncfiles)))
    target_temp[:] = np.nan
    target_salt = np.empty((len(depth),len(ncfiles)))
    target_salt[:] = np.nan
    target_depth = np.empty((len(depth),len(ncfiles)))
    target_depth[:] = np.nan
    target_time = []

    for x,file in enumerate(ncfiles):
        print(x)
        model = xr.open_dataset(file)
        t = model[time_name][:]
        tstamp_model = mdates.date2num(t)[0]
        target_time.append(mdates.num2date(tstamp_model))

        # Interpolating latg and longlider onto RTOFS grid
        sublon = np.interp(tstamp_model,tstamp_glider,long)
        sublat = np.interp(tstamp_model,tstamp_glider,latg)

        oklat, oklon = find_grid_position_hycom(lat,lon,sublat,sublon)

        target_temp[:,x] = np.asarray(model[temp_name][0,:,oklat,oklon])
        target_salt[:,x] = np.asarray(model[salt_name][0,:,oklat,oklon])

    return target_time, target_temp, target_salt

#####################################################################
def get_glider_transect_from_RTOFS_archv(afiles,lat,lon,nz,long,latg,tstamp_glider,jdm,idm):

    #target_t, target_temp_hafs_oc, target_salt_hafs_oc, depth = \
    #get_glider_transect_from_RTOFS_archv(afiles,lat,lon,nz,lon_glid,lat_glid,tstamp_glider,jdm,idm)
    #long = lon_glid
    #latg = lat_glid 
    
    layers = np.arange(0,nz)

    target_temp = np.empty((nz,len(afiles)))
    target_temp[:] = np.nan
    target_salt = np.empty((nz,len(afiles)))
    target_salt[:] = np.nan
    target_depth = np.empty((nz,len(afiles)))
    target_depth[:] = np.nan
    target_time = []

    for x,file in enumerate(afiles):
        print(file)
        lines = [line.rstrip() for line in open(file[:-1]+'b')]
        time_stamp = lines[-1].split()[2]
        hycom_days = lines[-1].split()[3]
        tzero = datetime(1901,1,1,0,0)
        timeRT = tzero+timedelta(float(hycom_days)-1)
        tstamp_model = mdates.date2num(timeRT)
        target_time.append(timeRT)

        # Interpolating latg and longlider onto RTOFS grid
        sublon = np.interp(tstamp_model,tstamp_glider,long)
        sublat = np.interp(tstamp_model,tstamp_glider,latg)

        oklat, oklon = find_grid_position_hycom(lat,lon,sublat,sublon)

        #ztmp = readVar(file[:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
        #target_ztmp = ztmp[oklat,oklon]
        target_ztmp = reatrive_one_value_rtofs_archv(file[:-2],'srfhgt',str(0),jdm,idm,oklat,oklon)*0.01
        for lyr in tuple(layers):
            print(lyr)
            target_temp[lyr,x] = reatrive_one_value_rtofs_archv(file[:-2],'temp',str(lyr+1),jdm,idm,oklat,oklon)
            target_salt[lyr,x] = reatrive_one_value_rtofs_archv(file[:-2],'salin',str(lyr+1),jdm,idm,oklat,oklon)
            #dp = readVar(file[:-2],'archive','thknss',[lyr+1])/2/9806
            target_dp = reatrive_one_value_rtofs_archv(file[:-2],'thknss',str(lyr+1),jdm,idm,oklat,oklon)/9806
            target_ztmp = np.append(target_ztmp,target_dp)

        target_depth[:,x] = np.cumsum(target_ztmp[0:-1]) + np.diff(np.cumsum(target_ztmp))/2

    return target_time, target_temp, target_salt, target_depth

#####################################################################
def reatrive_one_value_rtofs_archv(rtofs_file,var_name,klayer,jdm,idm,oklat,oklon):

    lines = [line.rstrip() for line in open(rtofs_file+'.b')]
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)
    fld = ma.array([],fill_value=1.2676506002282294e+30)

    inFile = rtofs_file + '.a'

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == 'field':
                nheading = n + 1
                #print(line.split()[0])

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == var_name and line.split()[4] == klayer:
                nvar = n - nheading
                #print(nvar)
                #print(n)

    fid = open(inFile,'rb')
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    #fld = np.array(fld)
    #fld = ma.reshape(fld,(jdm,idm))
    target_fld = fld[idm*(oklat)+oklon]

    #fid.seek((nvar)*4*(npad+idm*(oklat)+oklon),0)
    #fld = fid.read(4*(idm*(oklat)+oklon))
    #fld = struct.unpack('>'+str(idm*(oklat)+oklon)+'f',fld)
    #target_fld = fld[idm*oklat+oklon-50:idm*oklat+oklon+50]

    return target_fld

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
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import sys
import os
import glob
import struct

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_OCEAN,\
                            figure_transect_time_vs_depth,\
                            glider_data_vector_to_array,grid_glider_data,\
                            get_var_from_model_following_trajectory

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_dens_crit, OHC_from_profile

from eos80 import dens

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#####################################################################
#%% Time window
date_ini = yyyymmddhh_ini[0:4]+'/'+yyyymmddhh_ini[4:6]+'/'+yyyymmddhh_ini[6:8]+'/'+yyyymmddhh_ini[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(days=48)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#####################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = np.asarray(ncbath.variables['lat'][:])
bath_lon = np.asarray(ncbath.variables['lon'][:])
bath_elev = np.asarray(ncbath.variables['elevation'][:])

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

######################################################################
#%% Read glider data
files_gliders = sorted(glob.glob(os.path.join(folder_glider,'*2024*.nc')))
files_nesdis = sorted(glob.glob(os.path.join(folder_NESDIS,'ohc_na*2024*.nc')))

time_glider = np.empty((len(files_gliders),8000))
time_glider[:] = np.nan
ohc_glider = np.empty((len(files_gliders),8000))
ohc_glider[:] = np.nan
mlt_glider = np.empty((len(files_gliders),8000))
mlt_glider[:] = np.nan
mld_glider = np.empty((len(files_gliders),8000))
mld_glider[:] = np.nan
lon_glider = np.empty((len(files_gliders),8000))
lon_glider[:] = np.nan
lat_glider = np.empty((len(files_gliders),8000))
lat_glider[:] = np.nan
dataset_id = []
target_time_nesdis = np.empty((len(files_gliders),len(files_nesdis)))
target_time_nesdis[:] = np.nan
target_ohc_nesdis = np.empty((len(files_gliders),len(files_nesdis)))
target_ohc_nesdis[:] = np.nan
target_sst_nesdis = np.empty((len(files_gliders),len(files_nesdis)))
target_sst_nesdis[:] = np.nan
target_mld_nesdis = np.empty((len(files_gliders),len(files_nesdis)))
target_mld_nesdis[:] = np.nan

for fg,file_glider in enumerate(files_gliders):
    print(file_glider)
    gdata = xr.open_dataset(file_glider)#,decode_times=False)
    
    dataset_id.append(gdata.id)
    Temperature = np.asarray(gdata.variables['temperature'][:])
    Salinity = np.asarray(gdata.variables['salinity'][:])
    Density = np.asarray(gdata.variables['density'][:])
    Latitude = np.asarray(gdata.latitude)
    Longitude = np.asarray(gdata.longitude)
    Depth = np.asarray(gdata.depth)
    
    Time = np.asarray(gdata.time)
    okt = np.isfinite(Time)
    
    time = Time[okt]
    temperature = Temperature[okt]
    salinity = Salinity[okt]
    density = Density[okt]
    latitude = Latitude[okt]
    longitude = Longitude[okt]
    depth = Depth[okt]

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
    
    if len(timee) != 0:
        ################
        # glider_data_vector_to_array2
        depth = depthh
        time = timee
        var = temperat
        lat = latitud
        lon = longitud
        depthg, timeg, tempg, latg, long = glider_data_vector_to_array2(depth,time,var,lat,lon)
    
        var = salinit
        _, _, saltg, _, _ = glider_data_vector_to_array2(depth,time,var,lat,lon)
    
        ok = np.where(np.isfinite(timeg[0,:]))[0]
        timegg = timeg[0,ok]
        tempgg = tempg[:,ok]
        saltgg = saltg[:,ok]
        depthgg = depthg[:,ok]
        longg = long[0,ok]
        latgg = latg[0,ok]
        
        delta_z = 1     # default value is 0.3
        tempg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(tempgg,timegg,depthgg,delta_z)
        saltg_gridded, timeg_gridded, depthg_gridded = grid_glider_data(saltgg,timegg,depthgg,delta_z)
        
        # Conversion from glider longitude and latitude to HYCOM convention
        target_lonG, target_latG = geo_coord_to_HYCOM_coord(long[0,ok],latg[0,ok])
        lon_glide = target_lonG
        lat_glide = target_latG
    
        #tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
        tstamp_glider = timeg_gridded
        
        #%% OHC time series
        drho = 0.125
        ref_depth = 10
        ohcgg = np.empty((len(timegg)))
        ohcgg[:] = np.nan
        mltgg = np.empty((len(timegg)))
        mltgg[:] = np.nan
        mldgg = np.empty((len(timegg)))
        mldgg[:] = np.nan
        for t in np.arange(len(timeg_gridded)):
            temp = tempg_gridded[:,t]
            salt = saltg_gridded[:,t]
            dept = depthg_gridded
            density = dens(salt,temp,dept)
            mld, mlt, mls = MLD_dens_crit(drho,ref_depth,dept,temp,salt,density)
            mltgg[t] = mlt
            mldgg[t] = mld
            temp[dept<mld] = mlt
            if len(temp[np.isfinite(temp)]) != 0:
                if np.min(temp[np.isfinite(temp)]) < 26:
                    salt[dept<mld] = mls
                    density = dens(salt,temp,dept)
                    ohcgg[t] = OHC_from_profile(dept,temp,density)
                else:
                    ohcgg[t] = np.nan
            else:
                ohcgg[t] = np.nan
        
        time_glider[fg,0:len(timegg)] = timegg
        ohc_glider[fg,0:len(ohcgg)] = ohcgg
        mlt_glider[fg,0:len(mltgg)] = mltgg
        mld_glider[fg,0:len(mldgg)] = mldgg
        lon_glider[fg,0:len(longg)] = longg
        lat_glider[fg,0:len(latgg)] = latgg
    
        #############################################################
        #%% Read NESDIS OHC file 
        
        #%% Read time
        time_nesdis = []
        timestamp_nesdis = []
        for n,file in enumerate(files_nesdis):
            print(file)
            NESDIS = xr.open_dataset(file)
            t = NESDIS['time'][:]
            timestamp = mdates.date2num(t)[0]
            time_nesdis.append(mdates.num2date(timestamp))
            timestamp_nesdis.append(timestamp)
        
        lon_obs = longg
        lat_obs = latgg
        timestamp_obs = timegg
        files_model = files_nesdis 
        var_name = 'ohc'
        time_name = 'time'
        lon_name = 'longitude'
        lat_name = 'latitude'
        target_time_nesdiss, target_ohc_nesdiss = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,var_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name)
    
        var_name = 'sst'
        target_time_nesdiss, target_sst_nesdiss = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,var_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name)

        var_name = 'omld'
        target_time_nesdiss, target_mld_nesdiss = get_var_from_model_following_trajectory(lon_obs,lat_obs,timestamp_obs,files_model,var_name,time_name=time_name,lon_name=lon_name,lat_name=lat_name)

        tstamp = mdates.date2num(target_time_nesdiss)
        okt = np.logical_and(tstamp >= mdates.date2num(tini),tstamp <= mdates.date2num(tend))
    
        target_time_nesdis[fg,0:len(tstamp[okt])] = tstamp[okt]
        target_ohc_nesdis[fg,0:len(target_ohc_nesdiss[okt])] = target_ohc_nesdiss[okt]
        target_sst_nesdis[fg,0:len(target_sst_nesdiss[okt])] = target_sst_nesdiss[okt]
        target_mld_nesdis[fg,0:len(target_mld_nesdiss[okt])] = target_mld_nesdiss[okt]
    
##################################################################
# Calculate stats

# NESDIS vs Gliders
ohc_nesdis = []
sst_nesdis = []
mld_nesdis = []
ohc_glider_int_nesdis = []
mlt_glider_int_nesdis = []
mld_glider_int_nesdis = []
for fg in np.arange(len(files_gliders)):
    #okg = np.isfinite(ohc_glider[fg,:])
    okg = np.logical_and(np.isfinite(ohc_glider[fg,:]),ohc_glider[fg,:]<300)
    okn = np.isfinite(target_ohc_nesdis[fg,:])
    ohcg = ohc_glider[fg,:][okg]
    ohcn = target_ohc_nesdis[fg,:][okn]
    if len(ohcg)>=2 and len(ohcn)>=2:
        tg = time_glider[fg,:][okg]
        tn = target_time_nesdis[fg,:][okn]
        ohc_glider_int_nesdis.append(np.interp(tn,tg,ohcg))
        ohc_nesdis.append(ohcn)

    okg = np.isfinite(mlt_glider[fg,:])
    okn = np.isfinite(target_sst_nesdis[fg,:])
    mltg = mlt_glider[fg,:][okg]
    sstn = target_sst_nesdis[fg,:][okn]
    if len(mltg)>=2 and len(sstn)>=2:
        tg = time_glider[fg,:][okg]
        tn = target_time_nesdis[fg,:][okn]
        mlt_glider_int_nesdis.append(np.interp(tn,tg,mltg))
        sst_nesdis.append(sstn)

    okg = np.isfinite(mld_glider[fg,:])
    okn = np.isfinite(target_mld_nesdis[fg,:])
    mldg = mld_glider[fg,:][okg]
    mldn = target_mld_nesdis[fg,:][okn]
    if len(mldg)>=2 and len(mldn)>=2:
        tg = time_glider[fg,:][okg]
        tn = target_time_nesdis[fg,:][okn]
        mld_glider_int_nesdis.append(np.interp(tn,tg,mldg))
        mld_nesdis.append(mldn)

ohc_nesdis_flat = [j for sub in ohc_nesdis for j in sub]
ohc_glider_int_nesdis_flat = [j for sub in ohc_glider_int_nesdis for j in sub]
bias_ohc_nesdis = np.mean(np.asarray(ohc_glider_int_nesdis_flat) - np.asarray(ohc_nesdis_flat))
std_ohc_nesdis = np.std(ohc_nesdis_flat)
std_ohc_glidern = np.std(ohc_glider_int_nesdis_flat)
corr_ohc_glider_nesdis = np.corrcoef(ohc_nesdis_flat,ohc_glider_int_nesdis_flat)[0,1]

sst_nesdis_flat = [j for sub in sst_nesdis for j in sub]
mlt_glider_int_nesdis_flat = [j for sub in mlt_glider_int_nesdis for j in sub]
bias_sst_nesdis = np.mean(np.asarray(mlt_glider_int_nesdis_flat) - np.asarray(sst_nesdis_flat))
std_sst_nesdis = np.std(sst_nesdis_flat)
std_mlt_glidern = np.std(mlt_glider_int_nesdis_flat)
corr_sst_glider_nesdis = np.corrcoef(sst_nesdis_flat,mlt_glider_int_nesdis_flat)[0,1]

mld_nesdis_flat = [j for sub in mld_nesdis for j in sub]
mld_glider_int_nesdis_flat = [j for sub in mld_glider_int_nesdis for j in sub]
bias_mld_nesdis = np.mean(np.asarray(mld_glider_int_nesdis_flat) - np.asarray(mld_nesdis_flat))
std_mld_nesdis = np.std(mld_nesdis_flat)
std_mld_glidern = np.std(mld_glider_int_nesdis_flat)
corr_mld_glider_nesdis = np.corrcoef(mld_nesdis_flat,mld_glider_int_nesdis_flat)[0,1]

###################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
#okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])
fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
for fg in np.arange(len(files_gliders)):
    if fg == 0:
        plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange',label='Glider Track')
    else:
        plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend()
plt.title(date_ini[0:13]+'-'+date_end[0:13],fontsize=18)
plt.axis('scaled')
plt.xlim(lon_lim)
plt.ylim(lat_lim)
#plt.savefig('All_glider_tracks',bbox_inches='tight')

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

ngliders = len(ohc_glider[np.isfinite(ohc_glider)])

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Gliders = '+str(ngliders),color='blue',markersize=10,markeredgecolor='k')

theta = np.arccos(corr_ohc_glider_nesdis)
rr = std_ohc_nesdis/std_ohc_glidern
ax1.plot(theta,rr,'s',color = 'orange',markersize=8,markeredgecolor='k',label='SOHCS = '+str(len(ohc_nesdis_flat)))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('OHC ',fontsize=18)
plt.text(1.4,1.1,'Bias',fontsize=14)
plt.text(1.4,1.0,'SOHCS = '+str(np.round(bias_ohc_nesdis,2)),fontsize=14)
#plt.savefig('OHC_Taylor_diagram',bbox_inches='tight')    

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

ngliders = len(mlt_glider[np.isfinite(mlt_glider)])

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Gliders = '+str(ngliders),color='blue',markersize=10,markeredgecolor='k')

theta = np.arccos(corr_sst_glider_nesdis)
rr = std_sst_nesdis/std_mlt_glidern
ax1.plot(theta,rr,'s',color = 'orange',markersize=8,markeredgecolor='k',label='SOHCS = '+str(len(sst_nesdis_flat)))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('MLT ',fontsize=18)
plt.text(1.4,1.1,'Bias',fontsize=14)
plt.text(1.4,1.0,'SOHCS = '+str(np.round(bias_sst_nesdis,3)),fontsize=14)
#plt.savefig('MLT_Taylor_diagram',bbox_inches='tight')    

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 2.5

ngliders = len(mld_glider[np.isfinite(mld_glider)])

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Gliders = '+str(ngliders),color='blue',markersize=10,markeredgecolor='k')

theta = np.arccos(corr_mld_glider_nesdis)
rr = std_mld_nesdis/std_mld_glidern
ax1.plot(theta,rr,'s',color = 'orange',markersize=8,markeredgecolor='k',label='SOHCS = '+str(len(mld_nesdis_flat)))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('MLD ',fontsize=18)
plt.text(2.2,2.2,'Bias',fontsize=14)
plt.text(2.2,2.1,'SOHCS = '+str(np.round(bias_mld_nesdis,1)),fontsize=14)
#plt.savefig('MLD_Taylor_diagram',bbox_inches='tight')

###################################################################
#%% Figure time series OHC
#for fg in np.arange(len(files_gliders)):
for fg in [1,21,25]:
    lev = np.arange(-9000,9100,100)
    fig,ax = plt.subplots()
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange',label=dataset_id[fg].split('-')[0])
    plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
    #plt.legend()
    plt.title(date_ini[0:13]+'-'+date_end[0:13],fontsize=18)
    plt.axis('scaled')
    if np.isfinite(np.nanmin(lon_glider[fg,:])):
        plt.xlim([np.nanmin(lon_glider[fg,:])-4,np.nanmax(lon_glider[fg,:])+4])
        plt.ylim([np.nanmin(lat_glider[fg,:])-4,np.nanmax(lat_glider[fg,:])+4])
    #plt.savefig('Glider_track_'+dataset_id[fg].split('-')[0],bbox_inches='tight',dpi=350)    

    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot(time_glider[fg,:],ohc_glider[fg,:],'o-',color='k',label=dataset_id[fg].split('-')[0],markeredgecolor='k')
    plt.plot(target_time_nesdis[fg,:],target_ohc_nesdis[fg,:],'o-',color='orange',label='SOHCS',markeredgecolor='k')
    ax.fill_between(target_time_nesdis[fg,:],target_ohc_nesdis[fg,:]+13.5,target_ohc_nesdis[fg,:]-13.5,color='orange',alpha=0.3)
    #plt.legend()
    #plt.legend(loc='upper right',bbox_to_anchor=[1.1,1.0])
    plt.ylabel('OHC ($kJ/cm^2$)',fontsize=14)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.title('Ocean Heat Content',fontsize=16)
    plt.grid(True)
    ax.set_xlim(tini,tend)
    #plt.savefig('OHC_time_series_'+dataset_id[fg].split('-')[0],bbox_inches='tight')    

    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot(time_glider[fg,:],mlt_glider[fg,:],'o-',color='k',label=dataset_id[fg].split('-')[0],markeredgecolor='k')
    plt.plot(target_time_nesdis[fg,:],target_sst_nesdis[fg,:],'o-',color='orange',label='NESDIS',markeredgecolor='k')
    ax.fill_between(target_time_nesdis[fg,:],target_sst_nesdis[fg,:],target_sst_nesdis[fg,:],color='orange',alpha=0.3)
    #plt.legend()
    #plt.legend(loc='upper right',bbox_to_anchor=[1.1,0.5])
    plt.ylabel('MLT ($^oC$)',fontsize=14)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.title('Mixed Layer Temperature',fontsize=16)
    plt.grid(True)
    ax.set_xlim(tini,tend)
    #plt.savefig('OHC_time_series_'+dataset_id[fg].split('-')[0],bbox_inches='tight')    

    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot(time_glider[fg,:],mld_glider[fg,:],'o-',color='k',label=dataset_id[fg].split('-')[0],markeredgecolor='k')
    plt.plot(target_time_nesdis[fg,:],target_mld_nesdis[fg,:],'o-',color='orange',label='SOHCS',markeredgecolor='k')
    ax.fill_between(target_time_nesdis[fg,:],target_mld_nesdis[fg,:],target_mld_nesdis[fg,:],color='orange',alpha=0.3)
    #plt.legend()
    #plt.legend(loc='upper right',bbox_to_anchor=[1.1,0.5])
    plt.ylabel('MLD (m)',fontsize=14)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.title('Mixed Layer Depth',fontsize=16)
    plt.grid(True)
    ax.set_xlim(tini,tend)

##################################################################
