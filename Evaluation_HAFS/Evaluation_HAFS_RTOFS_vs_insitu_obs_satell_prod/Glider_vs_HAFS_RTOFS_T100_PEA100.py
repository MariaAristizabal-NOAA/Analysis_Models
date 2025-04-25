# User input
# Helene
#cycle = '2024092312'
#storm_num = '09'
#basin = 'al'
#storm_id = '09l'
#storm_name = 'Helene'

# Milton
cycle = '2024100800'
storm_num = '14'
basin = 'al'
storm_id = '14l'
storm_name = 'Milton'

exp_names = ['RTOFS','RTOFS_v2.5.test01','HFSA_oper','HFSB_oper','hafsv2p0p1a_2024rt_cplinit_AOBS']
exp_labels = ['RTOFS_oper','RTOFSv2.5','HFSA_oper','HFSB_oper','HAFS_cplinit_AOBS']
exp_colors = ['magenta','salmon','purple','lime','cyan']
hafs = [' ',' ','hfsa','hfsb','hfsa']
ocean = ['rtofs','rtofs','mom6','hycom','mom6']
markers = ['p','p','o','o','o']

folder_glider = '/scratch2/NCEPDEV/hurricane/noscrub/Maria.Aristizabal/Data/Gliders/2024/' 

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

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_myutils= '/home/Maria.Aristizabal/Utils/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'

################################################################################
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

##################################################################################
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

##################################################################################
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

##################################################################################
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

################################################################################
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
from Upper_ocean_metrics import T100, Potential_energy_anomaly100

from eos80 import dens

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
folder_exps = []
for i in np.arange(len(exp_names)):
    if exp_names[i][0:5] == 'RTOFS':
        folder_exps.append(scratch_folder + exp_names[i] + '/rtofs.' + cycle[:-2] + '/')  
    else:
        folder_exps.append(scratch_folder + exp_names[i] + '/' + cycle + '/' + storm_num + basin[-1] + '/')

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
bath_lat = np.asarray(ncbath.variables['lat'][:])
bath_lon = np.asarray(ncbath.variables['lon'][:])
bath_elev = np.asarray(ncbath.variables['elevation'][:])

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
files_gliders = sorted(glob.glob(os.path.join(folder_glider,'*2024*.nc')))

time_glider = np.empty((len(files_gliders),2000))
time_glider[:] = np.nan
T100_glider = np.empty((len(files_gliders),2000))
T100_glider[:] = np.nan
PEA100_glider = np.empty((len(files_gliders),2000))
PEA100_glider[:] = np.nan
lon_glider = np.empty((len(files_gliders),2000))
lon_glider[:] = np.nan
lat_glider = np.empty((len(files_gliders),2000))
lat_glider[:] = np.nan
dataset_id = []
target_T100_model = np.empty((len(files_gliders),len(folder_exps),len(time_fv3)))
target_T100_model[:] = np.nan
target_PEA100_model = np.empty((len(files_gliders),len(folder_exps),len(time_fv3)))
target_PEA100_model[:] = np.nan
target_time_model = np.empty((len(files_gliders),len(folder_exps),len(time_fv3)))
target_time_model[:] = np.nan

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
        t100gg = np.empty((len(timegg)))
        t100gg[:] = np.nan
        pea100gg = np.empty((len(timegg)))
        pea100gg[:] = np.nan
        for t in np.arange(len(timeg_gridded)):
            temp = tempg_gridded[:,t]
            salt = saltg_gridded[:,t]
            dept = depthg_gridded
            density = dens(salt,temp,dept)
            t100 = T100(dept,temp)
            t100gg[t] = t100
            if len(temp[np.isfinite(temp)]) != 0:
                density = dens(salt,temp,dept)
                pea100gg[t] = Potential_energy_anomaly100(dept,density)
            else:
                pea100gg[t] = np.nan
        
        time_glider[fg,0:len(timegg)] = timegg
        T100_glider[fg,0:len(t100gg)] = t100gg
        PEA100_glider[fg,0:len(pea100gg)] = pea100gg
        lon_glider[fg,0:len(longg)] = longg
        lat_glider[fg,0:len(latgg)] = latgg
    
        ####################################################################
        #%% Loop the experiments
        for i,folder in enumerate(folder_exps):
            print(i)
            print(folder)    
            #%% Get list files
            if ocean[i] == 'hycom':
                files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
                hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
                lon = np.asarray(hafs_ocean['Longitude'][:])
                lat = np.asarray(hafs_ocean['Latitude'][:])
                depth = np.asarray(hafs_ocean['Z'][:])
        
            if ocean[i] == 'mom6':
                files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*mom6*.nc')))
                hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
                lon = np.asarray(hafs_ocean['xh'][:])
                lat = np.asarray(hafs_ocean['yh'][:])
                depth = np.asarray(hafs_ocean['z_l'][:])

            if ocean[i] == 'rtofs':
                files_hafs_ocean = sorted(glob.glob(os.path.join(folder,'*3dz*.nc')))
                hafs_ocean = xr.open_dataset(files_hafs_ocean[0],decode_times=False)
                lon = np.asarray(hafs_ocean['Longitude'][:])
                lat = np.asarray(hafs_ocean['Latitude'][:])
                depth = np.asarray(hafs_ocean['Depth'][:])
        
            #%% Retrieve glider transect
            if np.min(lon) < 0:
                lon_glid = longg
            else:
                lon_glid = lon_glide
            lat_glid = lat_glide

            if ocean[i] == 'hycom':
                ncfiles = files_hafs_ocean
                time_name = 'MT'
                temp_name = 'temperature'
                salt_name = 'salinity'
            
                target_t, target_temp_model_oc, target_salt_model_oc = \
                get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

            if ocean[i] == 'mom6':
                ncfiles = files_hafs_ocean
                time_name = 'time'
                temp_name = 'temp'
                salt_name = 'so'
        
                target_t, target_temp_model_oc, target_salt_model_oc = \
                get_glider_transect_from_HAFS_OCEAN(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

            if ocean[i] == 'rtofs':
                ncfiles = files_hafs_ocean
                time_name = 'MT'
                temp_name = 'temperature'
                salt_name = 'salinity'
            
                target_t, target_temp_model_oc, target_salt_model_oc = \
                get_glider_transect_from_RTOFS_ncfiles(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,lon_glid,lat_glid,tstamp_glider)

            #########################
            # OHC, MLT, MLD
            for t in np.arange(len(target_t)):
                temp = target_temp_model_oc[:,t]
                salt = target_salt_model_oc[:,t]
                if np.ndim(depth) == 1:
                    density = dens(salt,temp,depth)
                    target_T100_model[fg,i,t] = T100(depth,temp)
                    target_PEA100_model[fg,i,t] = Potential_energy_anomaly100(depth,density)

                else:
                    depthh = depth[:,t]
                    density = dens(salt,temp,depthh)
                    target_T100_model[fg,i,t] = T100(depthh,temp)
                    target_PEA100_model[fg,i,t] = Potential_energy_anomaly100(depthh,density)
        
            timestamp = mdates.date2num(target_t)
            target_time_model[fg,i,0:len(target_t)] = timestamp 

    else:
        continue

##################################################################
# Calculate stats

# HAFS vs Gliders
T100_model = [[] for _ in range(len(exp_names))]
T100_glider_int_model = [[] for _ in range(len(exp_names))]
PEA100_model = [[] for _ in range(len(exp_names))]
PEA100_glider_int_model = [[] for _ in range(len(exp_names))]

for i in np.arange(len(exp_names)):
    for fg in np.arange(len(files_gliders)):
        okg = np.isfinite(T100_glider[fg,:])
        t100g = T100_glider[fg,:][okg]
        sort = np.argsort(target_T100_model[fg,i,:])
        okh = np.isfinite(target_T100_model[fg,i,sort])
        t100h = target_T100_model[fg,i,sort][okh]
        if len(t100g)>=2 and len(t100h)>=2:
            tg = time_glider[fg,:][okg]
            th = target_time_model[fg,i,sort][okh]
            okhh = np.logical_and(th >= np.min(tg), th <= np.max(tg))
            T100_glider_int_model[i].append(np.interp(th[okhh],tg,t100g))
            T100_model[i].append(t100h[okhh])

        okg = np.isfinite(PEA100_glider[fg,:])
        peag = PEA100_glider[fg,:][okg]
        sort = np.argsort(target_PEA100_model[fg,i,:])
        okh = np.isfinite(target_PEA100_model[fg,i,sort])
        peah = target_PEA100_model[fg,i,sort][okh]
        if len(peag)>=2 and len(peah)>=2:
            tg = time_glider[fg,:][okg]
            th = target_time_model[fg,i,sort][okh]
            okhh = np.logical_and(th >= np.min(tg), th <= np.max(tg))
            PEA100_glider_int_model[i].append(np.interp(th[okhh],tg,peag))
            PEA100_model[i].append(peah[okhh])

T100_model_flat = [[] for _ in range(len(exp_names))]
T100_glider_int_model_flat = [[] for _ in range(len(exp_names))]
PEA100_model_flat = [[] for _ in range(len(exp_names))]
PEA100_glider_int_model_flat = [[] for _ in range(len(exp_names))]
for i in np.arange(len(exp_names)):
    T100_model_flat[i] = [j for sub in T100_model[i] for j in sub]
    T100_glider_int_model_flat[i] = [j for sub in T100_glider_int_model[i] for j in sub]
    PEA100_model_flat[i] = [j for sub in PEA100_model[i] for j in sub]
    PEA100_glider_int_model_flat[i] = [j for sub in PEA100_glider_int_model[i] for j in sub]

bias_T100_model = [[] for _ in range(len(exp_names))]
std_T100_model = [[] for _ in range(len(exp_names))]
std_T100_glider_model = [[] for _ in range(len(exp_names))]
corr_T100_glider_model = [[] for _ in range(len(exp_names))]
bias_PEA100_model = [[] for _ in range(len(exp_names))]
std_PEA100_model = [[] for _ in range(len(exp_names))]
std_PEA100_glider_model = [[] for _ in range(len(exp_names))]
corr_PEA100_glider_model = [[] for _ in range(len(exp_names))]
for i in np.arange(len(exp_names)):
    bias_T100_model[i] = np.mean(np.asarray(T100_glider_int_model_flat[i]) - np.asarray(T100_model_flat[i]))
    std_T100_model[i] = np.std(T100_model_flat[i])
    std_T100_glider_model[i] = np.std(T100_glider_int_model_flat[i])
    corr_T100_glider_model[i] = np.corrcoef(T100_model_flat[i],T100_glider_int_model_flat[i])[0,1]

    bias_PEA100_model[i] = np.mean(np.asarray(PEA100_glider_int_model_flat[i]) - np.asarray(PEA100_model_flat[i]))
    std_PEA100_model[i] = np.std(PEA100_model_flat[i])
    std_PEA100_glider_model[i] = np.std(PEA100_glider_int_model_flat[i])
    corr_PEA100_glider_model[i] = np.corrcoef(PEA100_model_flat[i],PEA100_glider_int_model_flat[i])[0,1]

##################################################################
#%% Loop the experiments
lon_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),len(time_fv3)))
lat_forec_track[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get storm track from trak atcf files
    if hafs[i] != ' ':
        if hafs[i] == 'hfsa':
            file_track = folder + storm_id+'.' + cycle + '.hfsa.trak.atcfunix'
        if hafs[i] == 'hfsb':
            file_track = folder + storm_id +'.' + cycle + '.hfsb.trak.atcfunix'
        print(file_track)
        okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
        lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], _, _, _ = get_storm_track_and_int(file_track,storm_num)

###################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])
fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
for fg in np.arange(len(files_gliders)):
    if fg == 0:
        plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange',label='Glider Track')
    else:
        plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange')
plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
#plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim(lon_lim)
plt.ylim(lat_lim)
#plt.savefig('All_glider_tracks',bbox_inches='tight')    

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

ngliders = len(T100_glider[np.isfinite(T100_glider)])

fig,ax1 = taylor_template(angle_lim,std_lim)
ax1.plot(0,1,'o',label='Gliders = '+str(ngliders),color='blue',markersize=10,markeredgecolor='k')

for i in np.arange(len(exp_names)):
    theta = np.arccos(corr_T100_glider_model[i])
    rr = std_T100_model[i]/std_T100_glider_model[i]
    ax1.plot(theta,rr,markers[i],color = exp_colors[i],markersize=8,markeredgecolor='k',label=exp_labels[i]+' = '+str(len(T100_model_flat[i])))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('T100 ',fontsize=18)

y0 = 1.1
plt.text(1.4,y0,'Bias',fontsize=14)
for m in np.arange(len(folder_exps)):  # loop the models
    plt.text(1.4,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias_T100_model[m],3)),fontsize=14)
#plt.savefig('OHC_Taylor_diagram',bbox_inches='tight')    

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

ngliders = len(PEA100_glider[np.isfinite(PEA100_glider)])

fig,ax1 = taylor_template(angle_lim,std_lim)
ax1.plot(0,1,'o',label='Gliders = '+str(ngliders),color='blue',markersize=10,markeredgecolor='k')

for i in np.arange(len(exp_names)):
    theta = np.arccos(corr_PEA100_glider_model[i])
    rr = std_PEA100_model[i]/std_PEA100_glider_model[i]
    ax1.plot(theta,rr,markers[i],color = exp_colors[i],markersize=8,markeredgecolor='k',label=exp_labels[i]+' = '+str(len(PEA100_model_flat[i])))

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('PEA100 ',fontsize=18)

y0 = 1.1
plt.text(1.4,y0,'Bias',fontsize=14)
for m in np.arange(len(folder_exps)):  # loop the models
    plt.text(1.4,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias_PEA100_model[m],3)),fontsize=14)
#plt.savefig('MLT_Taylor_diagram',bbox_inches='tight')    

###################################################################
#%% Figure time series OHC
#for fg in np.arange(len(files_gliders)):
for fg in [1,21,25]:
    lev = np.arange(-9000,9100,100)
    okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])
    fig,ax = plt.subplots()
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='silver')
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    for i in np.arange(len(exp_names)):
        plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
    plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
    plt.plot(lon_glider[fg,:], lat_glider[fg,:],'.-',color='orange',label=dataset_id[fg].split('-')[0])
    plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.0])
    #plt.legend()
    plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
    plt.axis('scaled')
    if np.isfinite(np.nanmin(lon_glider[fg,:])):
        plt.xlim([np.nanmin(lon_glider[fg,:])-4,np.nanmax(lon_glider[fg,:])+4])
        plt.ylim([np.nanmin(lat_glider[fg,:])-4,np.nanmax(lat_glider[fg,:])+4])
    #plt.savefig('Glider_track_'+dataset_id[fg].split('-')[0],bbox_inches='tight',dpi=350)    
for fg in [1,21,25]:
    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot(time_glider[fg,:],T100_glider[fg,:],'o-',color='k',label=dataset_id[fg].split('-')[0],markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        sort = np.argsort(target_time_model[fg,i,:])
        plt.plot(target_time_model[fg,i,sort],target_T100_model[fg,i,sort],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
    #plt.legend()
    #plt.legend(loc='upper right',bbox_to_anchor=[1.1,1.0])
    plt.ylabel('T100 ($^oC$)',fontsize=14)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.title('Average Temp. Top 100 meters',fontsize=16)
    plt.grid(True)
    ax.set_xlim(tini,tend)
    #plt.savefig('OHC_time_series_'+dataset_id[fg].split('-')[0],bbox_inches='tight')    

    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot(time_glider[fg,:],PEA100_glider[fg,:],'o-',color='k',label=dataset_id[fg].split('-')[0],markeredgecolor='k')
    for i in np.arange(len(exp_names)):
        sort = np.argsort(target_time_model[fg,i,:])
        plt.plot(target_time_model[fg,i,sort],target_PEA100_model[fg,i,sort],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k')
    #plt.legend()
    #plt.legend(loc='upper right',bbox_to_anchor=[1.1,0.5])
    plt.ylabel('PEA100 ($J/m^3$)',fontsize=14)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.title('Potential Energy Anomaly Top 100 meters',fontsize=16)
    plt.grid(True)
    ax.set_xlim(tini,tend)
    #plt.savefig('OHC_time_series_'+dataset_id[fg].split('-')[0],bbox_inches='tight')    

##################################################################
