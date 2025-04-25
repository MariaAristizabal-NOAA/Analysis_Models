#%% User input

# Danielle
YYYYMMDDs = ['20220923','20220924','20220925','20220926','20220927','20220928','20220930']
HHs = ['00','06','12','18']
storm_num = '09'
basin = 'al'

exp_names = ['RTOFS','hafsv0p3a_2022rt_cpl','HWRF','HMON']
exp_labels = ['RTOFS','HAFSv0.3a','HWRF','HMON']
exp_colors = ['b','c','purple','g']

#xlim = [-178.0,-14.96]
#ylim = [-23.0,45.5]

# For OHC calculation
#xlim2 = [-87,-80]
#ylim2 = [15,22]
xlim2 = [-178.0,-14.96]
ylim2 = [-23.0,45.5]

home_folder = '/home/Maria.Aristizabal/'
scratch1_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch2_folder = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch2_folder + exp_names[0] + '/',
               scratch2_folder + exp_names[1] + '/',
               scratch2_folder + exp_names[2] + '/',
               scratch2_folder + exp_names[3] + '/']

# RTOFS grid file name
rtofs_grid_depth_folder = scratch1_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

#best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

#folder_fig = home_folder + 'Analysis/Evaluation_HAFS/Evaluation_HAFSv0p2a_phase3/Laura_2020/Figures/'
#####################################################################
def read_hycom_time(bfile):
    lines = [line.rstrip() for line in open(bfile)]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    time = tzero+timedelta(float(hycom_days))
    timestamp = mdates.date2num(time)
    return time, timestamp

#####################################################################
def read_hycom_field_klayer(rtofs_file,var_name,klayer):

# fhycom.a is assumed to contain idm*jdm 32-bit IEEE real values for
#   each array, in standard f77 element order, followed by padding
#   to a multiple of 4096 32-bit words, but otherwise with no control
#   bytes/words, and input values of 2.0**100 indicating a data void.

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
    #fid.seek((nvar-1)*4*(npad+ijdm),0)
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)

    fld = np.array(fld)
    fld = ma.reshape(fld,(jdm,idm))

    return fld

#####################################################################
def get_profile_rtofs_ab_file_desn_layers(nz,oklon,oklat,ab_file):

# fhycom.a is assumed to contain idm*jdm 32-bit IEEE real values for
#   each array, in standard f77 element order, followed by padding
#   to a multiple of 4096 32-bit words, but otherwise with no control
#   bytes/words, and input values of 2.0**100 indicating a data void.

    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz))
    target_temp_rtofs[:] = np.nan
    target_salt_rtofs = np.empty((nz))
    target_salt_rtofs[:] = np.nan
    target_z_rtofs = np.empty((nz))
    target_z_rtofs[:] = np.nan
    timeRTOFS = []

    lines = [line.rstrip() for line in open(ab_file + '.b')]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days))
    tfmeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT)

    idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
    jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)

    #fld = ma.array([],fill_value=1.2676506002282294e+30)

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == 'field':
                nheading = n + 1
                #print(line.split()[0])

    I = oklat
    J = oklon
    pos = (J-1) * idm + I

    fid = open(ab_file + '.a','rb')

    # Read ztmp
    var_name = 'srfhgt'
    klayer = '0'
    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == var_name and line.split()[4] == klayer:
                nvar = n - nheading 

    print('srfhgt ',nvar)
    fid.seek((nvar)*4*(npad+ijdm+pos),0)
    fld = fid.read(1*4)
    fld1 = struct.unpack('>'+'f',fld)
    target_ztmp = fld1[0]

    #ztmp = readVar(ab_file,'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    #target_ztmp = ztmp[oklat,oklon]

    for lyr in tuple(layers):
        print(lyr)

        # Read temp
        var_name = 'temp'
        for n,line in enumerate(lines):
            if len(line) > 0:
                if line.split()[0] == var_name and line.split()[4] == str(lyr+1):
                    nvar = n - nheading

        print('temp lyr ',lyr, ' ',nvar)
        fid.seek((nvar)*4*(npad+ijdm+pos),0)
        fld = fid.read(1*4)
        if len(fld) != 0:
            fld1 = struct.unpack('>'+'f',fld)
            target_temp_rtofs[lyr] = fld1[0]
        else: 
            target_temp_rtofs[lyr] = np.nan

        # Read salt
        var_name = 'salin'
        for n,line in enumerate(lines):
            if len(line) > 0:
                if line.split()[0] == var_name and line.split()[4] == str(lyr+1):
                    nvar = n - nheading

        print('salt lyr ',lyr, ' ',nvar)
        fid.seek((nvar)*4*(npad+ijdm+pos),0)
        fld = fid.read(1*4)
        if len(fld) != 0:
            fld1 = struct.unpack('>'+'f',fld)
            target_salt_rtofs[lyr] = fld1[0]
        else: 
            target_salt_rtofs[lyr] = np.nan

        # Read dp
        var_name = 'thknss'
        for n,line in enumerate(lines):
            if len(line) > 0:
                if line.split()[0] == var_name and line.split()[4] == str(lyr+1):
                    nvar = n - nheading

        print('dp lyr ',lyr, ' ',nvar)
        fid.seek((nvar)*4*(npad+ijdm+pos),0)
        fld = fid.read(1*4)
        if len(fld) != 0:
            fld1 = struct.unpack('>'+'f',fld)
            target_ztmp = np.append(target_ztmp,fld1[0]/9806)
        else: 
            target_ztmp = np.append(target_ztmp,np.nan)

    target_z_rtofs = np.cumsum(target_ztmp[0:-1]) + np.diff(np.cumsum(target_ztmp))/2

    time_rtofs = np.asarray(timeRTOFS)

    return target_temp_rtofs, target_temp_rtofs, target_z_rtofs, time_rtofs

#####################################################################
'''
temp_rtofs, salt_rtofs, z_rtofs, _ = get_profile_rtofs_ab_file_desn_layers(nz,oklonr[0],oklatr[0],artofs_file[:-2])
ab_file = artofs_file[:-2]
lyr = 1
temp_RTOFS = readVar(ab_file,'archive','temp',[lyr+1])
temp_RTOFS[1209,1348]

########

afile = artofs_file[:-2]
var_name = 'temp'
klayer = '1'

lines = [line.rstrip() for line in open(afile+'.b')]
ijdm = idm*jdm
npad = 4096-(ijdm%4096)
fld = ma.array([],fill_value=1.2676506002282294e+30)
#fld2 = ma.array([],fill_value=1e30)

inFile = afile + '.a'
#inFile = afile

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

I = 1348
J = 1209
pos = (J-1) * idm + I
fid = open(inFile,'rb')
fid.seek((nvar)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld[pos]
#fld = np.array(fld)
#fld = ma.reshape(fld,(jdm,idm))
#fld[1209,1348]

I = 1348
J = 1209
pos = (J-1) * idm + I
fid = open(inFile,'rb')
fid.seek((nvar)*4*(npad+ijdm+pos),0)
fld = fid.read(1*4)
fld2 = struct.unpack('>'+'f',fld)

'''
#####################################################################
def depth_aver_top_100(depth,var):

    varmean100 = np.empty(var.shape[1])
    varmean100[:] = np.nan
    if depth.ndim == 1:
        okd = np.abs(depth) <= 100
        if len(depth[okd]) != 0:
            for t in np.arange(var.shape[1]):
                if len(np.where(np.isnan(var[okd,t]))[0])>10:
                    varmean100[t] = np.nan
                else:
                    varmean100[t] = np.nanmean(var[okd,t],0)
    else:
        for t in np.arange(depth.shape[1]):
            okd = np.abs(depth[:,t]) <= 100
            if len(depth[okd,t]) != 0:
                if len(np.where(np.isnan(var[okd,t]))[0])>10:
                    varmean100[t] = np.nan
                else:
                    varmean100[t] = np.nanmean(var[okd,t])
            else:
                varmean100[t] = np.nan

    return varmean100array

#####################################################################
def taylor_template(angle_lim,std_lim):

    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure(figsize=(7,7))
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
import sys
import os
import glob
import numpy as np
import numpy.ma as ma
import xarray as xr
import struct
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import seawater as sw

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine 

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar, parse_b

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#################################################################################
# Reading RTOFS grid from ab files

print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(rtofs_grid_file+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lonn_rtofs = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
latt_rtofs = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth_rtofs = np.asarray(readdepth(rtofs_depth_file,'depth'))

#oklon_rtofs = np.logical_and(lonn_rtofs[0,:] >= xlim[0]+360,lonn_rtofs[0,:] <= xlim[1]+360)
#oklat_rtofs = np.logical_and(latt_rtofs[:,0] >= ylim[0],latt_rtofs[:,0] <= ylim[1])

#lon_rtofs = lonn_rtofs[:,oklon_rtofs][oklat_rtofs,:]
#lat_rtofs = latt_rtofs[:,oklon_rtofs][oklat_rtofs,:]

oklon_rtofs2 = np.logical_and(lonn_rtofs[0,:] >= xlim2[0]+360,lonn_rtofs[0,:] <= xlim2[1]+360)
oklat_rtofs2 = np.logical_and(latt_rtofs[:,0] >= ylim2[0],latt_rtofs[:,0] <= ylim2[1])

lon_rtofs2 = lonn_rtofs[:,oklon_rtofs2][oklat_rtofs2,:]
lat_rtofs2 = latt_rtofs[:,oklon_rtofs2][oklat_rtofs2,:]

#####################################################################
# Reading HAFS-HYCOM grid

print('Retrieving coordinates from HAFS-HYCOM')
hafs_hycom_folder = scratch2_folder + exp_names[1] + '/' + YYYYMMDDs[-1] + HHs[0]
files_hafs_hycom = sorted(glob.glob(os.path.join(hafs_hycom_folder,'*3z*.nc')))
hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
lonn_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
latt_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

xlimh2, ylimh2  = geo_coord_to_HYCOM_coord(xlim2,ylim2)

if np.min(lonn_hafs_hycom) < 0:
    oklon_hafs2 = np.where(np.logical_and(lonn_hafs_hycom>=xlim2[0],lonn_hafs_hycom<=xlim2[1]))[0]
else:
    oklon_hafs2 = np.where(np.logical_and(lonn_hafs_hycom>=xlimh2[0],lonn_hafs_hycom<=xlimh2[1]))[0]

oklat_hafs2 = np.where(np.logical_and(latt_hafs_hycom>=ylimh2[0],latt_hafs_hycom<=ylimh2[1]))[0]

lon_hafs_hycom2 = lonn_hafs_hycom[oklon_hafs2]
lat_hafs_hycom2 = latt_hafs_hycom[oklat_hafs2]

#####################################################################
'''
# Reading files

cycle_info = np.empty((len(YYYYMMDDs),len(HHs),len(files_hafs_hycom)))
cycle_info[:] = np.nan

ohc_skill = np.empty((len(YYYYMMDDs),len(HHs),len(files_hafs_hycom),3,5))
ohc_skill[:] = np.nan

ohc_spatial_bias = np.empty((len(YYYYMMDDs),len(HHs),len(files_hafs_hycom),len(lat_hafs_hycom2),len(lon_hafs_hycom2)))
ohc_spatial_bias[:] = np.nan

# cycle though RTOFS daily cycle
#for daily_cycle,yyyymmdd in enumerate(YYYYMMDDs[0:1]):
for daily_cycle,yyyymmdd in enumerate(YYYYMMDDs):

    print('reading RTOFS forecast cycle ',yyyymmdd)
    rtofs_folder = scratch2_folder + exp_names[0] + '/' + yyyymmdd + '/'
    artofs_files = np.asarray(sorted(glob.glob(os.path.join(rtofs_folder,'*archv.a'))))
    ncrtofs_files = np.asarray(sorted(glob.glob(os.path.join(rtofs_folder,'*TS.nc'))))
    
    # Read RTOFS time
    time_rtof = []
    timestamp_rtof = []
    for ncrtofs_file in ncrtofs_files:
        rtofs = xr.open_dataset(ncrtofs_file)
        t = rtofs['MT'][:]
        timestamp_rt = mdates.date2num(t)[0]
        timestamp_rtof.append(mdates.date2num(t)[0])
        time_rtof.append(mdates.num2date(timestamp_rt))

    time_rtof = np.asarray(time_rtof)
    timestamp_rtof = np.asarray(timestamp_rtof)

    # Need to sort time and files accordingly
    sort_time = np.argsort(timestamp_rtof)
    time_rtofs = time_rtof[sort_time]
    timestamp_rtofs = timestamp_rtof[sort_time]
    artofs_files_sorted = artofs_files[sort_time]
    ncrtofs_files_sorted = ncrtofs_files[sort_time]

    # Cycle though HAFS forecast 6-hourly cycles
    #for hour_cycle,hh in enumerate(HHs[1:2]):
    for hour_cycle,hh in enumerate(HHs):
        print('reading HAFS forecast cycle ',yyyymmdd+hh)
        hafs_hycom_folder = scratch2_folder + exp_names[1] + '/' + yyyymmdd + hh + '/'
        files_hafs_hycom = sorted(glob.glob(os.path.join(hafs_hycom_folder,'*3z*.nc')))
    
        if len(files_hafs_hycom) != 0:
            for fhour,file in enumerate(files_hafs_hycom):
                print(file)
                hycom = xr.open_dataset(file)
                t = hycom['MT'][:]
                timestamp = mdates.date2num(t)[0]
                time_hafs_hycom = mdates.num2date(timestamp)

                ok_rtofs = np.where(timestamp_rtofs == timestamp)[0][0]

                ncrtofs_file = ncrtofs_files_sorted[ok_rtofs]
                print(ncrtofs_file)
                rtofs = xr.open_dataset(ncrtofs_file)
                zlev_rtofs = np.asarray(rtofs['Depth'])
                temp_rtofs = np.asarray(rtofs['pot_temp'])[0,:,:,:][:,oklat_rtofs2,:][:,:,oklon_rtofs2]
                salt_rtofs = np.asarray(rtofs['salinity'])[0,:,:,:][:,oklat_rtofs2,:][:,:,oklon_rtofs2]
                zlev_array = np.reshape(np.tile(zlev_rtofs,(temp_rtofs.shape[1]*temp_rtofs.shape[2],1)).T,(temp_rtofs.shape[0],temp_rtofs.shape[1],temp_rtofs.shape[2]))
                cp = 3985 #Heat capacity in J/(kg K)
                no26 = temp_rtofs < 26
                temp_rtofs[no26] = np.nan
                salt_rtofs[no26] = np.nan
                dens_rtofs = sw.dens(salt_rtofs,temp_rtofs,zlev_array)
                rho0 = np.nanmean(dens_rtofs,axis=0)
                zlev_array_fac = (zlev_array[0:-1,:,:] + zlev_array[1:,:,:])/2
                zero_array = np.zeros((1,temp_rtofs.shape[1],temp_rtofs.shape[2]))
                bott_array = np.ones((1,temp_rtofs.shape[1],temp_rtofs.shape[2]))* zlev_array_fac[-1,0,0] + (zlev_rtofs[-1] - zlev_rtofs[-2])
                zlev_array_face = np.vstack((zero_array,zlev_array_fac,bott_array))
                dz_array = np.diff(zlev_array_face,axis=0)
                ohc_rtofs = np.abs(cp * rho0 * np.nansum((temp_rtofs-26)*dz_array,axis=0))
                ohc_rtofs = ohc_rtofs * 10**(-7) # in kJ/cm^2

                ohc_hafs_hycom = np.asarray(hycom['ocean_heat_content'][0,oklat_hafs2,:][:,oklon_hafs2])

                # Store cycle info 
                ff = fhour*6
                if len(str(ff))==1:
                    fh = '00' + str(ff)
                if len(str(ff))==2:
                    fh = '0' + str(ff)
                if len(str(ff))==3:
                    fh = str(ff)
                cycle_info[daily_cycle,hour_cycle,fhour] = int(yyyymmdd + hh + fh)

                # Spatial bias
                ohc_spatial_bias[daily_cycle,hour_cycle,fhour,:,:] = ohc_rtofs - ohc_hafs_hycom
                  
                # Define date frames
                df_rtofs_hafs_ohc = pd.DataFrame(data=np.array([np.ravel(ohc_rtofs,order='F'),
                                          np.ravel(ohc_hafs_hycom,order='F')]).T,
                                        columns=['ohc_rtofs','ohc_hafs'])                

                #%% SST and SSS statistics.
                df_hafs = df_rtofs_hafs_ohc.dropna()
                Nhafs = len(df_hafs)-1  #For Unbiased estimmator.

                cols = ['CORRELATION','RTOFS_STD','HAFS_STD','CRMSE','BIAS']

                #CORR
                ohc_skill[daily_cycle,hour_cycle,fhour,0,0] = df_hafs.corr()['ohc_rtofs']['ohc_hafs']

                #OSTD
                ohc_skill[daily_cycle,hour_cycle,fhour,0,1] = df_hafs.std().ohc_rtofs

                #MSTD
                ohc_skill[daily_cycle,hour_cycle,fhour,0,2] = df_hafs.std().ohc_hafs

                #CRMSE
                ohc_skill[daily_cycle,hour_cycle,fhour,0,3] = np.sqrt(np.nansum(((df_hafs.ohc_rtofs-df_hafs.mean().ohc_rtofs)-(df_hafs.ohc_hafs-df_hafs.mean().ohc_hafs))**2)/Nhafs)

                #BIAS
                ohc_skill[daily_cycle,hour_cycle,fhour,0,4] = df_hafs.mean().ohc_rtofs - df_hafs.mean().ohc_hafs

                ohc_skillscores = pd.DataFrame(ohc_skill[daily_cycle,hour_cycle,fhour,:,:],index=['HAFS','HWRF','HMON'],columns=cols)
                print(ohc_skillscores)

ohc_spatial_bias_mean =  np.nanmean(np.nanmean(np.nanmean(ohc_spatial_bias,axis=0),axis=0),axis=0)

np.save('ohc_skill_array',ohc_skill)
np.save('ohc_spatial_bias_array',ohc_spatial_bias_mean)
'''

ohc_skill = np.load('ohc_skill_array.npy')
ohc_spatial_bias_mean = np.load('ohc_spatial_bias_array.npy')

################################ Figures #####################

# Plot OHC all domain RTOFS
kw = dict(levels=np.arange(0,161,20))
plt.figure(figsize=(5, 4))
plt.contour(lon_rtofs2-360,lat_rtofs2,ohc_rtofs,[100],colors='grey',alpha=0.5)
plt.contourf(lon_rtofs2-360,lat_rtofs2,ohc_rtofs,cmap='Spectral_r',extend='max',**kw)
cbar = plt.colorbar(pad=0.008,extendrect=True)
cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
plt.axis('scaled')
#plt.ylim([np.min(lat_hafs_hycom),np.max(lat_hafs_hycom)])
#plt.xlim([np.min(lon_hafs_hycom),np.max(lon_hafs_hycom)])
plt.title('OHC RTOFS ' + ncrtofs_file.split('/')[-1])
#plt.title('OHC RTOFS \n ' + ncrtofs_file.split('/')[-2] + ' ' + artofs_file.split('/')[-1] )
plt.savefig('fig1',bbox_inches = 'tight',pad_inches = 0.1)

# Plot OHC all domain HAFS
kw = dict(levels=np.arange(0,161,20))
plt.figure(figsize=(5,4))
plt.contour(lon_hafs_hycom2,lat_hafs_hycom2,ohc_hafs_hycom,[100],colors='grey',alpha=0.5)
plt.contourf(lon_hafs_hycom2,lat_hafs_hycom2,ohc_hafs_hycom,cmap='Spectral_r',extend='max',**kw)
cbar = plt.colorbar(pad=0.008,extendrect=True)
cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
plt.axis('scaled')
plt.ylim([np.min(lat_hafs_hycom2),np.max(lat_hafs_hycom2)])
plt.xlim([np.min(lon_hafs_hycom2),np.max(lon_hafs_hycom2)])
plt.title('SST HAFS ' + yyyymmdd + hh + ' f' + str(fhour*6))
plt.savefig('fig2',bbox_inches = 'tight',pad_inches = 0.1)

# Plot OHC spatial bias all domain HAFS vs RTOFS
kw = dict(levels=np.arange(-10,11,2))
plt.figure(figsize=(10,4))
#plt.contour(lon_hafs_hycom2,lat_hafs_hycom2,ohc_spatial_bias_mean,[0],colors='grey')
plt.contourf(lon_hafs_hycom2,lat_hafs_hycom2,ohc_spatial_bias_mean,cmap='bwr',**kw)
cbar = plt.colorbar(pad=0.008,extend='both')
cbar.ax.set_ylabel('$Kg/cm^2$',fontsize = 14)
plt.axis('scaled')
plt.ylim([np.min(lat_hafs_hycom2),np.max(lat_hafs_hycom2)])
plt.xlim([np.min(lon_hafs_hycom2),np.max(lon_hafs_hycom2)])
#plt.title('OHC Bias RTOFS - HAFS ' + yyyymmdd + hh + ' f' + str(fhour*6))
plt.title('OHC Bias RTOFS - HAFS all cycles')
plt.savefig('fig3',bbox_inches = 'tight',pad_inches = 0.1)

###########################################################
n_daily_cycle = 0
n_6hourly_cycle = 1
n_fhour_forecast = 0

xticks = np.arange(0,23,4)
xtick_labels = ['f00','f24','f48','f72','f94','f118']

# Plot correlation of each model vs RTOFS, one cycle
plt.figure(figsize=(10,5))
for m in np.arange(ohc_skill.shape[3]):
    plt.plot(np.arange(22),ohc_skill[n_daily_cycle,n_6hourly_cycle,:,m,0],'o-',label=exp_labels[m+1],color=exp_colors[m+1],markeredgecolor='k')
plt.xticks(xticks,xtick_labels)
plt.grid(True)
plt.legend(loc='upper left',bbox_to_anchor=[0.85,1.15])
plt.title('OHC Correlation vs RTOFS cycle ' + str(int(cycle_info[n_daily_cycle,n_6hourly_cycle,0]))[0:-3],fontsize=18)
plt.savefig('fig4',bbox_inches = 'tight',pad_inches = 0.1)

# Plot OHC for each model, one cycle
plt.figure(figsize=(10,5))
# STD RTOFS
plt.plot(np.arange(22),ohc_skill[n_daily_cycle,n_6hourly_cycle,:,0,1],'o-',label=exp_labels[0],color=exp_colors[0],markeredgecolor='k')
for m in np.arange(ohc_skill.shape[3]):
    plt.plot(np.arange(22),ohc_skill[n_daily_cycle,n_6hourly_cycle,:,m,2],'o-',label=exp_labels[m+1],color=exp_colors[m+1],markeredgecolor='k')
plt.xticks(xticks,xtick_labels)
plt.grid(True)
plt.ylabel('$^oC$',fontsize=16)
plt.legend(loc='upper left',bbox_to_anchor=[0.85,1.15])
plt.title('OHC STD cycle ' + str(int(cycle_info[n_daily_cycle,n_6hourly_cycle,0]))[0:-3],fontsize=18)
plt.savefig('fig5',bbox_inches = 'tight',pad_inches = 0.1)


'''
# Plot CRMSE vs fhour, one cycle
plt.figure(figsize=(10,5))
for m in np.arange(sst_skill.shape[3]):
    plt.plot(np.arange(22),sst_skill[n_daily_cycle,n_6hourly_cycle,:,m,3],'o-',label=exp_labels[m+1],color=exp_colors[m+1],markeredgecolor='k')
plt.xticks(xticks,xtick_labels)
plt.grid(True)
plt.legend(loc='upper left',bbox_to_anchor=[0.85,1.15])
plt.title('SST CRMSE vs RTOFS cycle ' + str(int(cycle_info[n_daily_cycle,n_6hourly_cycle,0]))[0:-3],fontsize=18)
'''

'''
# Plot Bias for each model vs RTOFS, one cycle
plt.figure(figsize=(10,5))
for m in np.arange(sst_skill.shape[3]):
    plt.plot(np.arange(22),sst_skill[n_daily_cycle,n_6hourly_cycle,:,m,4],'o-',label=exp_labels[m+1],color=exp_colors[m+1],markeredgecolor='k')
plt.plot(np.arange(22),np.tile(0,22),'--k',linewidth=2)
plt.xticks(xticks,xtick_labels)
plt.grid(True)
plt.ylabel('$^oC$',fontsize=16)
plt.legend(loc='upper left',bbox_to_anchor=[0.85,1.15])
plt.title('SST Bias vs RTOFS cycle ' + str(int(cycle_info[n_daily_cycle,n_6hourly_cycle,0]))[0:-3],fontsize=18)
'''

# Plot OHC Bias for HAFS vs RTOFS, all cycles
bias_max = np.nanmax(np.nanmax(ohc_skill[:,:,:,0,4],axis=0),axis=0)
bias_min = np.nanmin(np.nanmin(ohc_skill[:,:,:,0,4],axis=0),axis=0)

plt.figure(figsize=(10,5))
plt.plot(np.arange(22),np.tile(0,22),'--k',linewidth=2)
for d1 in np.arange(ohc_skill.shape[0]): # loop the daily cycles
    for c6 in np.arange(ohc_skill.shape[1]): # loop 6 hourly cycles 
        plt.plot(np.arange(22),ohc_skill[d1,c6,:,0,4],'o',color=exp_colors[1],markeredgecolor='k')
plt.fill_between(np.arange(22),bias_min,bias_max,color=exp_colors[1],alpha=0.3)
plt.plot(np.arange(22),ohc_skill[d1,c6,:,0,4],'o',label=exp_labels[1],color=exp_colors[1],markeredgecolor='k')
plt.xticks(xticks,xtick_labels)
plt.grid(True)
plt.ylabel('$Kg/cm^2$',fontsize=16)
plt.legend(loc='upper left',bbox_to_anchor=[0.85,1.15])
plt.title('OHC Bias vs RTOFS all cycles',fontsize=18)
plt.savefig('fig6',bbox_inches = 'tight',pad_inches = 0.1)

#%% Combine OHC into one normalized Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='RTOFS',color='red',markersize=10,markeredgecolor='k')

scores = ohc_skill
for m in np.arange(scores.shape[3]):  # loop the models
    print('m=',m)
    for d1 in np.arange(scores.shape[0]): # loop the daily cycles
        for c6 in np.arange(scores.shape[1]): # loop 6 hourly cycles 
            for f in np.arange(scores.shape[2]):  # loop the forecast hours
                print('f=',f)
                theta = np.arccos(scores[d1,c6,f,m,0])
                rr = scores[d1,c6,f,m,1]/scores[d1,c6,f,m,2]
                if (d1==0 and c6==0) and f==0:
                    ax1.plot(theta,rr,markers[m],color = exp_colors[m+1],markersize=8,markeredgecolor='k',label=exp_labels[m+1])
                else:
                    ax1.plot(theta,rr,markers[m],color = exp_colors[m+1],markersize=8,markeredgecolor='k')

plt.legend(loc='upper left',bbox_to_anchor=[0.75,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('OHC All cycles ',fontsize=18)
#file = folder_fig + 'Taylor_norm_2019082800_v2'
plt.savefig('fig7',bbox_inches = 'tight',pad_inches = 0.1)

#####################################################################

