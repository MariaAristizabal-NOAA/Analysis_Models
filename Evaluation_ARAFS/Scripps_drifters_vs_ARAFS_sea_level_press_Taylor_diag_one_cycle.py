#%% User input
# forecasting cycle to be used
cycle = '2026010300'
url_drifter = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/Data/Lagrangian_Drifters_Scripps/LDL_sea_level_press_Jan_2026.nc'

exp_labels = ['UnCoupled','Coupled']
exp_colors = ['blue','orange']

lon_lim = [-180,-80]
lat_lim = [0,70]

folder_exps = ['/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_uncoupled/','/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/ARAFS_Exp4_alaska_4_a_coupled/']

# folder utils for Hycom
#folder_myutils= '/home/Maria.Aristizabal/Utils/'

################################################################################
def get_corresponding_model_array_from_obs_array_grib2_file(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_model_name,lon_model,lat_model,var_model_name):

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
        model = grib2io.open(file,mode='r')
        t = model.select(shortName=var_model_name)[0].validDate
        timestamp = mdates.date2num(t)
        #timestamp = t.timestamp()
        target_time.append(t)

        okt = int(np.interp(timestamp,timestamp_obs_sort,np.arange(len(timestamp_obs_sort))))
        oktt = timestamp_obs_sort[0:okt+1] >= timestamp
        ind = np.arange(okt+1)[oktt]
        sublon_obs = lon_obs_sort[ind]
        sublat_obs = lat_obs_sort[ind]
        subvar_obs.append(var_obs_sort[ind].tolist())

        if lon_model.ndim == 1:
            oklon = np.round(np.interp(sublon_obs,lon_model,np.arange(len(lon_model)))).astype(int)
        if lon_model.ndim == 2:
            oklon = np.round(np.interp(sublon_obs,lon_model[0,:],np.arange(len(lon_model[0,:])))).astype(int)
        if lat_model.ndim == 1:
            oklat = np.round(np.interp(sublat_obs,lat_model,np.arange(len(lat_model)))).astype(int)
        if lat_model.ndim == 2:
            oklat = np.round(np.interp(sublat_obs,lat_model[:,0],np.arange(len(lat_model[:,0])))).astype(int)
        var = model.select(shortName=var_model_name)[0].data[oklat,oklon]
        var[var==0] = np.nan        
        subvar_model.append(var.tolist())        

    return subvar_obs,subvar_model

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
import grib2io
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

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
#%% Read drifter data
url = url_drifter

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
wmo_id = np.asarray(gdata.wmo_ID)
sea_level_pressure = np.asarray(gdata.sea_level_pressure)

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
wmo_idDr = wmo_id[oktimeg]
slpDr = sea_level_pressure[oktimeg]

# Find the different drifter within lat, lon and time window
oklat = np.logical_and(latDr >= lat_lim[0], latDr <= lat_lim[1])
lonDD = lonDr[oklat]
oklon = np.logical_and(lonDD >= lon_lim[0], lonDD <= lon_lim[1])

# Fields within lat and lon window
timeD = timeDr[oklat][oklon]
timestampD = timestampDr[oklat][oklon]
latD = latDr[oklat][oklon]
lonD = lonDr[oklat][oklon]
slpD = slpDr[oklat][oklon]
wmo_idD = wmo_idDr[oklat][oklon]

codes = np.unique(wmo_idD)

########################################################################
number_obs = np.empty(len(folder_exps))
number_obs[:] = np.nan
number_obs_used = np.empty(len(folder_exps))
number_obs_used[:] = np.nan
std_subsst_obs = np.empty(len(folder_exps))
std_subsst_obs[:] = np.nan
std_subsst_model = np.empty(len(folder_exps))
std_subsst_model[:] = np.nan
bias = np.empty(len(folder_exps))
bias[:] = np.nan
corr = np.empty(len(folder_exps))
corr[:] = np.nan

# Loop the experiments
for i,folder in enumerate(folder_exps):
    print(folder)

    #%% Get list files
    files_fv3 = sorted(glob.glob(os.path.join(folder+cycle+'/00E/','arafs*.grb2')))
    model = grib2io.open(files_fv3[0],mode='r')

    lon_fv3 = model.select(shortName='ELON')[0].data
    lat_fv3 = model.select(shortName='NLAT')[0].data

    #################################################################################
    #%% Retrieve HAFS_HYCOM temp. following saildrone trajectory

    files_model = files_fv3
    time_model_name = 'time'
    lat_model = lat_fv3
    lon_model = lon_fv3-360
    var_model_name = 'PRMSL'
    timestamp_obss = timestampD
    lon_obss = lonD
    lat_obss = latD
    var_obss = slpD

    subsst_obs, subsst_model = get_corresponding_model_array_from_obs_array_grib2_file(timestamp_obss,lon_obss,lat_obss,var_obss,files_model,time_model_name,lon_model,lat_model,var_model_name)

    #######################################################################        
    subsst_obs_arr = np.asarray([item for sublist in subsst_obs for item in sublist])
    subsst_model_arr = np.asarray([item for sublist in subsst_model for item in sublist])/100
    ok1 = np.isfinite(subsst_obs_arr)
    subsst_obs_arra = subsst_obs_arr[ok1]
    subsst_model_arra = subsst_model_arr[ok1]
    ok2 = np.isfinite(subsst_model_arra)
    subsst_obs_array = subsst_obs_arra[ok2]
    subsst_model_array = subsst_model_arra[ok2]
    number_obs[i] = lonD.shape[0]
    number_obs_used[i] = len(subsst_obs_array)
    std_subsst_obs[i] = np.nanstd(subsst_obs_array)
    std_subsst_model[i] = np.nanstd(subsst_model_array)
    bias[i] = np.nanmean(subsst_model_array) - np.nanmean(subsst_obs_array)
    corr[i] = np.corrcoef(subsst_obs_array,subsst_model_array)[0,1]
    
    plt.figure(figsize=(8,6))
    plt.plot(subsst_obs_array,subsst_model_array,'.',color=exp_colors[i],markersize=7,markeredgecolor='k')
    plt.plot(np.arange(950,1050),np.arange(950,1050),'-',color='silver',linewidth=2)
    plt.ylim([950,1050])
    plt.xlim([950,1050])
    plt.xlabel('Drifters Observations',fontsize=14)
    plt.ylabel(exp_labels[i],fontsize=14)
    plt.title('Sea Level Pressure',fontsize=18)
    plt.text(1010,980,'Total Observations = ' + str(number_obs[i]))
    plt.text(1010,975,'Observations used = ' + str(number_obs_used[i]))
    plt.text(1010,970,'Bias = ' + str(np.round(bias[i],2)))
    plt.text(1010,965,'STD obs = ' + str(np.round(std_subsst_obs[i],2)))
    plt.text(1010,960,'STD model = ' + str(np.round(std_subsst_model[i],2)))
    plt.text(1010,955,'Corr = ' + str(np.round(corr[i],2)))

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Drifters',color='orangered',markersize=10,markeredgecolor='k')

for m in np.arange(len(folder_exps)):  # loop the models
    print('m=',m)
    theta = np.arccos(corr[m])
    rr = std_subsst_model[m]/std_subsst_obs[m]
    ax1.plot(theta,rr,'o',color = exp_colors[m],markersize=8,markeredgecolor='k',label=exp_labels[m])

plt.legend(loc='upper left',bbox_to_anchor=[0.72,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax1.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('Sea Level Pressure',fontsize=18)
#plt.savefig('fig7',bbox_inches = 'tight',pad_inches = 0.1)

