#!/usr/bin/env python3

"""This script is to compare ARAFS with Argo floats"""

import os
import glob
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

# Parse the yaml config file
print('Parse the config file: plot_atmos4.yml:')
with open('plot_atmos4.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
cycle = conf['ymdh']
experiments = conf['Experiments']
exp_markers = conf['Exp_markers']
exp_colors = conf['Exp_colors']
data_prod = conf['Data_prod']
fhhhs = conf['fhhhs']
commodels = conf['COMmodels']
comobs = conf['COMobs']
fileobs = conf['Fileobs']

xlim = conf['xlim']
ylim = conf['ylim']

cartopy.config['data_dir'] = conf['cartopyDataDir']

print(conf)

################################################################################
def get_var_from_model_for_Argo_float_profile(files_model,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs):

    target_var = np.empty((len(depth_levels)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        model = xr.open_dataset(file)
        if file.split('.')[-1] == 'nc':
            t = model[time_name][:]
            timestamp = mdates.date2num(t)[0]
            


            target_time.append(mdates.num2date(timestamp))



            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])

        # Interpolating lat_obs and lon_obs into model grid
        sublon = np.interp(timestamp,timestamp_obs,lon_obs)
        sublat = np.interp(timestamp,timestamp_obs,lat_obs)
        if lon_model.ndim == 1:
            oksort_lon = np.argsort(lon_model)
            oksort_lat = np.argsort(lat_model)
            lonn = lon_model[oksort_lon]
            latt = lat_model[oksort_lat]
        if lon_model.ndim == 2:
            oksort_lon = np.argsort(lon_model[0,:])
            oksort_lat = np.argsort(lat_model[:,0])
            lonn = lon_model[0,oksort_lon]
            latt = lat_model[oksort_lat,0]
        oklon = int(np.round(np.interp(sublon,lonn,np.arange(len(lonn)))))
        oklat = int(np.round(np.interp(sublat,latt,np.arange(len(latt)))))
        if file.split('.')[-1] == 'nc':
            if model[var_name].ndim == 4:
                var = np.asarray(model[var_name][0,depth_levels,oksort_lat,oksort_lon])
            if model[var_name].ndim == 3:
                var = np.asarray(model[var_name][0,oksort_lat,oksort_lon])
        if len(depth_levels) > 1:
            target_var[:,x] = var[:,oklat,oklon]
        else:
            target_var[x] = var[oklat,oklon]

    return target_time, target_var

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

####################################################################
# Read Surface temp from observations products
ncfile = comobs[0] + '/' + fileobs[0]
OBS = xr.open_dataset(ncfile)

time_obs = np.asarray(OBS['time'])
lat_obs = np.asarray(OBS['latitude'])
lon_obs = np.asarray(OBS['longitude'])
temp_obs = np.asarray(OBS['temp']) 
temp_qc_obs = np.asarray(OBS['temp_qc']) 
salt_obs = np.asarray(OBS['psal']) 
salt_qc_obs = np.asarray(OBS['psal_qc']) 
depth_obs = np.asarray(OBS['pres']) 
depth_qc_obs = np.asarray(OBS['pres_qc']) 

depth_qc_obs = np.asarray(OBS['pres_qc']) 

#timestamp_obs = time_obs.astype('int')/1e9
timestamp_obs = mdates.date2num(time_obs)

total_nobs = len(np.asarray(OBS['platform_number']))
total_plat = len(np.unique(np.asarray(OBS['platform_number'])))
plat_ids = np.unique(np.asarray(OBS['platform_number']))

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.axis('scaled')
plt.plot(lon_obs,lat_obs,'.',label='ARGO = '+ str(total_nobs))
plt.legend()
plt.title(str(time_obs[0]).split('T')[0])

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

'''
okplat = np.asarray(OBS['platform_number']) == plat_ids[0]

timep = time_obs[okplat]
latp = lat_obs[okplat]
lonp = lon_obs[okplat]
tempp = temp_obs[okplat]
tempp_qc = temp_qc_obs[okplat]
saltp = salt_obs[okplat]
saltp_qc = salt_qc_obs[okplat]
depthp = depth_obs[okplat]
depthp_qc = depth_qc_obs[okplat]

fig = plt.figure(figsize=(6,8))
plt.plot(tempp,-depthp,'.-',label='ARGO')
plt.legend()
plt.xlabel('Temperature $^oC$')
plt.ylabel('Depth (m)')
plt.title('Lon = '+str(np.round(lonp[0],2))+' Lat = '+str(np.round(latp[0],2)))
'''

########################################################################
number_obs_temp = np.empty((len(experiments)))
number_obs_temp[:] = np.nan
number_obs_temp_used = np.empty((len(experiments)))
number_obs_temp_used[:] = np.nan

bias_temp = np.empty((len(experiments)))
bias_temp[:] = np.nan
corr_temp = np.empty((len(experiments)))
corr_temp[:] = np.nan
std_temp_obs = np.empty((len(experiments)))
std_temp_obs[:] = np.nan
std_temp_model = np.empty((len(experiments)))
std_temp_model[:] = np.nan

temp_Obs = []
temp_Model = []

for ex,exp in enumerate(experiments):
    print(exp)
    ################################################################
    if exp == 'ARAFS_OCN':
        files_model = np.sort(glob.glob(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/*mom6*.nc'))
         
    times_model = []
    timestamps_model = []
    for x,file in enumerate(files_model):
        model = xr.open_dataset(file)
        t = np.asarray(model['time'][:])
        timestamp = mdates.date2num(t)[0]
        times_model.append(mdates.num2date(timestamp))        
        timestamps_model.append(timestamp)        

    for plat_id in plat_ids:
        print(plat_id)
        okid = np.asarray(OBS['platform_number']) == plat_id
        timestamp_obs_unique = np.unique(timestamp_obs[okid])
        for tt in timestamp_obs_unique:
            oktime = timestamp_obs[okid] == tt 
            timepp = time_obs[okid][oktime]
            timestamppp = timestamp_obs[okid][oktime]
            latpp = lat_obs[okid][oktime]
            lonpp = lon_obs[okid][oktime]
            depthpp = depth_obs[okid][oktime]
            depthpp_qc = depth_qc_obs[okid][oktime]
            temppp = temp_obs[okid][oktime]
            temppp_qc = temp_qc_obs[okid][oktime]
            saltpp = salt_obs[okid][oktime]
            saltpp_qc = salt_qc_obs[okid][oktime]
            
            # Use only observations with qc flag = '1', good data
            good = temppp_qc == '1'
            timep = timepp[good] 
            timestampp = timestamppp[good]
            latp = latpp[good] 
            lonp = lonpp[good] 
            depthp = depthpp[good]
            tempp = temppp[good] 
            saltp = saltpp[good] 

            if len(timep) == 0:
                print('Bad data')
                break

            ok_file_model = int(np.round(np.interp(tt,timestamps_model,np.arange(len(timestamps_model)))))

            file_model = files_model[ok_file_model]
            nc_model = xr.open_dataset(file_model)
            time_model = np.asarray(nc_model['time'])
            lon_model = np.asarray(nc_model['xh'])
            lat_model = np.asarray(nc_model['yh'])

            latp_target = np.unique(latp)[0]
            lonp_target = np.unique(lonp)[0]
            oklat = int(np.round(np.interp(latp_target,lat_model,np.arange(len(lat_model)))))
            oklon = int(np.round(np.interp(lonp_target,lon_model,np.arange(len(lon_model)))))
            temp_model =  np.asarray(nc_model['temp'][0,:,oklat,oklon])
            depth_model = np.asarray(nc_model['z_l'])

            #################################################################
            # Interpolate Argo temp onto model temp
            tempp_int = np.interp(depth_model,depthp,tempp,left=np.nan,right=np.nan)

            '''
            fig = plt.figure(figsize=(6,8))
            plt.plot(tempp,-depthp,'.-',label='ARGO')
            plt.plot(temp_model,-depth_model,'.-',label=exp)
            plt.plot(tempp_int,-depth_model,'o-',label='ARGO int')
            plt.legend()
            plt.xlabel('Temperature $^oC$')
            plt.ylabel('Depth (m)')
            plt.title('Lon = '+str(np.round(lonp[0],2))+' Lat = '+str(np.round(latp[0],2)))
            plt.ylim([np.min(-depthp),0])
            '''

            temp_Obs.append(tempp_int.tolist())
            temp_Model.append(temp_model.tolist())
     
    #################################################################
    # Calculate statistics
    temp_Obs_flat = np.asarray([item for sublist in temp_Obs for item in sublist])
    temp_Model_flat = np.asarray([item for sublist in temp_Model for item in sublist])

    okt_obs = np.isfinite(temp_Obs_flat)
    temp_Obs_flat_fin = temp_Obs_flat[okt_obs]
    temp_Model_flat_fin = temp_Model_flat[okt_obs]

    okt_mod = np.isfinite(temp_Model_flat_fin)
    temp_Obs_flat_finit = temp_Obs_flat_fin[okt_mod]
    temp_Model_flat_finit = temp_Model_flat_fin[okt_mod]

    number_obs_temp[ex] = len(time_obs)
    number_obs_temp_used[ex] = len(temp_Obs_flat_finit)

    bias_temp[ex] = np.mean(temp_Model_flat_finit)-np.nanmean(temp_Obs_flat_finit)
    corr_temp[ex] = np.corrcoef(temp_Model_flat_finit,temp_Obs_flat_finit)[0,1]
    std_temp_obs[ex] = np.std(temp_Obs_flat_finit)
    std_temp_model[ex] = np.std(temp_Model_flat_finit)

    ##################################################################
    # Scatter plots
    fig = plt.figure(figsize=(6,5))
    plt.plot(temp_Obs_flat_finit,temp_Model_flat_finit,'.',label='Total Obs used= '+ str(number_obs_temp_used[ex]))
    plt.plot(np.arange(32),np.arange(32),'-k')
    plt.legend()
    #plt.title(str(time_obs[0]).split('T')[0])
    plt.xlabel('ARGO Floats')
    plt.ylabel(exp)


##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax = taylor_template(angle_lim,std_lim)
ax.plot(0,1,'o',label='ARGO'+' =  '+str(np.round(number_obs_temp[ex],2)),color='red',markersize=10,markeredgecolor='k')

for ex,exp in enumerate(experiments):  # loop the models
    print(exp)
    theta = np.arccos(corr_temp[ex])
    rr = std_temp_model[ex]/std_temp_obs[ex]
    print('corr = ',corr_temp[ex])
    print('normalized std = ',rr)
    ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k',label=experiments[ex]+' =  '+str(np.round(number_obs_temp_used[ex],2)))

plt.legend(loc='upper left',bbox_to_anchor=[0.7,1.0])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
plt.title('Temperature Cycle='+cycle,fontsize=18)

y0 = 1.1
x0 = 1.4
plt.text(x0,y0,'Bias (Model - Obs)',fontsize=14)
for m in np.arange(len(experiments)):  # loop the models
    plt.text(x0,y0-(m+1)*0.1,experiments[m]+' = '+str(np.round(bias_temp[m],3)),fontsize=14)
#plt.savefig('Taylor_diag_Temp',bbox_inches = 'tight',pad_inches = 0.1)




