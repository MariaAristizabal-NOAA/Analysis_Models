#!/usr/bin/env python3

"""This script is to compare ARAFS surface temperature with OSTIA SST"""

import os
import glob
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from datetime import datetime, timedelta, timezone
import grib2io
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.interpolate import griddata
import time

# Increase fontsize of labels globally
plt.rcParams.update({
    'font.size': 14,          # General font size
    'axes.titlesize': 20,     # Title size
    'axes.labelsize': 16,     # X and Y label size
    'xtick.labelsize': 14,    # X-axis tick size
    'ytick.labelsize': 14,    # Y-axis tick size
    'legend.fontsize': 14     # Legend size
})

# Parse the yaml config file
print('Parse the config file: plot_atmos_ERA5.yml:')
with open('plot_atmos_ERA5.yml', 'rt') as f:
    conf = yaml.safe_load(f)
#conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
exp_labels = conf['exp_labels']
exp_colors = conf['exp_colors']
exp_markers = conf['exp_markers']
data_prod = conf['Data_prod']
ymdhs = conf['ymdhs']
fhhhs = conf['fhhhs']
commodels = conf['COMmodels']
comobs = conf['COMobs']
fileobs_prefix = conf['Fileobs_prefix']

xlim = conf['xlim']
ylim = conf['ylim']

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
"""
mean_tmp_model = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
mean_tmp_model[:] = np.nan
std_tmp_model = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
std_tmp_model[:] = np.nan

mean_skt_obs = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
mean_skt_obs[:] = np.nan
std_skt_obs = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
std_skt_obs[:] = np.nan

number_obs = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
number_obs[:] = np.nan

bias = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
bias[:] = np.nan
corr = np.empty((len(exp_labels),len(ymdhs),len(fhhhs)))
corr[:] = np.nan

bias_spatial = np.empty((len(exp_labels),len(ymdhs),len(fhhhs),721,1440))
bias_spatial[:] = np.nan

for ex,exp in enumerate(exp_labels):
    print(exp)
    for y,ymdh in enumerate(ymdhs):
        print(ymdh)
        start = time.perf_counter()
        for f,fhhh in enumerate(fhhhs):
            print(fhhh)

            valid_time = datetime.strptime(ymdh,'%Y%m%d%H').replace(tzinfo=timezone.utc) + timedelta(hours=int(fhhh[1:]))
            print('valid time ',valid_time)

            # Read the correct file
            ym_obs = valid_time.strftime('%Y%m') 
            ncfile = glob.glob(comobs+'/'+fileobs_prefix+ym_obs+'*.nc')[0]
            OBS = xr.open_dataset(ncfile)
    
            times_obs = np.asarray(OBS['time'])
            timestamps_obs = times_obs.astype('int')/1e9
            lat_obs = np.asarray(OBS['latitude'])
            lon_obs = np.asarray(OBS['longitude'])
            skts_obs = np.asarray(OBS['SKT']) - 273.15 # convert K to degC

            # Find the correct timestamps
            okt = timestamps_obs == valid_time.timestamp()
            time_obs = times_obs[okt]
            skt_obs = skts_obs[okt,:,:][0,:,:] 
            timestamp_obs = timestamps_obs[okt]
    
            fname = conf['stormModel'].lower()+'.'+ymdh+'.'+fhhh+'.grb2'
            grib2file = os.path.join(conf['COMmodels'][ex]+ymdh+'/00E/', fname)
            print(f'grib2file: {grib2file}')
            grb = grib2io.open(grib2file,mode='r')
     
            print('Extracting lat, lon')
            lat = grb.select(shortName='NLAT')[0].data
            lon = grb.select(shortName='ELON')[0].data
     
            print('Extracting skin temperature')
            level = 'surface'
            tmp = grb.select(shortName='TMP',level=level)[0].data
            land = grb.select(shortName='LAND',level=level)[0].data
            tmp = tmp - 273.15 # convert K to degC

            land[land==1] = np.nan
            land[land==0] = 1
            tmp = tmp * land
     
            '''
            # Transform to geographic coordinates
            lon_geo = np.copy(lon)
            lon_geo[lon_geo>180] = lon_geo[lon_geo>180] - 360

            lon_geo_sort = np.empty((lon.shape[0],lon.shape[1]))
            lon_geo_sort[:] = np.nan
            wtmp_sort = np.empty((lon.shape[0],lon.shape[1]))
            wtmp_sort[:] = np.nan
            for yy in np.arange(lon_geo.shape[0]):
                print(yy)
                sort = np.argsort(lon_geo[yy,:])
                lon_geo_sort[yy,:] = lon_geo[yy,sort]
                wtmp_sort[yy,:] = wtmp[yy,sort]
            '''

            # ERA5 coordinates equal to FV3 coordinates
            '''
            lon_geo = np.copy(lon_obs)
            lon_geo[lon_geo<0] = lon_geo[lon_geo<0] + 360
            sort = np.argsort(lon_geo)
            lon_obs_sort = lon_geo[sort]
            sst_obs_sort = sst_obs[:,sort]
            '''
            lon_obs_sort = lon_obs
            skt_obs_sort = skt_obs

            '''
            # subset
            oklon = np.logical_and(lon[0,:]>=xlim[0]+360,lon[0,:]<xlim[1]+360)
            oklat = np.logical_and(lat[:,0]>=ylim[0],lat[:,0]<ylim[1])
     
            lon_sub = lon[oklat,:][:,oklon]
            lat_sub = lat[oklat,:][:,oklon]
            tmp_sub = wtmp[oklat,:][:,oklon]
            '''

            lon_sub = lon
            lat_sub = lat
            tmp_sub = tmp

            '''
            levels = np.arange(-5,36,5)
            plt.figure()
            plt.contourf(lon_sub,lat_sub,tmp_sub,levels=levels)
            plt.colorbar()

            levels = np.arange(-5,36,5)
            plt.figure()
            plt.contourf(lon_obs_sort,lat_obs,sst_obs_sort,levels=levels)
            plt.contourf(lon_obs_sort,lat_obs,tmp_sub_interp,levels=levels,cmap='jet')
            plt.colorbar()

            levels = np.arange(-5,36,5)
            plt.figure()
            plt.contourf(lon_obs_sort,lat_obs,tmp_sub_interp,levels=levels)
            plt.colorbar()
            '''

            ################################################################
            # interpolation from fv3 to CRW grid
            # --- Original grid ---
            if lon_sub.ndim == 2 and lat_sub.ndim == 2:
                lon_subb = lon_sub
                lat_subb = lat_sub
    
            if lon_sub.ndim == 1 and lat_sub.ndim == 1:
                lon_subb, lat_subb =  np.meshgrid(lon_sub, lat_sub)
    
            lon_orig = lon_subb
            lat_orig = lat_subb
    
            # Example data on original grid
            data = tmp_sub
    
            # --- Target grid ---
            # subset
            lon_targ, lat_targ =  np.meshgrid(lon_obs_sort, lat_obs)
    
            # --- Interpolate ---
            points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
            values = data.ravel()
            tmp_sub_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')
            
            bias_spatial[ex,y,f,:,:,] = tmp_sub_interp - skt_obs_sort
    
            #########################################################
            # Calculate differen between ARAFS and ERA5
            oksst1 = np.isfinite(skt_obs_sort[:,:])
            skt_obs_ff = skt_obs_sort[oksst1]
    
            tmp_sub_interp_ff = tmp_sub_interp[oksst1]
            oksst2 = np.isfinite(tmp_sub_interp_ff)
            tmp_sub_interp_fin = tmp_sub_interp_ff[oksst2]
            skt_obs_fin = skt_obs_ff[oksst2]
    
            number_obs[ex,y,f] = len(skt_obs_fin)
    
            mean_tmp_model[ex,y,f] = np.nanmean(tmp_sub_interp_fin)
            std_tmp_model[ex,y,f] = np.nanstd(tmp_sub_interp_fin)
    
            mean_skt_obs[ex,y,f] = np.nanmean(skt_obs_fin)
            std_skt_obs[ex,y,f] = np.nanstd(skt_obs_fin)
    
            bias[ex,y,f] = np.nanmean(tmp_sub_interp_fin)-np.nanmean(skt_obs_fin)
    
            corr[ex,y,f] = np.corrcoef(tmp_sub_interp_fin,skt_obs_fin)[0,1]
 
            #########################################################
            '''
            # Scatter plots
            fig = plt.figure(figsize=(6,5))
            plt.plot(skt_obs_fin,tmp_sub_interp_fin,'.',label='Total Obs used = '+str(number_obs[ex,y,f])+' \n Bias = '+str(np.round(bias[ex,y,f],2)))
            plt.plot(np.arange(-40,32),np.arange(-40,32),'-k')
            plt.legend()
            plt.title(ymdh + ' ' + fhhh)
            plt.xlabel('ERA5 Skin Temperature')
            plt.ylabel(exp)
            plt.axis('scaled')
 
            #########################################################
            # Spatial differences
            #levels = np.arange(-2,2.1,0.1)
            levels = np.arange(-10,10.1,1)
            plt.figure()
            plt.contourf(lon_obs_sort,lat_obs,tmp_sub_interp - skt_obs_sort,cmap='bwr',levels=levels,extend='both')
            plt.colorbar(shrink=0.7)
            plt.title(ymdh + ' ' + fhhh + '\n Bias (' + exp + ' - ERA5) = ' + str(np.round(bias[ex,y,f],2)),fontsize=16)
            plt.axis('scaled')
            
            # tmp model interpolated
            levels = np.arange(-60,61,10)
            plt.figure()
            plt.contourf(lon_obs_sort,lat_obs,tmp_sub_interp,cmap='jet',levels=levels,extend='both')
            plt.colorbar()
            plt.title(ymdh + ' ' + fhhh + ' ' + exp)
 
            levels = np.arange(-60,61,10)
            plt.figure()
            plt.contourf(lon_obs_sort,lat_obs,skt_obs_sort,cmap='jet',levels=levels,extend='both')
            plt.colorbar()
            plt.title(ymdh + ' ' + fhhh + ' ERA5')
            '''
 
        end = time.perf_counter()
        print(f"Elapsed time: {end - start:.6f} seconds")

np.savez('SKT_bias_ERA5_2023.npz', number_obs=number_obs, mean_tmp_model=mean_tmp_model, std_tmp_model=std_tmp_model, mean_skt_obs=mean_skt_obs, std_skt_obs=std_skt_obs, bias=bias, corr=corr, bias_spatial=bias_spatial, lon_obs=lon_obs_sort, lat_obs=lat_obs, exp_labels=exp_labels, ymdhs=ymdhs, fhhhs=fhhhs)
"""

loaded_data = np.load('SKT_bias_ERA5_2023.npz')
number_obs = loaded_data['number_obs']
mean_tmp_model = loaded_data['mean_tmp_model']
mean_skt_obs = loaded_data['mean_skt_obs']
std_tmp_model = loaded_data['std_tmp_model']
std_skt_obs = loaded_data['std_skt_obs']
bias = loaded_data['bias']
corr = loaded_data['corr']
bias_spatial = loaded_data['bias_spatial']
lon_obs = loaded_data['lon_obs']
lat_obs = loaded_data['lat_obs']
exp_labels = loaded_data['exp_labels']
ymdhs = loaded_data['ymdhs']
fhhhs = loaded_data['fhhhs']

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax = taylor_template(angle_lim,std_lim)
ax.plot(0,1,'o',label='ERA5',color='red',markersize=10,markeredgecolor='k')

for ex,exp in enumerate(exp_labels):  # loop the models
    print(exp)
    for y,ymdh in enumerate(ymdhs):
        print(ymdh)
        for f,fhhh in enumerate(conf['fhhhs']):
            print(fhhh)
            theta = np.arccos(corr[ex,y,f])
            rr = std_tmp_model[ex,y,f]/std_skt_obs[ex,y,f]
            print('corr = ',corr[ex,y,f])
            print('normalized std = ',rr)
            if f==0:
                #ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k',label=exp_labels[ex]+' =  '+str(np.round(number_obs[ex,f],2)))
                ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k',label=exp_labels[ex])
            else:
                ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k')

plt.legend(loc='upper left',bbox_to_anchor=[0.7,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
#plt.title('SST Cycle='+cycle+' fhour='+str(fhour),fontsize=18)
plt.title('Skin Temperature',fontsize=18)
plt.savefig('Taylor_diag_SKT',bbox_inches = 'tight',pad_inches = 0.1)

#y0 = 1.1
#plt.text(1.5,y0,'Bias',fontsize=14)
#for m in np.arange(len(folder_exps)):  # loop the models
#    plt.text(1.5,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias[m],3)),fontsize=14)

##################################################################
# time series mean skin temp
mean_tmp_model_all_cycles = np.mean(mean_tmp_model,axis=1)
mean_skt_obs_all_cycles = np.mean(mean_skt_obs,axis=1)
mean_bias = np.mean(np.mean(bias,axis=1),axis=1)
plt.figure(figsize=(11,6))
for ex,exp in enumerate(exp_labels):
    print(exp)
    plt.plot(fhhhs,mean_tmp_model_all_cycles[ex,:],'.-',markersize=15,markeredgecolor='k',label=exp+' TMP',color=exp_colors[ex])

plt.plot(fhhhs,mean_skt_obs_all_cycles[0,:],'.-',markersize=15,markeredgecolor='grey',label=conf['Data_prod']+' Skin Temp.',color='k')
plt.legend()
plt.text(0.05,0.1,'Bias ('+exp_labels[0]+') = '+str(np.round(mean_bias[0],2)),transform=ax.transAxes,fontsize=18)
plt.text(0.05,0.05,'Bias ('+exp_labels[1]+') = '+str(np.round(mean_bias[1],2)),transform=ax.transAxes,fontsize=18)
#plt.title('Domain Mean Temp: 30 cycles: Jan-Mar 2023-2026')
plt.title('Domain Mean Temp: '+ymdh)
plt.ylabel('Temperature ($^oC$)')
plt.xlabel('Forecast Hour')
plt.ylim([11,13.5])

##################################################################
# spatial map mean bias skin temp
#levels = np.arange(-10,10.1,1)
levels = np.arange(-5,5.1,0.5)
for ex,exp in enumerate(exp_labels):
    bias_spat = np.nanmean(np.nanmean(bias_spatial[ex,:,:,:,:],axis=0),axis=0)
    plt.figure()
    plt.contourf(lon_obs,lat_obs,bias_spat,cmap='bwr',levels=levels,extend='both')
    plt.colorbar(shrink=0.7)
    #plt.title(ymdh + ' ' + fhhh + '\n Bias (' + exp + ' - ERA5) = ' + str(np.round(bias[ex,y,f],2)),fontsize=16)
    plt.title(ymdh + '\n Bias (' + exp + ' - ERA5) = ' + str(np.round(np.nanmean(bias_spat),2)),fontsize=16)
    plt.axis('scaled')

