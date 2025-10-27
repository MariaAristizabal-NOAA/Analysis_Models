#!/usr/bin/env python3

"""This script is to compare ARAFS surface temperature with OSTIA SST"""

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import grib2io
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#from datetime import datetime, timedelta
import matplotlib.dates as mdates

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

# Parse the yaml config file
print('Parse the config file: plot_atmos3.yml:')
with open('plot_atmos3.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
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
sst_obs = np.asarray(OBS['analysed_sst']) - 273.15 # convert K to degC

timestamp_obs = time_obs.astype('int')/1e9

########################################################################
mean_sst_model = np.empty((len(experiments),len(conf['fhhhs'])))
mean_sst_model[:] = np.nan
std_sst_model = np.empty((len(experiments),len(conf['fhhhs'])))
std_sst_model[:] = np.nan

mean_sst_obs = np.empty((len(conf['fhhhs'])))
mean_sst_obs[:] = np.nan
std_sst_obs = np.empty((len(conf['fhhhs'])))
std_sst_obs[:] = np.nan

number_obs = np.empty((len(experiments),len(conf['fhhhs'])))
number_obs[:] = np.nan

bias = np.empty((len(experiments),len(conf['fhhhs'])))
bias[:] = np.nan
corr = np.empty((len(experiments),len(conf['fhhhs'])))
corr[:] = np.nan

for ex,exp in enumerate(experiments):
    print(exp)
    for f,fhhh in enumerate(conf['fhhhs']):
        print(fhhh)
        if exp == 'ARAFS' or exp == 'ARAFS_OCN':
            conf['fhours'] = int(fhhh[1:])
            conf['fcstTime'] = pd.to_timedelta(conf['fhours'], unit='h')
            conf['validTime'] = conf['initTime'] + conf['fcstTime']
            fname = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+fhhh+'.grb2'
            grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname)
            if os.path.exists(grib2file):
                print(f'grib2file: {grib2file}')
            else:
                fname = conf['ymdh']+'.'+conf['stormModel'].lower()+'.'+conf['stormDomain']+'.atm.'+fhhh+'.grb2'
                grib2file = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname)
                print(f'grib2file: {grib2file}')
            grb = grib2io.open(grib2file,mode='r')

            print('Extracting lat, lon')
            lat = grb.select(shortName='NLAT')[0].data
            lon = grb.select(shortName='ELON')[0].data

            print('Extracting surface temperature')
            levstr='surface'
            tmp = grb.select(shortName='TMP', level=levstr)[0].data
            tmp = tmp - 273.15 # convert K to degC

            # Transform to geographic coordinated
            lon[lon>180] = lon[lon>180] - 360

            # subset
            oklon = np.logical_and(lon[0,:]>=xlim[0],lon[0,:]<xlim[1])
            oklat = np.logical_and(lat[:,0]>=ylim[0],lat[:,0]<ylim[1])

            lon_sub = lon[oklat,:][:,oklon]
            lat_sub = lat[oklat,:][:,oklon]
            tmp_sub = tmp[oklat,:][:,oklon]

        ################################################################
        if exp == 'MOM6':
            # SST from ocean
            fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+fhhh+'.nc'
            ncfile = os.path.join(conf['COMmodels'][ex]+conf['ymdh']+'/'+conf['Stormid']+'/', fname_ocn)
            nc = xr.open_dataset(ncfile)

            tmp = np.asarray(nc['SST'][0,:,:])
            lon = np.asarray(nc.xh)
            lat = np.asarray(nc.yh)

            # subset
            oklon = np.logical_and(lon>=xlim[0],lon<xlim[1])
            oklat = np.logical_and(lat>=ylim[0],lat<ylim[1])

            lon_sub = lon[oklon]
            lat_sub = lat[oklat]
            tmp_sub = tmp[oklat,:][:,oklon]

        #################################################################
        # interpolation from fv3 and MOM6 to Ostia grid
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
        lon_targ, lat_targ =  np.meshgrid(lon_obs, lat_obs)

        # --- Interpolate ---
        points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
        values = data.ravel()
        tmp_sub_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')

        #################################################################
        # Choose correct timestamp from observations
        tt = conf['validTime'].date()
        timestamp_model = datetime.datetime.combine(tt, datetime.time()).timestamp()
        okt = np.where(timestamp_obs <  timestamp_model)[0][-1]

        print('Model time ',conf['validTime'])
        print('Obs time ',time_obs[okt])

        #################################################################
        # Calculate differen between tmp from ARAFS and Ostia
        oksst1 = np.isfinite(sst_obs[okt,:,:])
        sst_obs_ff = sst_obs[okt,oksst1]

        tmp_sub_interp_ff = tmp_sub_interp[oksst1]
        oksst2 = np.isfinite(tmp_sub_interp_ff)
        tmp_sub_interp_fin = tmp_sub_interp_ff[oksst2]
        sst_obs_fin = sst_obs_ff[oksst2]

        number_obs[ex,f] = len(sst_obs_fin)

        mean_sst_model[ex,f] = np.nanmean(tmp_sub_interp_fin)
        std_sst_model[ex,f] = np.nanstd(tmp_sub_interp_fin)

        mean_sst_obs[f] = np.nanmean(sst_obs_fin)
        std_sst_obs[f] = np.nanstd(sst_obs_fin)

        bias[ex,f] = np.nanmean(tmp_sub_interp_fin)-np.nanmean(sst_obs_fin)

        corr[ex,f] = np.corrcoef(tmp_sub_interp_fin,sst_obs_fin)[0,1]

##################################################################
# Taylor diagram
angle_lim = np.pi/2
std_lim = 1.5

fig,ax = taylor_template(angle_lim,std_lim)
ax.plot(0,1,'o',label='OSTIA',color='red',markersize=10,markeredgecolor='k')

for ex,exp in enumerate(experiments):  # loop the models
    print(exp)
    for f,fhhh in enumerate(conf['fhhhs']):
        print(fhhh)
        theta = np.arccos(corr[ex,f])
        rr = std_sst_model[ex,f]/std_sst_obs[f]
        print('corr = ',corr[ex,f])
        print('normalized std = ',rr)
        if f==0:
            ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k',label=experiments[ex]+' =  '+str(np.round(number_obs[ex,f],2)))
        else:
            ax.plot(theta,rr,exp_markers[ex],color = exp_colors[ex],markersize=8,markeredgecolor='k')

plt.legend(loc='upper left',bbox_to_anchor=[0.7,1.15])

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))

contours = ax.contour(ts,rs,rms,[0.5,1,1.5],colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)
#plt.title('SST Cycle='+cycle+' fhour='+str(fhour),fontsize=18)
plt.title('SST',fontsize=18)
plt.savefig('Taylor_diag_SST',bbox_inches = 'tight',pad_inches = 0.1)

#y0 = 1.1
#plt.text(1.5,y0,'Bias',fontsize=14)
#for m in np.arange(len(folder_exps)):  # loop the models
#    plt.text(1.5,y0-(m+1)*0.1,exp_labels[m]+' = '+str(np.round(bias[m],3)),fontsize=14)

