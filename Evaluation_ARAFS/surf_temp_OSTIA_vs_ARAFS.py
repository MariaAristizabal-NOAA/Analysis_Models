#!/usr/bin/env python3

"""This script is to compare ARAFS surface temperature with OSTIA SST"""

import os
import yaml
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import grib2io
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Parse the yaml config file
print('Parse the config file: plot_atmos3.yml:')
with open('plot_atmos3.yml', 'rt') as f:
    conf = yaml.safe_load(f)
conf['initTime'] = pd.to_datetime(conf['ymdh'], format='%Y%m%d%H', errors='coerce')
experiments = conf['Experiments']
data_prod = conf['Data_prod']
fhhhs = conf['fhhhs']
commodels = conf['COMmodels']
comobs = conf['COMobs']
fileobs = conf['Fileobs']

xlim = conf['xlim']
ylim = conf['ylim']

################################################################
# Read Surface temp from observations products
ncfile = comobs[0] + '/' + fileobs[0]
OBS = xr.open_dataset(ncfile)

time_obs = np.asarray(OBS['time']) 
lat_obs = np.asarray(OBS['latitude']) 
lon_obs = np.asarray(OBS['longitude']) 
sst_obs = np.asarray(OBS['analysed_sst']) - 273.15 # convert K to degC 

timestamp_obs = time_obs.astype('int')/1e9

################################################################
# Read Surface temp from experiment
mean_tmp_model = np.empty((len(experiments),len(conf['fhhhs'])))
mean_tmp_model[:] = np.nan
std_tmp_model = np.empty((len(experiments),len(conf['fhhhs'])))
std_tmp_model[:] = np.nan
mean_tmp_model_int = np.empty((len(experiments),len(conf['fhhhs'])))
mean_tmp_model_int[:] = np.nan
std_tmp_model_int = np.empty((len(experiments),len(conf['fhhhs'])))
std_tmp_model_int[:] = np.nan

mean_sst_obs = np.empty((len(conf['fhhhs'])))
mean_sst_obs[:] = np.nan
std_sst_obs = np.empty((len(conf['fhhhs'])))
std_sst_obs[:] = np.nan

mean_diff = np.empty((len(experiments),len(conf['fhhhs'])))
mean_diff[:] = np.nan
std_diff = np.empty((len(experiments),len(conf['fhhhs'])))
std_diff[:] = np.nan

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

            land = grb.select(shortName='LAND')[0].data
            land[land==1.0] = np.nan
            land[land==0.0] = 1.0
            
            # Transform to geographic coordinated
            lon[lon>180] = lon[lon>180] - 360
            
            # subset
            oklon = np.logical_and(lon[0,:]>=xlim[0],lon[0,:]<xlim[1]) 
            oklat = np.logical_and(lat[:,0]>=ylim[0],lat[:,0]<ylim[1]) 
            
            lon_sub = lon[oklat,:][:,oklon]
            lat_sub = lat[oklat,:][:,oklon]
            tmp_sub = tmp[oklat,:][:,oklon]
            land_sub = land[oklat,:][:,oklon]
            tmp_sub = tmp_sub*land_sub

            lon_arafs = lon_sub
            lat_arafs = lat_sub

            '''
            fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.'+fhhh+'.nc'
            ncfile = os.path.join(conf['COMmodels'][2]+conf['ymdh']+'/'+conf['Stormid']+'/', fname_ocn)
            nc = xr.open_dataset(ncfile)

            tmp_ocn = np.asarray(nc['SST'][0,:,:])
            lon_ocn = np.asarray(nc.xh)
            lat_ocn = np.asarray(nc.yh)
            ocean = np.isfinite(tmp_ocn)

            # subset
            oklon_ocn = np.logical_and(lon_ocn>=xlim[0],lon_ocn<xlim[1])
            oklat_ocn = np.logical_and(lat_ocn>=ylim[0],lat_ocn<ylim[1])

            lon_sub_ocn = lon_ocn[oklon_ocn]
            lat_sub_ocn = lat_ocn[oklat_ocn]

            ocean_sub = ocean[oklat_ocn,:][:,oklon_ocn]

            lon_mom6 = lon_sub_ocn
            lat_mom6 = lat_sub_ocn
            '''
        
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
        # interpolation from MOM6 sub to ARAFS sub in order to use ocean masking on ARAFS fields and get rid off lakes in the mainland.
        # --- Original grid ---

        fname_ocn = '00e.'+conf['ymdh']+'.'+conf['stormModel'].lower()+'.mom6.f000.nc'
        ncfile = os.path.join(conf['COMmodels'][2]+conf['ymdh']+'/'+conf['Stormid']+'/', fname_ocn)
        nc = xr.open_dataset(ncfile)
 
        tmp_ocn = np.asarray(nc['SST'][0,:,:])
        lon_ocn = np.asarray(nc.xh)
        lat_ocn = np.asarray(nc.yh)
        ocean = np.isfinite(tmp_ocn)
 
        # subset
        oklon_ocn = np.logical_and(lon_ocn>=xlim[0],lon_ocn<xlim[1])
        oklat_ocn = np.logical_and(lat_ocn>=ylim[0],lat_ocn<ylim[1])
 
        lon_sub_ocn = lon_ocn[oklon_ocn]
        lat_sub_ocn = lat_ocn[oklat_ocn]
 
        ocean_sub = ocean[oklat_ocn,:][:,oklon_ocn]
 
        lon_mom6 = lon_sub_ocn
        lat_mom6 = lat_sub_ocn

        lon_subb, lat_subb =  np.meshgrid(lon_mom6, lat_mom6)

        lon_orig = lon_subb
        lat_orig = lat_subb

        # Example data on original grid
        data = ocean_sub.astype(int)

        # --- Target grid ---
        lon_targ = lon_arafs
        lat_targ = lat_arafs

        # --- Interpolate ---
        points = np.array([lon_orig.ravel(), lat_orig.ravel()]).T
        values = data.ravel()
        ocean_sub_interp = griddata(points, values, (lon_targ, lat_targ), method='linear')
        ocean_sub_interp[ocean_sub_interp<0.99] = np.nan

        #################################################################
        # Choose correct timestamp from observations
        tt = conf['validTime'].date()
        timestamp_model = datetime.datetime.combine(tt, datetime.time()).timestamp()
        okt = np.where(timestamp_obs <  timestamp_model)[0][-1]

        print('Model time ',conf['validTime'])
        print('Obs time ',time_obs[okt])

        #################################################################
        # Calculate differen between tmp from ARAFS and Ostia
        mean_tmp_model_int[ex,f] = np.nanmean(tmp_sub_interp)
        std_tmp_model_int[ex,f] = np.nanstd(tmp_sub_interp)
        if exp == 'ARAFS' or exp == 'ARAFS_OCN':
            mean_tmp_model[ex,f] = np.nanmean(tmp_sub*ocean_sub_interp)
            std_tmp_model[ex,f] = np.nanstd(tmp_sub*ocean_sub_interp)
        if exp == 'MOM6':
            mean_tmp_model[ex,f] = np.nanmean(tmp_sub)
            std_tmp_model[ex,f] = np.nanstd(tmp_sub)

        mean_sst_obs[f] = np.nanmean(sst_obs[okt,:,:])
        std_sst_obs[f] = np.nanstd(sst_obs[okt,:,:])

        diff = tmp_sub_interp-sst_obs[okt,:,:]
        mean_diff[ex,f] = np.nanmean(diff)
        std_diff[ex,f] = np.nanstd(diff)

        #################################################################
        '''
        plt.figure(figsize=(12,8))
        cflevels = np.arange(8,29,2)
        cmap = plt.get_cmap('turbo', len(cflevels)-1)
        norm = mpl.colors.BoundaryNorm(cflevels, cmap.N)
        cf = plt.pcolor(lon_sub,lat_sub,tmp_sub*ocean_sub_interp, cmap=cmap,norm=norm)
        #cf = plt.pcolor(lon_sub,lat_sub,tmp_sub, cmap=cmap,norm=norm)
        #cmap = plt.get_cmap('turbo')
        #cf = plt.contourf(lon_sub,lat_sub,tmp_sub, levels=cflevels, cmap=cmap, extend='both')
        plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
        lb = plt.contour(lon_sub,lat_sub,tmp_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
        plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
        plt.axis('scaled')
    
        model_info = os.environ.get('TITLEgraph','').strip()
        var_info = 'Surface Temperature (${^o}$C, shaded)'
        title_left = """{0}{1}""".format(model_info,var_info)
        plt.title(title_left, loc='left', y=0.99)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        plt.title(title_right, loc='right', y=0.99)
    
        plt.figure(figsize=(12,8))
        cflevels = np.arange(8,29,2)
        cmap = plt.get_cmap('turbo')
        cf = plt.contourf(lon_obs,lat_obs,tmp_sub_interp, levels=cflevels, cmap=cmap, extend='both')
        plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
        lb = plt.contour(lon_sub,lat_sub,tmp_sub,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
        plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
        plt.axis('scaled')
    
        model_info = os.environ.get('TITLEgraph','').strip()
        var_info = 'Surface Temperature (${^o}$C, shaded)'
        title_left = """{0}{1}""".format(model_info,var_info)
        plt.title(title_left, loc='left', y=0.99)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        plt.title(title_right, loc='right', y=0.99)
    
        plt.figure(figsize=(12,8))
        cflevels = np.arange(8,29,2)
        cmap = plt.get_cmap('turbo')
        cf = plt.contourf(lon_obs,lat_obs,sst_obs[okt,:,:], levels=cflevels, cmap=cmap, extend='both')
        plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
        lb = plt.contour(lon_obs,lat_obs,sst_obs[okt,:,:],levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
        plt.clabel(lb, lb.levels, inline=True,fmt='%1.0f', fontsize=6,colors='grey')
        plt.axis('scaled')
    
        model_info = 'OSTIA '
        var_info = 'Surface Temperature (${^o}$C, shaded)'
        title_left = """{0}{1}""".format(model_info,var_info)
        plt.title(title_left, loc='left', y=0.99)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        plt.title(title_right, loc='right', y=0.99)
        plt.xlim([xlim[0],xlim[1]])
        '''
     
        '''
        plt.figure(figsize=(12,8))
        cflevels = np.arange(-0.3,0.31,0.05)
        cmap = plt.get_cmap('bwr')
        cf = plt.contourf(lon_obs,lat_obs,diff, levels=cflevels, cmap=cmap, extend='both')
        plt.colorbar(cf, orientation='vertical', pad=0.02, aspect=30, shrink=0.7, extendrect=True, ticks=cflevels)
        plt.axis('scaled')
        plt.xlim([xlim[0],xlim[1]])
    
        model_info = exp + ' - OSTIA '
        var_info = 'Surf Temp Diff (${^o}$C)'
        title_left = """{0}{1}""".format(model_info,var_info)
        plt.title(title_left, loc='left', y=0.99)
        title_right = conf['initTime'].strftime('Init: %Y%m%d%HZ ')+fhhh.upper()+conf['validTime'].strftime(' Valid: %Y%m%d%HZ')
        plt.title(title_right, loc='right', y=0.99)
        plt.savefig(exp+'_'+fhhh)
        '''
    #################################################################
plt.figure(figsize=(10,6))
for ex,exp in enumerate(experiments):
    plt.errorbar(fhhhs,mean_tmp_model[ex,:],std_tmp_model[ex,:],fmt='o-',capsize=5,label=exp)
#plt.errorbar(fhhhs,mean_sst_obs,std_sst_obs,fmt='o-',capsize=5,label='OSTIA')
plt.legend()
plt.ylabel('SST no interpolation ($^oC$)')
plt.xlabel('Forecast Hour')
plt.savefig('SST_mean')

plt.figure(figsize=(10,6))
for ex,exp in enumerate(experiments):
    plt.errorbar(fhhhs,mean_tmp_model_int[ex,:],std_tmp_model_int[ex,:],fmt='o-',capsize=5,label=exp)
plt.errorbar(fhhhs,mean_sst_obs,std_sst_obs,fmt='o-',capsize=5,label='OSTIA')
plt.legend()
plt.ylabel('SST interpolation onto OSTIA ($^oC$)')
plt.xlabel('Forecast Hour')
plt.savefig('SST_mean_interp')

plt.figure(figsize=(10,6))
for ex,exp in enumerate(experiments):
    plt.errorbar(fhhhs,mean_diff[ex,:],std_diff[ex,:],fmt='o-',capsize=5,label=exp)
    plt.legend()
plt.plot(fhhhs,np.zeros(len(conf['fhhhs'])),color='k')
plt.ylabel('SST diff ($^oC$)')
plt.xlabel('Forecast Hour')
plt.savefig('SST_mean_diff')


