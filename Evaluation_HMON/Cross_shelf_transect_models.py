#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 13:57:56 2020

@author: aristizabal
"""

#%%
#storm_id = '95l'
#cycle = '2020081018'

#%%
def HWRF_POM_cross_transect(storm_id,cycle,tempt_lim_MAB,tempt_lim_GoM):

    #%%
    from matplotlib import pyplot as plt
    import numpy as np
    import xarray as xr
    from matplotlib.dates import date2num, num2date
    from datetime import datetime, timedelta
    import os
    import os.path
    import glob
    import cmocean
    
    # Increase fontsize of labels globally
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('legend',fontsize=14)
    
    #%%
    
    ti = datetime.today()
    ffig = '/home/Maria.Aristizabal/Figures/'+ str(ti.year) + '/' + ti.strftime('%b-%d') 
    folder_fig =  ffig + '/' + storm_id + '_' + cycle + '/'
    
    os.system('mkdir ' +  ffig)
    os.system('mkdir ' +  folder_fig)
    
    
    #%% Bathymetry file
    
    bath_file = '/scratch2/NOS/nosofs/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
    
    #%% Reading bathymetry data
    ncbath = xr.open_dataset(bath_file)
    bath_lat = ncbath.variables['lat'][:]
    bath_lon = ncbath.variables['lon'][:]
    bath_elev = ncbath.variables['elevation'][:]
    
    #%% folder and file names
    ti = datetime.today() - timedelta(1)
    
    folder_hwrf_pom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF_POM_' + storm_id + '_' + str(ti.year) + '/' + 'HWRF_POM_' + storm_id + '_' + cycle + '/'
    
    #%% Reading POM grid file
    grid_file = sorted(glob.glob(os.path.join(folder_hwrf_pom,'*grid*.nc')))[0]
    pom_grid = xr.open_dataset(grid_file)
    lon_pom = np.asarray(pom_grid['east_e'][:])
    lat_pom = np.asarray( pom_grid['north_e'][:])
    zlevc = np.asarray(pom_grid['zz'][:])
    topoz = np.asarray(pom_grid['h'][:])
    
    #%% Getting list of POM files
    ncfiles = sorted(glob.glob(os.path.join(folder_hwrf_pom,'*pom.0*.nc')))
    
    # Reading POM time
    time_pom = []
    for i,file in enumerate(ncfiles):
        print(i)
        pom = xr.open_dataset(file)
        tpom = pom['time'][:]
        timestamp_pom = date2num(tpom)[0]
        time_pom.append(num2date(timestamp_pom))
    
    time_POM = np.asarray(time_pom)
    
    oktime_POM = np.where(time_POM == time_POM[1])[0][0] #first file
    
    #%% Figure temp transect along Endurance line
    
    x1 = -74.1
    y1 = 39.4
    x2 = -73.0
    y2 = 38.6
    # Slope
    m = (y1-y2)/(x1-x2)
    # Intercept
    b = y1 - m*x1
    
    X = np.arange(x1,-72,0.05)
    Y = b + m*X
    
    dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
    
    oklon = np.round(np.interp(X,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
    oklat = np.round(np.interp(Y,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)
    topoz_pom = np.asarray(topoz[oklat,oklon])
    zmatrix_POM = np.dot(topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
    dist_matrix = np.tile(dist,(zmatrix_POM.shape[0],1))
    
    trans_temp_POM = np.empty((zmatrix_POM.shape[0],zmatrix_POM.shape[1]))
    trans_temp_POM[:] = np.nan
    pom = xr.open_dataset(ncfiles[oktime_POM])
    for x in np.arange(len(X)):
        print(x)
        trans_temp_POM[:,x] = np.asarray(pom['t'][0,:,oklat[x],oklon[x]])
    
    #min_valt = 4
    #max_valt = 27
    nlevelst = int(np.ceil(tempt_lim_MAB[1] - tempt_lim_MAB[0] + 1))
    kw = dict(levels = np.linspace(tempt_lim_MAB[0],tempt_lim_MAB[1],nlevelst))
    
    fig, ax = plt.subplots(figsize=(9, 3))
    plt.contourf(dist_matrix,zmatrix_POM,trans_temp_POM,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(dist_matrix,zmatrix_POM,trans_temp_POM,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.title('HWRF-POM Endurance Line MAB ' + 'Storm ' + storm_id + ' Cycle ' + cycle + ' \n on ' + str(time_POM[oktime_POM])[0:13],fontsize=16)
    plt.ylim([-100,0])
    plt.xlim([0,200])
    ax.set_ylabel('Depth (m)',fontsize=14)
    ax.set_xlabel('Along Transect Distance (km)',fontsize=14)
    
    file = folder_fig + 'HWRF_POM_temp_MAB_endurance_line_cycle_'+cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Bathymetry GEBCO HYCOM domain
    kw = dict(levels =  np.arange(-5000,1,200))
    
    plt.figure()
    plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
    plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
    plt.plot(X,Y,'-k')
    plt.colorbar()
    plt.axis('scaled')
    plt.title('GEBCO Bathymetry')
    plt.xlim(-76,-70)
    plt.ylim(35,42)
    
    file = folder_fig + 'MAB_transect'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Figure temp transect across GoMex
    
    #x1 = -90
    #y1 = 20 + 52/60
    #x2 = -90
    #y2 = 30
    #Y = np.arange(y1,y2,0.05)
    #X = np.tile(x1,len(Y))

    x1 = -89.5
    y1 = 28.8
    x2 = -85.1
    y2 = 26.7

    # Slope
    m = (y1-y2)/(x1-x2)
    # Intercept
    b = y1 - m*x1

    X = np.arange(x1,x2,0.05)
    Y = b + m*X

    dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
    
    oklon = np.round(np.interp(X,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
    oklat = np.round(np.interp(Y,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)
    topoz_pom = np.asarray(topoz[oklat,oklon])
    zmatrix_POM = np.dot(topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
    dist_matrix = np.tile(dist,(zmatrix_POM.shape[0],1))
    
    trans_temp_POM = np.empty((zmatrix_POM.shape[0],zmatrix_POM.shape[1]))
    trans_temp_POM[:] = np.nan
    pom = xr.open_dataset(ncfiles[oktime_POM])
    for x in np.arange(len(X)):
        print(x)
        trans_temp_POM[:,x] = np.asarray(pom['t'][0,:,oklat[x],oklon[x]])
    
    #min_valt = 12
    #max_valt = 32
    #kw = dict(levels = np.arange(tempt_lim_GoM[0],tempt_lim_GoM[1],1))
    kw = dict(levels = np.arange(tempt_lim_GoM[0],tempt_lim_GoM[1],1))
    
    #fig, ax = plt.subplots(figsize=(9, 3))
    fig, ax = plt.subplots()
    plt.contourf(dist_matrix,zmatrix_POM,trans_temp_POM,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(dist_matrix,zmatrix_POM,trans_temp_POM,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.title('HWRF-POM  Across GoMex ' + 'Storm ' + storm_id + ' Cycle ' + ' \n on ' + str(time_POM[oktime_POM])[0:13],fontsize=16)
    plt.ylim([-300,0])
    #plt.xlim([0,200])
    ax.set_ylabel('Depth (m)',fontsize=14)
    ax.set_xlabel('Along Transect Distance (km)',fontsize=14)
    
    file = folder_fig + 'HWRF_POM_temp_GoMex_across_cycle_'+cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Bathymetry GEBCO HYCOM domain
    kw = dict(levels =  np.arange(-5000,1,200))
    
    plt.figure()
    plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
    plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
    plt.plot(X,Y,'-k')
    plt.colorbar()
    plt.axis('scaled')
    plt.title('GEBCO Bathymetry')
    plt.xlim(-98,-80)
    plt.ylim(18,32)
    
    file = folder_fig + 'GoMex_transect'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
#%%
def HMON_HYCOM_cross_transect(storm_id,cycle,tempt_lim_MAB,tempt_lim_GoM):    
    
    #%% 
    import numpy as np
    import matplotlib.pyplot as plt
    from datetime import datetime, timedelta
    import os
    import os.path
    import glob
    import cmocean
    from matplotlib.dates import date2num, num2date
    import xarray as xr
    
    import sys
    sys.path.append('/home/Maria.Aristizabal/NCEP_scripts/')
    from utils4HYCOM import readgrids
    #from utils4HYCOM import readdepth, readVar
    from utils4HYCOM2 import readBinz
    
    # Increase fontsize of labels globally
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('legend',fontsize=14)
    
    #%%
    
    ti = datetime.today()
    ffig = '/home/Maria.Aristizabal/Figures/'+ str(ti.year) + '/' + ti.strftime('%b-%d') 
    folder_fig =  ffig + '/' + storm_id + '_' + cycle + '/'
    
    os.system('mkdir ' +  ffig)
    os.system('mkdir ' +  folder_fig)
    
    #%% Bathymetry file
    
    bath_file = '/scratch2/NOS/nosofs/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
    
    #%% Reading bathymetry data
    ncbath = xr.open_dataset(bath_file)
    bath_lat = ncbath.variables['lat'][:]
    bath_lon = ncbath.variables['lon'][:]
    bath_elev = ncbath.variables['elevation'][:]
    
        #%% folder and file names    
    ti = datetime.today() - timedelta(1)
    
    folder_hmon_hycom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HMON_HYCOM_' + storm_id + '_' + str(ti.year) + '/' + 'HMON_HYCOM_' + storm_id + '_' + cycle + '/'
    
    #%% Reading RTOFS grid    
    grid_file = sorted(glob.glob(os.path.join(folder_hmon_hycom,'*regional.grid.*')))[0][:-2]
    
    #%% Reading RTOFS grid
    print('Retrieving coordinates from RTOFS')
    # Reading lat and lon
    #lines_grid = [line.rstrip() for line in open(grid_file+'.b')]
    lon_hycom = np.array(readgrids(grid_file,'plon:',[0]))
    lat_hycom = np.array(readgrids(grid_file,'plat:',[0]))
    
    #depth_HMON_HYCOM = np.asarray(readdepth(HMON_HYCOM_depth,'depth'))
    
    # Reading depths
    afiles = sorted(glob.glob(os.path.join(folder_hmon_hycom,'*hat10_3z'+'*.a')))
    lines=[line.rstrip() for line in open(afiles[0][:-2]+'.b')]
    z = []
    for line in lines[6:]:
        if line.split()[2]=='temp':
            #print(line.split()[1])
            z.append(float(line.split()[1]))
    depth_HYCOM = np.asarray(z)
    
    time_HYCOM = []
    for x, file in enumerate(afiles):
        print(x)
        #lines=[line.rstrip() for line in open(file[:-2]+'.b')]
    
        #Reading time stamp
        year = int(file.split('/')[-1].split('.')[1][0:4])
        month = int(file.split('/')[-1].split('.')[1][4:6])
        day = int(file.split('/')[-1].split('.')[1][6:8])
        hour = int(file.split('/')[-1].split('.')[1][8:10])
        dt = int(file.split('/')[-1].split('.')[-2][1:])
        timestamp_HYCOM = date2num(datetime(year,month,day,hour)) + dt/24
        time_HYCOM.append(num2date(timestamp_HYCOM))
    
    # Reading 3D variable from binary file
    oktime = 0 # first file 
    temp_HMON_HYCOM = readBinz(afiles[oktime][:-2],'3z','temp')
    #salt_HMON_HYCOM = readBinz(afiles[oktime][:-2],'3z','salinity')
    
    #%%
    
    x1 = -74.1
    y1 = 39.4
    x2 = -73.0
    y2 = 38.6
    # Slope
    m = (y1-y2)/(x1-x2)
    # Intercept
    b = y1 - m*x1
    
    X = np.arange(x1,-72,0.05)
    Y = b + m*X
    
    dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
    
    oklon = np.round(np.interp(X,lon_hycom[0,:]-360,np.arange(len(lon_hycom[0,:])))).astype(int)
    oklat = np.round(np.interp(Y,lat_hycom[:,0],np.arange(len(lat_hycom[:,0])))).astype(int)
    
    trans_temp_HYCOM = temp_HMON_HYCOM[oklat,oklon,:]
    
    #min_valt = 4
    #max_valt = 27
    nlevelst = tempt_lim_MAB[1] - tempt_lim_MAB[0] + 1
    kw = dict(levels = np.linspace(tempt_lim_MAB[0],tempt_lim_MAB[1],nlevelst))
    
    fig, ax = plt.subplots(figsize=(9, 3))
    plt.contourf(dist,-depth_HYCOM,trans_temp_HYCOM.T,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(dist,-depth_HYCOM,trans_temp_HYCOM.T,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Along Transect Distance (km)',fontsize=14)
    plt.title('HMON-HYCOM Endurance Line ' + 'Storm ' + storm_id + ' Cycle ' + cycle + ' \n on ' + str(time_HYCOM[oktime])[0:13],fontsize=16)
    plt.ylim([-100,0])
    plt.xlim([0,200])
    
    file = folder_fig + 'HMON_HYCOM_temp_MAB_endurance_line_cycle_'+cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Bathymetry GEBCO HYCOM domain
    kw = dict(levels =  np.arange(-5000,1,200))
    
    plt.figure()
    plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
    plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
    plt.plot(X,Y,'-k')
    plt.colorbar()
    plt.axis('scaled')
    plt.title('GEBCO Bathymetry')
    plt.xlim(-76,-70)
    plt.ylim(35,42)
    
    file = folder_fig + 'MAB_transect'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%%
    
    #x1 = -90
    #y1 = 20 + 52/60
    #x2 = -90
    #y2 = 30
    #Y = np.arange(y1,y2,0.05)
    #X = np.tile(x1,len(Y))

    x1 = -89.5
    y1 = 28.8
    x2 = -85.1
    y2 = 26.7

    # Slope
    m = (y1-y2)/(x1-x2)
    # Intercept
    b = y1 - m*x1

    X = np.arange(x1,x2,0.05)
    Y = b + m*X

    
    dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
    
    oklon = np.round(np.interp(X,lon_hycom[0,:]-360,np.arange(len(lon_hycom[0,:])))).astype(int)
    oklat = np.round(np.interp(Y,lat_hycom[:,0],np.arange(len(lat_hycom[:,0])))).astype(int)
    
    trans_temp_HYCOM = temp_HMON_HYCOM[oklat,oklon,:]
    
    #min_valt = 12
    #max_valt = 32
    kw = dict(levels = np.arange(tempt_lim_GoM[0],tempt_lim_GoM[1]))
    
    #fig, ax = plt.subplots(figsize=(9, 3))
    fig, ax = plt.subplots()
    plt.contourf(dist,-depth_HYCOM,trans_temp_HYCOM.T,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(dist,-depth_HYCOM,trans_temp_HYCOM.T,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Along Transect Distance (km)',fontsize=14)
    plt.title('HMON-HYCOM Across GoMex ' + 'Storm ' + storm_id + ' Cycle ' + cycle + ' \n on ' + str(time_HYCOM[oktime])[0:13],fontsize=16)
    plt.ylim([-300,0])
    #plt.xlim([0,200])
    
    file = folder_fig + 'HMON_HYCOM_temp_GoMex_across_cycle_'+cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Bathymetry GEBCO HYCOM domain
    kw = dict(levels =  np.arange(-5000,1,200))
    
    plt.figure()
    plt.contour(bath_lon,bath_lat,bath_elev,levels=[0],colors='k')
    plt.contourf(bath_lon,bath_lat,bath_elev,cmap=cmocean.cm.topo,**kw)
    plt.plot(X,Y,'-k')
    plt.colorbar()
    plt.axis('scaled')
    plt.title('GEBCO Bathymetry')
    plt.xlim(-98,-80)
    plt.ylim(18,32)
    
    file = folder_fig + 'GoMex_transect'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
