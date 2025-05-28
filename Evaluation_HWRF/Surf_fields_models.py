
"""
Created on Thu Jun 25 14:08:36 2020

@author: aristizabal
"""

def HWRF_POM_fields(cycle,storm_id,lon_forec_track,lat_forec_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig):

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
    
    #%% Bathymetry file
    
    bath_file = '/scratch2/NOS/nosofs/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
    
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
    
    #%% folder and file names
    ti = datetime.today() - timedelta(1)
    
    #folder_hwrf_pom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF_POM_' + storm_id + '_' + str(ti.year) + '/' + 'HWRF_POM_' + storm_id + '_' + cycle + '/'
    folder_hwrf_pom = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF_POM_' + storm_id + '_' + str(ti.year) + '/' + 'HWRF_POM_' + storm_id + '_' + cycle + '/'
    
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
    
    #%%
    oklon_POM = np.where(np.logical_and(lon_pom[0,:] >= lon_lim[0],lon_pom[0,:] <= lon_lim[1]))[0]
    oklat_POM = np.where(np.logical_and(lat_pom[:,0] >= lat_lim[0],lat_pom[:,0] <= lat_lim[1]))[0]
    #oktime_POM = np.where(time_POM == time_POM[0])[0][0]
    oktime_POM = np.where(time_POM == time_POM[1])[0][0] #second file
    
    lon_POM =lon_pom[0,oklon_POM]
    lat_POM =lat_pom[oklat_POM,0]
    
    #%% Getting POM variables
    #pom = xr.open_dataset(ncfiles[0]) #firts file is cycle
    pom = xr.open_dataset(ncfiles[oktime_POM]) #second file in cycle
    sst_POM = np.asarray(pom['t'][0,0,oklat_POM,oklon_POM])
    sst_POM[sst_POM==0] = np.nan
    sss_POM = np.asarray(pom['s'][0,0,oklat_POM,oklon_POM])
    sss_POM[sss_POM==0] = np.nan
    ssh_POM = np.asarray(pom['elb'][0,oklat_POM,oklon_POM])
    ssh_POM[ssh_POM==0] = np.nan
    su_POM = np.asarray(pom['u'][0,0,oklat_POM,oklon_POM])
    su_POM[su_POM==0] = np.nan
    sv_POM = np.asarray(pom['v'][0,0,oklat_POM,oklon_POM])
    sv_POM[sv_POM==0] = np.nan
    
    #%% SST
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5),vmin=temp_lim[0],vmax=temp_lim[1])
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)
    plt.contour(lon_POM,lat_POM,sst_POM,levels=[26],colors='white')
    plt.contourf(lon_POM,lat_POM,sst_POM,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HWRF-POM SST ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_POM[oktime_POM])[0:13],fontsize=16)
    
    file_name = folder_fig + 'HWRF_POM_SST_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

   #%% SST and velocity
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)
    plt.contourf(lon_POM,lat_POM,sst_POM,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HWRF-POM SST ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_POM[oktime_POM])[0:13],fontsize=16)

    q = plt.quiver(lon_POM[::7], lat_POM[::7],su_POM[::7,::7],sv_POM[::7,::7])
    plt.quiverkey(q,np.max(lon_POM)+5,np.max(lat_POM)+0.2,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})

    file_name = folder_fig + 'HWRF_POM_SST_UV_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% SSS
    kw = dict(levels = np.arange(salt_lim[0],salt_lim[1],0.5))    
 
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_POM,lat_POM,sss_POM,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HWRF-POM SSS ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_POM[oktime_POM])[0:13],fontsize=16)
        
    file_name = folder_fig + 'HWRF_POM_SSS_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% SSH    
    kw = dict(levels = np.arange(-1,1.1,0.1))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_POM,lat_POM,ssh_POM,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HWRF-POM SSH ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_POM[oktime_POM])[0:13],fontsize=16)

    file_name = folder_fig + 'HWRF_POM_SSH_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% SSH   
    kw = dict(levels = np.arange(-1,1.1,0.1))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,[0,10000],colors='seashell')
    plt.contourf(lon_POM,lat_POM,ssh_POM,cmap=cmocean.cm.curl,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('meters',fontsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HWRF-POM SSH ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_POM[oktime_POM])[0:13],fontsize=16)
        
    q = plt.quiver(lon_POM[::7], lat_POM[::7],su_POM[::7,::7],sv_POM[::7,::7])
    plt.quiverkey(q,np.max(lon_POM)+5,np.max(lat_POM)+0.2,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file_name = folder_fig + 'HWRF_POM_SSH_UV_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Figure temp transect along storm path 
    
    lat_forec_tracku, ind = np.unique(lat_forec_track,return_index=True)
    lon_forec_tracku = lon_forec_track[ind]
    lon_forec_track_interp = np.interp(lat_pom[:,0],lat_forec_tracku,lon_forec_tracku,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_pom[:,0])
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
    oklon = np.round(np.interp(lon_forec_track_int,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)
    topoz_pom = np.asarray(topoz[oklat,oklon])
    zmatrix_POM = np.dot(topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
    lat_matrix = np.tile(lat_pom[oklat,0],(zmatrix_POM.shape[0],1))
    
    trans_temp_POM = np.empty((zmatrix_POM.shape[0],zmatrix_POM.shape[1]))
    trans_temp_POM[:] = np.nan
    for x in np.arange(len(lon_forec_track_int)):
        trans_temp_POM[:,x] = np.asarray(pom['t'][0,:,oklat[x],oklon[x]])
    
    kw = dict(levels = np.arange(tempt_lim[0],tempt_lim[1],1))    

    plt.figure()
    plt.contourf(lat_matrix,zmatrix_POM,trans_temp_POM,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(lat_matrix,zmatrix_POM,trans_temp_POM,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Latitude ($^o$)',fontsize=14)
    plt.title('HWRF-POM Temperature ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n along Forecasted Storm Track on '+str(time_POM[oktime_POM])[0:13],fontsize=12)
    plt.ylim([-300,0])
    
    file = folder_fig + 'HWRF-POM_temp_along_forecasted_track_' + cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% Figure temp transect along chosen transect
'''    
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
    
    X = lon_forec_track
    Y = lat_forec_track
    
    dist = np.sqrt((X-x1)**2 + (Y-y1)**2)*111 # approx along transect distance in km
    
    oklon = np.round(np.interp(X,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
    oklat = np.round(np.interp(Y,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)
    topoz_pom = np.asarray(topoz[oklat,oklon])
    zmatrix_POM = np.dot(topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
    dist_matrix = np.tile(dist,(zmatrix_POM.shape[0],1))
    
    trans_temp_POM = np.empty((zmatrix_POM.shape[0],zmatrix_POM.shape[1]))
    trans_temp_POM[:] = np.nan
    for x in np.arange(len(lon_forec_track)):
        trans_temp_POM[:,x] = np.asarray(pom['t'][oktime_POM,:,oklat[x],oklon[x]])
       
    max_valt = 26
    min_valt = 8  
    nlevelst = max_valt - min_valt + 1
    
    kw = dict(levels = np.linspace(min_valt,max_valt,nlevelst))
    
    fig, ax = plt.subplots(figsize=(9, 3))
    plt.contourf(dist_matrix,zmatrix_POM,trans_temp_POM,cmap='gnuplot',**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(dist_matrix,zmatrix_POM,trans_temp_POM,[12],color='k')
    cbar.ax.set_ylabel('Temperature ($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.title('HWRF-POM Transect on ' + str(time_POM[oktime_POM])[0:13],fontsize=16)
    plt.ylim([-100,0])
    plt.xlim([0,200])
    ax.set_ylabel('Depth (m)',fontsize=14)
    ax.set_xlabel('Along Transect Distance (km)',fontsize=14)
    
    file = folder_fig + 'HWRF-POM_cross_shelf_transect'
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
'''
    
#%%    

def HMON_HYCOM_fields(cycle,storm_id,lon_forec_track,lat_forec_track,lon_lim,lat_lim,temp_lim,salt_lim,temp200_lim,salt200_lim,tempb_lim,tempt_lim,folder_fig):
        
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
    
    #%% Bathymetry file
    
    bath_file = '/scratch2/NOS/nosofs/Maria.Aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
    
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
        
    #%% folder and file names    
    ti = datetime.today() - timedelta(1)
    
    #folder_hmon_hycom = '/scratch2/NOS/nosofs/Maria.Aristizabal/HMON_HYCOM_' + storm_id + '_' + str(ti.year) + '/' + 'HMON_HYCOM_' + storm_id + '_' + cycle + '/'
    folder_hmon_hycom = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HMON_HYCOM_' + storm_id + '_' + str(ti.year) + '/' + 'HMON_HYCOM_' + storm_id + '_' + cycle + '/'
    
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

    #%%
    
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
    salt_HMON_HYCOM = readBinz(afiles[oktime][:-2],'3z','salinity')
    uvel_HMON_HYCOM = readBinz(afiles[oktime][:-2],'3z','u-veloc.')
    vvel_HMON_HYCOM = readBinz(afiles[oktime][:-2],'3z','v-veloc.')
    
    #%% 
    oklon_HYCOM = np.where(np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360))[0]
    oklat_HYCOM = np.where(np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1]))[0]
    
    lon_HYCOM =lon_hycom[0,oklon_HYCOM]-360
    lat_HYCOM =lat_hycom[oklat_HYCOM,0]
    
    sst_HYCOM = temp_HMON_HYCOM[oklat_HYCOM,:,0][:,oklon_HYCOM]
    sss_HYCOM = salt_HMON_HYCOM[oklat_HYCOM,:,0][:,oklon_HYCOM]
    su_HYCOM = uvel_HMON_HYCOM[oklat_HYCOM,:,0][:,oklon_HYCOM]
    sv_HYCOM = vvel_HMON_HYCOM[oklat_HYCOM,:,0][:,oklon_HYCOM]
    
    #okdepth200 = np.where(depth_HYCOM >= 200)[0][0]    
    #temp200_HYCOM = temp_HMON_HYCOM[oklat_HYCOM,:,okdepth200][:,oklon_HYCOM]
    #salt200_HYCOM = salt_HMON_HYCOM[oklat_HYCOM,:,okdepth200][:,oklon_HYCOM]
    
    #%% SST    
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))    

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)    
    plt.contour(lon_HYCOM,lat_HYCOM,sst_HYCOM,levels=[26],colors='white')
    plt.contourf(lon_HYCOM,lat_HYCOM,sst_HYCOM,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HMON-HYCOM SST ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_HYCOM[oktime])[0:13],fontsize=16)
    
    file_name = folder_fig + 'HMON_HYCOM_SST_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% SSS
    kw = dict(levels = np.arange(salt_lim[0],salt_lim[1],0.5))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)    
    plt.contourf(lon_HYCOM,lat_HYCOM,sss_HYCOM,cmap=cmocean.cm.haline,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HMON-HYCOM SSS ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_HYCOM[oktime])[0:13],fontsize=16)
    
    file_name = folder_fig + 'HMON_HYCOM_SSS_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    
    #%% SST and velocity vectors
    kw = dict(levels = np.arange(temp_lim[0],temp_lim[1],0.5))

    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)    
    plt.contourf(lon_HYCOM,lat_HYCOM,sst_HYCOM,cmap=cmocean.cm.thermal,**kw)
    plt.plot(lon_forec_track,lat_forec_track,'.-',color='grey')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    plt.axis('scaled')
    plt.xlim(lon_lim[0],lon_lim[1])
    plt.ylim(lat_lim[0],lat_lim[1])
    plt.title('HMON-HYCOM SST ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n on '+str(time_HYCOM[oktime])[0:13],fontsize=16)
    
    q = plt.quiver(lon_HYCOM[::10], lat_HYCOM[::10],su_HYCOM[::10,::10],sv_HYCOM[::10,::10])
    #plt.quiverkey(q,np.min(lon_HYCOM)-5,np.max(lat_HYCOM)-5,1,"1 m/s",coordinates='data',color='k',fontproperties={'size': 14})
    
    file_name = folder_fig + 'HMON_HYCOM_SST_UV_' + cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
        
    #%% Temp 200 meters
    '''
    kw = dict(levels = np.linspace(10,25,31))
    
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)  
    plt.contourf(lon_HYCOM,lat_HYCOM,temp200_HYCOM,cmap=cmocean.cm.thermal) #,**kw)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    cbar.ax.tick_params(labelsize=14)
    plt.title('HMON-HYCOM Temperatute at 200 m \n on '+str(time_HYCOM[oktime])[0:13],fontsize=16)
    
    file_name = folder_fig + 'HWRF_POM_temp200_' + str(time_HYCOM[0])[0:10]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    '''
    
    #%% Salt 200 meters
    '''
    plt.figure()
    plt.contour(bath_lonsub,bath_latsub,bath_elevsub,levels=[0],colors='k')
    plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,levels=[0,10000],colors='papayawhip',alpha=0.5)  
    plt.contourf(lon_HYCOM,lat_HYCOM,salt200_HYCOM,cmap=cmocean.cm.haline) #,**kw)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14,labelpad=15)
    cbar.ax.tick_params(labelsize=14)
    plt.title('HMON-HYCOM Salinity at 200 m \n on '+str(time_HYCOM[oktime])[0:13],fontsize=16)
    
    file_name = folder_fig + 'HWRF_POM_salt200_' + str(time_HYCOM[0])[0:10]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
    '''
    
    #%% Figure temp transect along storm path
        
    lat_forec_tracku, ind = np.unique(lat_forec_track,return_index=True)
    lon_forec_tracku = lon_forec_track[ind]
    lon_forec_track_interp = np.interp(lat_hycom[:,0],lat_forec_tracku,lon_forec_tracku,left=np.nan,right=np.nan)
    lat_forec_track_interp = np.copy(lat_hycom[:,0])
    lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan
    
    lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
    lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
    
    oklon = np.round(np.interp(lon_forec_track_int,lon_hycom[0,:]-360,np.arange(len(lon_hycom[0,:])))).astype(int)
    oklat = np.round(np.interp(lat_forec_track_int,lat_hycom[:,0],np.arange(len(lat_hycom[:,0])))).astype(int)
    
    trans_temp_HYCOM = temp_HMON_HYCOM[oklat,oklon,:]
    
    kw = dict(levels = np.arange(tempt_lim[0],tempt_lim[1],1))
    
    plt.figure()
    plt.contourf(lat_hycom[oklat,0],-depth_HYCOM,trans_temp_HYCOM.T,cmap=cmocean.cm.thermal,**kw)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.contour(lat_hycom[oklat,0],-depth_HYCOM,trans_temp_HYCOM.T,[26],color='k')
    cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.xlabel('Latitude ($^o$)',fontsize=14)
    plt.title('HMON-HYCOM Temperature ' + 'Storm ' + storm_id + ' Cycle ' + cycle + '\n along Forecasted Storm Track on '+ str(time_HYCOM[oktime])[0:13],fontsize=12)
    plt.ylim([-300,0])
    
    file = folder_fig + 'HMON_HYCOM_temp_along_forecasted_track_' + cycle
    plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)    
    

