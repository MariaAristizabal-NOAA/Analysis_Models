#%% User input
# forecasting cycle to be used

# Milton
cycle = '2024100606'
storm_num = '14'
storm_id = '14l'
basin = 'al'
storm_name = 'Milton'

exp_names = ['RTOFS','RTOFS_v2.5.test01','HFSA_oper','HFSB_oper','hafsv2p0p1a_2024rt_cplinit_AOBS']
exp_labels = ['RTOFS_oper','RTOFSv2.5','HFSA_oper','HFSB_oper','HAFS_cplinit_AOBS']
exp_colors = ['magenta','salmon','purple','lime','cyan']
hafs = [' ',' ','hfsa','hfsb','hfsa']
ocean = ['rtofs','rtofs','mom6','hycom','mom6']
markers = ['p','p','o','o','o']

lon_lim = [-100,-60.0]
lat_lim = [5.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'
best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom
folder_myutils= '/home/Maria.Aristizabal/Utils/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'

#================================================================
# Calculate ocean heat content
def OHC(tempp,saltt,zl):
    cp = 3985 #Heat capacity in J/(kg K)
    zl_array = np.reshape(np.tile(zl,(tempp.shape[1]*tempp.shape[2],1)).T,(tempp.shape[0],tempp.shape[1],tempp.shape[2]))
    no26 = temp < 26
    tempp[no26] = np.nan
    saltt[no26] = np.nan
    density = dens(saltt,tempp,zl_array)
    rho0 = np.nanmean(density,axis=0)
    zl_array_fac = (zl_array[0:-1,:,:] + zl_array[1:,:,:])/2
    zero_array = np.zeros((1,tempp.shape[1],tempp.shape[2]))
    bott_array = np.ones((1,tempp.shape[1],tempp.shape[2]))* zl_array_fac[-1,0,0] + (zl[-1] - zl[-2])
    zl_array_face = np.vstack((zero_array,zl_array_fac,bott_array))
    dz_array = np.diff(zl_array_face,axis=0)
    ohc = np.abs(cp * rho0 * np.nansum((tempp-26)*dz_array,axis=0)) * 10**(-7) # in kJ/cm^2

    return ohc

################################################################################
def T100(depth,temp):
    # This function calculates the depth average temperature in the top 100
    # meters
    # Inputs
    # depth, temp: 1D vectors depth and temperature
    # Output
    # T100: depth average temperature in the top 100 meters

    okd = np.abs(depth) <= 100
    okv = np.isfinite(temp[okd])
    tempp = temp[okd][okv]
    dept = np.abs(depth[okd][okv])
    if len(tempp) != 0:
        T100 = np.trapz(tempp,dept)/np.max(dept)
    else:
        T100 = np.nan

    return T100

################################################################################
def Potential_energy_anomaly100(depth,dens):
    # This function calculates the potential energy anomaly
    # (Simpson J, Brown J, Matthews J, Allen G (1990) Tidal straining, density
    # currents and stirring in the control of estuarine stratification.
    # Estuaries 13(2):125▒~@~S132), in the top 100 meters
    # Inputs
    # depth, dens: 1D vectors depth and density
    # Output
    # PEA: potential energy anomaly in J/m^3

    g = 9.8 #m/s
    dindex = np.fliplr(np.where(np.asarray(np.abs(depth)) <= 100))[0]
    if len(dindex) == 0:
        PEA = np.nan
    else:
        zz = np.asarray(np.abs(depth[dindex]))
        denss = np.asarray(dens[dindex])
        ok = np.isfinite(denss)
        z = zz[ok]
        densi = denss[ok]
        if len(z)==0 or len(densi)==0 or np.min(zz) > 10 or np.max(zz) < 30:
            PEA = np.nan
        else:
            if z[-1] - z[0] > 0:
                # So PEA is < 0
                # sign = -1
                # Adding 0 to sigma integral is normalized
                z = np.append(0,z)
            else:
                # So PEA is < 0
                # sign = 1
                # Adding 0 to sigma integral is normalized
                z = np.flipud(z)
                z = np.append(0,z)
                densit = np.flipud(densi)

            # adding density at depth = 0
            densitt = np.interp(z,z[1:],densit)
            density = np.flipud(densitt)

            # defining sigma
            max_depth = np.nanmax(zz[ok])
            sigma = -1*z/max_depth
            sigma = np.flipud(sigma)

            rhomean = np.trapz(density,sigma,axis=0)
            drho = rhomean - density
            torque = drho * sigma
            PEA = g * max_depth * np.trapz(torque,sigma,axis=0)

    return PEA

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

################################################################################
import sys
import os
import glob
import xarray as xr
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, GOFS_coor_to_geo_coord 

#sys.path.append(folder_uom)
#from Upper_ocean_metrics import T100, Potential_energy_anomaly100

from eos80 import dens

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
'''
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')
'''

#################################################################################
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

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

ohc_forec_track_hafs = np.empty((len(folder_exps),300))
ohc_forec_track_hafs[:] = np.nan
ohc_forec_track_hafs_mean = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_mean[:] = np.nan
ohc_forec_track_hafs_max = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_max[:] = np.nan
ohc_forec_track_hafs_min = np.empty((len(folder_exps),300))
ohc_forec_track_hafs_min[:] = np.nan
lon_forec_track_hafs = np.empty((len(folder_exps),300)) 
lon_forec_track_hafs[:] = np.nan
lat_forec_track_hafs = np.empty((len(folder_exps),300)) 
lat_forec_track_hafs[:] = np.nan
#%% Loop the experiments
for i,folder in enumerate(folder_exps[0:1]):
    print(folder)

    #%% Get list files
    if ocean == 'mom6':
        file_ocean = folder + '/' + storm_id + '.' + cycle + '.' + hafs[i] + '.mom6.f000.nc'
        nc = xr.open_dataset(file_ocean)
        temp = np.asarray(nc['temp'][0,:,:])
        salt = np.asarray(nc['so'][0,:,:])
        lon = np.asarray(nc['xh'])
        lat = np.asarray(nc['yh'])
        zl = np.asarray(nc['z_l'])
        time = np.asarray(nc['time'])
        t100_ocean= T100(temp,zl)
        ohc_ocean = OHC(temp,salt,zl)

    if ocean == 'hycom':
        file_ocean = folder + '/' + storm_id + '.' + cycle + '.' + hafs[i] + '.hycom.3z.f000.nc'
        nc = xr.open_dataset(file_ocean)
        temp = np.asarray(nc['temperature'][0,:,:,:])
        salt = np.asarray(nc['salinity'][0,:,:,:])
        lon = np.asarray(nc['Longitude'])
        lat = np.asarray(nc['Latitude'])
        depth = np.asarray(nc['Z'][:])
        time = np.asarray(nc['MT'][:])
        t100_ocean= T100(temp,zl)
        ohc_ocean = OHC(temp,salt,depth)

    if ocean[i] == 'rtofs':
        file_ocean = folder + '/rtofs_glo_3dz_n024_6hrly_hvr_US_east.nc'
        nc = xr.open_dataset(file_ocean)
        temp = np.asarray(nc['temperature'][0,:,:,:])
        salt = np.asarray(nc['salinity'][0,:,:,:])
        lon = np.asarray(nc['Longitude'][:])
        lat = np.asarray(nc['Latitude'][:])
        zl = np.asarray(nc['Depth'][:])
        zl_array = np.tile(zl,(lat.shape[1],lat.shape[0],1)).T
        #t100_ocean, t100b = T100(temp,zl_array)

        density = dens(salt,temp,zl_array)
        t100_ocean = np.empty((lat.shape[0],lat.shape[1]))
        t100_ocean[:] = np.nan
        pea100 = np.empty((lat.shape[0],lat.shape[1]))
        pea100[:] = np.nan
        for y in np.arange(lat.shape[0]):
            print(y)
            for x in np.arange(lat.shape[1]):
                t100_ocean[y,x] = T100(zl,temp[:,y,x])
                pea100[y,x] = Potential_energy_anomaly100(zl,density[:,y,x])

        ohc_ocean = OHC(temp,salt,zl)

    #================================================================
    # Constrain lon limits between -180 and 180 so it does not conflict with the cartopy projection PlateCarree
    '''
    lon[lon>180] = lon[lon>180] - 360
    lon[lon<-180] = lon[lon<-180] + 360
    sort_lon = np.argsort(lon)
    lon = lon[sort_lon]

    # define grid boundaries
    lonmin_new = np.min(lon)
    lonmax_new = np.max(lon)
    latmin = np.min(lat)
    latmax = np.max(lat)
    print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))

    # Sort OHC according with sorted lon 
    ohc_ocean = ohc_ocean[:,sort_lon]
    '''

    '''
    # Spatial interpolation of ohc_hafs to the NESDIS spatial resolution
    lono,lato = np.meshgrid(lon_nesdis,lat_nesdis)
    lonh,lath = np.meshgrid(lon,lat)
    
    interpolator = LinearNDInterpolator(list(zip(np.ravel(lonh),np.ravel(lath))),np.ravel(ohc_hafs))
    ohc_from_hafs_to_nesdis = interpolator((lono,lato))
    '''

    ##############################################
    # vertical profiles in areas of interest
    # Warm core eddy: High OHC, High T100 and mid-range PEA (easily mixed) 
    target_lon = -87
    target_lat = 26
    oklat1,oklon1 = find_grid_position_hycom(lat,lon,target_lat,target_lon)

    # Caribbean warm pool: high OHC, high T100 and mid-range PEA (easily mixed) 
    target_lon = -83
    target_lat = 20
    oklat2,oklon2 = find_grid_position_hycom(lat,lon,target_lat,target_lon)

    # West Florida shell: low OHC, high T100 and low PEA (easily mixed) 
    target_lon = -83
    target_lat = 26
    oklat3,oklon3 = find_grid_position_hycom(lat,lon,target_lat,target_lon)

    # Caribbean barrier layer: low OHC, low T100 and high PEA (difficult mixed) 
    target_lon = -65
    target_lat = 14
    oklat4,oklon4 = find_grid_position_hycom(lat,lon,target_lat,target_lon)

    # North Atlantic: Low OHC, low T100 and low PEA (easily mixed) 
    target_lon = -70
    target_lat = 25
    oklat5,oklon5 = find_grid_position_hycom(lat,lon,target_lat,target_lon)

    temp = np.asarray(nc['temperature'][0,:,:,:])
    plt.figure(figsize=(4,8))
    plt.plot(temp[:,oklat1,oklon1],-zl,'.-',label='Warm Core Eddy')
    plt.plot(temp[:,oklat2,oklon2],-zl,'.-',label='Caribbean')
    plt.plot(temp[:,oklat3,oklon3],-zl,'.-',label='Florida Shell')
    plt.plot(temp[:,oklat4,oklon4],-zl,'.-',label='Barrier Layer')
    plt.plot(temp[:,oklat5,oklon5],-zl,'.-',label='North Atlantic')
    plt.ylim([-100,0])
    plt.xlim([20,32])
    plt.xlabel('Temperature ($^0C$)',fontsize=14)
    plt.legend()

    salt = np.asarray(nc['salinity'][0,:,:,:])
    density = dens(salt,temp,zl_array)
    plt.figure(figsize=(4,8))
    plt.plot(density[:,oklat1,oklon1],-zl,'.-',label='Warm Core Eddy')
    plt.plot(density[:,oklat2,oklon2],-zl,'.-',label='Caribbean')
    plt.plot(density[:,oklat3,oklon3],-zl,'.-',label='Florida Shell')
    plt.plot(density[:,oklat4,oklon4],-zl,'.-',label='Barrier Layer')
    plt.plot(density[:,oklat5,oklon5],-zl,'.-',label='North Atlantic')
    plt.ylim([-100,0])
    plt.xlim([1020,1026])
    plt.xlabel('Density ($kg/m^3$)',fontsize=14)
    plt.legend()

    #%% Figure OHC all domain HAFS
    lev = np.arange(0,161,20)
    kw = dict(levels=np.arange(0,161,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    #plt.contour(lon,lat,ohc_ocean,lev,colors='grey',alpha=0.5)
    plt.contour(lon,lat,ohc_ocean,[100],colors='k')
    plt.contourf(lon,lat,ohc_ocean,cmap='Spectral_r',extend='max',**kw)
    cbar = plt.colorbar(extendrect=True)
    #plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    cbar.ax.set_ylabel('$kJ/cm^2$',fontsize = 14)
    plt.plot(lon[oklat1,oklon1],lat[oklat1,oklon1],'*',color='limegreen')
    plt.plot(lon[oklat2,oklon2],lat[oklat2,oklon2],'*',color='limegreen')
    plt.plot(lon[oklat3,oklon3],lat[oklat3,oklon3],'*',color='limegreen')
    plt.plot(lon[oklat4,oklon4],lat[oklat4,oklon4],'*',color='limegreen')
    plt.plot(lon[oklat5,oklon5],lat[oklat5,oklon5],'*',color='limegreen')
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('OHC HAFS '+ cycle)
    plt.legend()
    #plt.text(-95,5,file.split('/')[-1].split('.')[-2],fontsize=14)
    #plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

    #%% Figure T100 all domain HAFS
    lev = np.arange(25.5,31,0.5)
    kw = dict(levels=np.arange(25.5,31,0.5))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon,lat,t100_ocean,[29],colors='k')
    plt.contour(lon,lat,t100_ocean,[26],colors='grey')
    plt.contourf(lon,lat,t100_ocean,cmap='Spectral_r',extend='both',**kw)
    cbar = plt.colorbar(extendrect=True)
    #plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    plt.plot(lon[oklat1,oklon1],lat[oklat1,oklon1],'*',color='limegreen')
    plt.plot(lon[oklat2,oklon2],lat[oklat2,oklon2],'*',color='limegreen')
    plt.plot(lon[oklat3,oklon3],lat[oklat3,oklon3],'*',color='limegreen')
    plt.plot(lon[oklat4,oklon4],lat[oklat4,oklon4],'*',color='limegreen')
    plt.plot(lon[oklat5,oklon5],lat[oklat5,oklon5],'*',color='limegreen')
    cbar.ax.set_ylabel('$^oC$',fontsize = 14)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('T100 '+ cycle)
    plt.legend()
    
    #%% Figure pea100 all domain HAFS
    kw = dict(levels=np.arange(0,301,20))
    plt.figure(figsize=(8,5))
    plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
    plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
    plt.contour(lon,lat,pea100,[200],colors='k',alpha=0.5)
    plt.contourf(lon,lat,pea100,cmap='Spectral_r',extend='both',**kw)
    cbar = plt.colorbar(extendrect=True)
    plt.plot(lon[oklat1,oklon1],lat[oklat1,oklon1],'*',color='limegreen')
    plt.plot(lon[oklat2,oklon2],lat[oklat2,oklon2],'*',color='limegreen')
    plt.plot(lon[oklat3,oklon3],lat[oklat3,oklon3],'*',color='limegreen')
    plt.plot(lon[oklat4,oklon4],lat[oklat4,oklon4],'*',color='limegreen')
    plt.plot(lon[oklat5,oklon5],lat[oklat5,oklon5],'*',color='limegreen')
    #plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=1)
    cbar.ax.set_ylabel('$J/m^3$',fontsize = 14)
    plt.axis('scaled')
    plt.ylim(lat_lim)
    plt.xlim(lon_lim)
    plt.title('PEA '+ cycle)
    plt.legend()

#################################################################################
