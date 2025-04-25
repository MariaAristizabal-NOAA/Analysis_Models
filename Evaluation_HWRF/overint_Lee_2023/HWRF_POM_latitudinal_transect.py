#%% User input
# forecasting cycle to be used

# Lee 13L
yyyymmddhh = '2023090706'
stormname = 'Lee'
stormid = '13l'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
folder_hwrf = scratch_folder + 'HWRF_2023/' + yyyymmddhh + '/' + stormid + '/'

# Transect lon and lat limits
lon_lim = [-59.8,-59.8]
lat_lim = [15.8,22.4]

################################################################################
#%% Get storm track from trak atcf files
def get_storm_track_and_int(file_track,storm_num):

    ff = open(file_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    lead_time = []
    intt = []
    rmww = []
    for l in f:
        if l.split(',')[1].strip() == storm_num:
            lat = float(l.split(',')[6][0:4])/10
            if l.split(',')[6][4] == 'N':
                lat = lat
            else:
                lat = -lat
            lon = float(l.split(',')[7][0:5])/10
            if l.split(',')[7][4] == 'E':
                lon = lon
            else:
                lon = -lon
            latt.append(lat)
            lont.append(lon)
            lead_time.append(int(l.split(',')[5][1:4]))
            intt.append(float(l.split(',')[8]))
            rmww.append(float(l.split(',')[19]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    intt = np.asarray(intt)
    rmww = np.asarray(rmww)
    lead_time, ind = np.unique(lead_time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    rmw_track = rmww[ind]

    return lon_track, lat_track, lead_time, int_track, rmw_track

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#################################################################################
#%% Get list files
files_hwrf_pom = sorted(glob.glob(os.path.join(folder_hwrf,'*pom*00*.nc')))
file_hwrf_pom_grid = sorted(glob.glob(os.path.join(folder_hwrf,'*pom*grid*.nc')))[0]

################################################################################
#%% Reading HWRF/POM grid
print('Retrieving coordinates from POM')
pom_grid = xr.open_dataset(file_hwrf_pom_grid,decode_times=False)
lon_pom = np.asarray(pom_grid['east_e'][:])
lat_pom = np.asarray(pom_grid['north_e'][:])
zlev_pom = np.asarray(pom_grid['zz'][:])
hpom = np.asarray(pom_grid['h'][:])
zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
zmatrix_pom = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])

################################################################################
#%% Get storm track from trak atcf files

#file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*')))[0]
file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*hwrf*atcfunix')))[0]
#file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*hwrf*3hourly')))[0]

lon_forec_track_hwrf, lat_forec_track_hwrf, lead_time_hwrf, int_track_hwrf, _ = \
get_storm_track_and_int(file_track_hwrf,stormid[0:-1])

#################################################################################
#%% temp along storm path
#lat_forec_track = lat_forec_track_hwrf
#lon_forec_track = lon_forec_track_hwrf
lat_forec_track = np.arange(lat_lim[0],lat_lim[1],0.5)
lon_forec_track = np.tile(lon_lim[0],len(lat_forec_track))

lon_forec_track_interp = np.interp(lat_pom[:,0],lat_forec_track,lon_forec_track,left=np.nan,right=np.nan)
lat_forec_track_interp = np.copy(lat_pom[:,0])
lat_forec_track_interp[np.isnan(lon_forec_track_interp)] = np.nan

lon_forec_track_int = lon_forec_track_interp[np.isfinite(lon_forec_track_interp)]
lat_forec_track_int = lat_forec_track_interp[np.isfinite(lat_forec_track_interp)]
lon_forec_track_int[0:len(lon_forec_track_int)] = lon_forec_track_int
lat_forec_track_int[0:len(lat_forec_track_int)] = lat_forec_track_int

lon_forec_track_int_g = lon_forec_track_int
lat_forec_track_int_g = lat_forec_track_int

oklon = np.round(np.interp(lon_forec_track_int,lon_pom[0,:],np.arange(len(lon_pom[0,:])))).astype(int)
oklat = np.round(np.interp(lat_forec_track_int,lat_pom[:,0],np.arange(len(lat_pom[:,0])))).astype(int)

MODEL0 = xr.open_dataset(files_hwrf_pom[0])
temp0 = temp = np.asarray(MODEL0['t'][0,:,oklat,:][:,:,oklon])[:,:,0] 
temp0[temp0 == 0] = np.nan

depth_pom = zmatrix_pom[:,oklat,:][:,:,oklon][:,:,0]
latt_pom =  np.tile(lat_pom[oklat,0],(depth_pom.shape[0],1))

for n,file in enumerate(files_hwrf_pom[0:2]):
    print(file)
    pom = xr.open_dataset(file)
    ff = str(n*6)
    if len(ff)==1:
        fhour = '00' + ff
    if len(ff)==2:
        fhour = '0' + ff
    if len(ff)==3:
        fhour = ff 

    MODEL = xr.open_dataset(file)
    temp = np.asarray(MODEL['t'][0,:,oklat,:][:,:,oklon])[:,:,0]
    temp[temp == 0] = np.nan

    print(np.nanmin(temp-temp0))
    #############################################################
    # Temp        

    kw = dict(levels=np.arange(15,31.1,0.5))
    fig,ax = plt.subplots(figsize=(8,4))
    ctr = ax.contourf(latt_pom,depth_pom,temp,cmap='Spectral_r',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(latt_pom,depth_pom,temp,[26],colors='k')
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.set_ylim([-300,0])
    ax.set_title('HWRF: POM Forecast for '+ stormname.upper() + ' Init: ' + yyyymmddhh + ' F' + fhour + '\n' + 'Temperature')
    fname = stormid + '.' + yyyymmddhh + '.pom.temp_trans.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()

    #####################################################################
    # Figure temp-temp_f000

    kw = dict(levels=np.arange(-4,4.1,0.2))
    fig,ax = plt.subplots(figsize=(8,4))
    ctr = ax.contourf(latt_pom,depth_pom,temp-temp0,cmap='seismic',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(latt_pom,depth_pom,temp-temp0,[0],colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.set_ylim([-300,0])
    ax.set_title('HWRF: POM Forecast for '+ stormname.upper() + ' Init: ' + yyyymmddhh + ' F' + fhour + '\n' + 'Temperature difference')
    ax.set_xlabel('dt min = ' + str(np.round(np.nanmin(temp-temp0),2)),fontsize=14) 
    fname = stormid + '.' + yyyymmddhh + '.pom.temp_diff_trans.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()

    ############################################################################
