#%% User input
# forecasting cycle to be used

yyyymmddhh = '2023090706'
stormname = 'lee'
stormid = '13l'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
folder_hwrf = scratch_folder + 'HWRF_2023/' + yyyymmddhh + '/' + stormid + '/'

#yyyymmddhh = '2021070312'
#stormname = 'elsa'
#stormid = '05l'
#scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
#folder_hwrf = scratch_folder + 'HWRF2021_oper/' + stormid + '/' + yyyymmddhh

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

################################################################################
#%% Reading HWRF/POM grid
print('Retrieving coordinates from POM')
pom_grid = xr.open_dataset(files_hwrf_pom[0],decode_times=False)
lon_pom = np.asarray(pom_grid['east_e'][:])
lat_pom = np.asarray(pom_grid['north_e'][:])
#zlev_pom = np.asarray(pom_grid['zz'][:])
#hpom = np.asarray(pom_grid['h'][:])
#zmatrix = np.dot(hpom.reshape(-1,1),zlev_pom.reshape(1,-1)).T
#zmatrix_pom = zmatrix.reshape(zlev_pom.shape[0],hpom.shape[0],hpom.shape[1])

################################################################################
#%% Get storm track from trak atcf files

#file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*')))[0]
file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*hwrf*atcfunix')))[0]
#file_track_hwrf = sorted(glob.glob(os.path.join(folder_hwrf,'*trak*hwrf*3hourly')))[0]

lon_forec_track_hwrf, lat_forec_track_hwrf, lead_time_hwrf, int_track_hwrf, _ = \
get_storm_track_and_int(file_track_hwrf,stormid[0:-1])

#################################################################################
#%% Read POM
lon_forec_track = lon_forec_track_hwrf
lat_forec_track = lat_forec_track_hwrf
lon = lon_pom
lat = lat_pom

for n,file in enumerate(files_hwrf_pom):
    print(file)
    pom = xr.open_dataset(file)
    #t = np.asarray(pom['time'][:])
    #timestamp = mdates.date2num(t)[0]
    ff = str(n*6)
    #ff = str(int(file.split('/')[-1].split('.')[-2]))
    if len(ff)==1:
        fhour = '00' + ff
    if len(ff)==2:
        fhour = '0' + ff
    if len(ff)==3:
        fhour = ff 

    #fhour = ['0'+ff if len(ff)==2 else '00'+ff][0]

    xlim = [lon_forec_track[::2][n]-5,lon_forec_track[::2][n]+5]
    ylim = [lat_forec_track[::2][n]-5,lat_forec_track[::2][n]+5]

    #xlim = [-70,-45]
    #ylim = [12,28]

    oklon = np.where(np.logical_and(lon[0,:]>xlim[0],lon[0,:]<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat[:,0]>ylim[0],lat[:,0]<ylim[1]))[0]

    lonn = lon_pom[oklat,:][:,oklon]
    latt = lat_pom[oklat,:][:,oklon]

    MODEL = xr.open_dataset(file)
    sst = np.asarray(MODEL['t'][0,0,oklat,:][:,oklon])
    sst[sst == 0] = np.nan

    MODEL0 = xr.open_dataset(files_hwrf_pom[0])
    sst_f0 = np.asarray(MODEL0['t'][0,0,oklat,:][:,oklon])
    sst_f0[sst_f0 == 0] = np.nan

    u0 = np.asarray(MODEL['u'][0,0,oklat,:][:,oklon])
    u0[u0 == 0] = np.nan
    v0 = np.asarray(MODEL['v'][0,0,oklat,:][:,oklon])
    v0[v0 == 0] = np.nan
    V0 = np.sqrt(u0**2 + v0**2)

    print(np.nanmin(sst-sst_f0))
    ############################################################################
    # Figure SST        

    #kw = dict(levels=np.arange(19,31.1,0.5))
    kw = dict(levels=np.linspace(19, 32, 27))
    fig,ax = plt.subplots()
    ctr = ax.contourf(lonn,latt,sst,cmap='jet',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(lonn,latt,sst,np.linspace(19, 32, 27),colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.plot(lon_forec_track[::2], lat_forec_track[::2],'.-',markersize=10,color='k',alpha=0.5)
    ax.plot(lon_forec_track[::2][n], lat_forec_track[::2][n],'o-',color='k',markeredgecolor='green',markersize=8)
    ax.axis('scaled')
    ax.set_ylim([np.min(latt),np.max(latt)])
    ax.set_xlim([np.min(lonn),np.max(lonn)])
    #ax.set_xlim([-70,-45])
    #ax.set_ylim([12,28])
    ax.set_title('HWRF: POM Forecast for '+ stormname.upper() + ' Init: ' + yyyymmddhh + ' F' + fhour + '\n' + 'Sea Surface Temperature')
    fname = stormid + '.' + yyyymmddhh + '.pom.sst.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()

    ############################################################################
    # Figure SST-SST_f000

    kw = dict(levels=np.arange(-4,4.1,0.5))
    fig,ax = plt.subplots()
    ctr = ax.contourf(lonn,latt,sst-sst_f0,cmap='RdBu_r',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$^oC$',fontsize=14)
    cs = ax.contour(lonn,latt,sst-sst_f0,[0],colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ax.plot(lon_forec_track[::2], lat_forec_track[::2],'.-',markersize=10,color='k',alpha=0.5)
    ax.plot(lon_forec_track[::2][n], lat_forec_track[::2][n],'o-',color='k',markeredgecolor='green',markersize=8)
    ax.axis('scaled')
    ax.set_ylim([np.min(latt),np.max(latt)])
    ax.set_xlim([np.min(lonn),np.max(lonn)])
    ax.set_xlabel('dT min = ' + str(np.round(np.nanmin(sst-sst_f0),2)) +' $^oC$',fontsize=14)
    ax.set_title('HWRF: POM Forecast for ' + stormname.upper() + ' Init: ' + yyyymmddhh + ' F' + fhour + '\n' + 'Sea Surface Temperature Difference')
    fname = stormid + '.' + yyyymmddhh + '.pom.sst_diff.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    plt.close()

    ############################################################################
    # Figure surface velocity

    kw = dict(levels=np.arange(0,2.1,0.1))
    fig,ax = plt.subplots()
    ctr = ax.contourf(lonn,latt,V0,cmap='Spectral_r',**kw,extend='both')
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$m/s$',fontsize=14)
    q = ax.quiver(lonn[::2,::2],latt[::2,::2],u0[::2,::2],v0[::2,::2],scale=15)
    #ax.quiverkey(q,-72,20,1,'1 m/s')
    ax.plot(lon_forec_track[::2], lat_forec_track[::2],'.-',markersize=10,color='k',alpha=0.5)
    ax.plot(lon_forec_track[::2][n], lat_forec_track[::2][n],'o-',color='k',markeredgecolor='green',markersize=8)
    ax.set_ylim([np.min(latt),np.max(latt)])
    ax.set_xlim([np.min(lonn),np.max(lonn)])
    ax.set_title('HWRF: POM Forecast for '+ stormname.upper() + ' Init: ' + yyyymmddhh + ' F' + fhour + '\n' + 'Sea Surface Currents')
    fname = stormid + '.' + yyyymmddhh + '.pom.ssc.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    plt.close()

#################################################################################
