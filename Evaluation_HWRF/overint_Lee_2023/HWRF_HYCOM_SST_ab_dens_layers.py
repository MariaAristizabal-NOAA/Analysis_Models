#%% User input
cycle = '2023090706'
storm_id = '13l'
storm_name = 'lee'

# folder ab files HYCOM
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
rtofs_folder = scratch_folder + 'HWRF_HYCOM/' + cycle + '/' + storm_id + '/'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/fix/'
rtofs_grid_file = rtofs_grid_depth_folder + 'hafs_hycom_hat10.basin.regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'hafs_hycom_hat10.basin.regional.depth'

var_name = 'temp'
klayer = '1.00' #In this case is z in meters
#colormap='Spectral_r'
colormap='jet'
min_val = 18
max_val = 31
delta_val = 0.2  # delta in colorbar
delta_contour = 1 # delta in contour plot
units = '$^oC$'
lon_lim = [-54,-44]
lat_lim = [11,22]

'''
var_name = 'srfhgt'
klayer = '0'
colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = -110
max_val = 110
delta_val = 10
delta_contour = 30
units = 'Cm'
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]
'''

'''
var_name = 'salin'
klayer = '1'
colormap='GnBu_r'
#colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = 30
max_val = 40
delta_val = 0.5
delta_contour = 2
units = ' '
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]
'''

##################################################################################
def read_field_klayer(rtofs_file,var_name,klayer):

    lines = [line.rstrip() for line in open(rtofs_file+'.b')]
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)
    fld = ma.array([],fill_value=1.2676506002282294e+30)

    inFile = rtofs_file + '.a'

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == 'field':
                nheading = n + 1
                print(line.split()[0])

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[2] == var_name and line.split()[1] == klayer:
                nvar = n - nheading
                print(nvar)
                print(n)

    fid = open(inFile,'rb')
    #fid.seek((nvar-1)*4*(npad+ijdm),0)
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    fld = np.array(fld)
    fld = ma.reshape(fld,(jdm,idm))

    return fld

##################################################################################
#%% 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import glob
import sys
import numpy.ma as ma
import struct

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

##################################################################################
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

#%% Get storm track from trak atcf files

##################################################################################
file_track_hwrf = sorted(glob.glob(os.path.join(rtofs_folder,'*trak*hwrf*atcfunix')))[0]
file_track_hwrf = rtofs_folder + storm_name + storm_id + '.' + cycle + '.trak.hwrf.atcfunix'

lon_forec_track, lat_forec_track, lead_time, int_track, _ = \
get_storm_track_and_int(file_track_hwrf,storm_id[:-1])

##################################################################################
#%% finding list of hycom files
artofs_files = sorted(glob.glob(os.path.join(rtofs_folder,'*hwrf_rtofs_hat10_3z*.a')))

##################################################################################
#%% Reading HYCOM grid from ab files
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(rtofs_grid_file+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth_hycom = np.asarray(readdepth(rtofs_depth_file,'depth'))

##################################################################################
print('Reading ab files')
#%% file for refence cooling
n = 0
rtofs_file = artofs_files[n][:-2]
fld = read_field_klayer(rtofs_file,var_name,klayer)
fld = fld
if var_name == 'srfhgt':
    fld = fld*10

fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan
Field0 = fld2

for n,artofs_file in enumerate(artofs_files):
    print(artofs_file)
    fld = read_field_klayer(artofs_file[:-2],var_name,klayer)
    if var_name == 'srfhgt':
        fld = fld*10
    fld2 = np.copy(fld)
    mask = fld2 > 10**5
    fld2[mask] = np.nan
    xlim = [lon_forec_track[::2][n]-5,lon_forec_track[::2][n]+5]
    ylim = [lat_forec_track[::2][n]-5,lat_forec_track[::2][n]+5]
    oklon = np.logical_and(lon_hycom[0,:] >= xlim[0]+360,lon_hycom[0,:] <= xlim[1]+360)
    oklat = np.logical_and(lat_hycom[:,0] >= ylim[0],lat_hycom[:,0] <= ylim[1])
    #oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
    #oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
    maxval = np.nanmax(fld2[:,oklon][oklat,:])
    minval = np.nanmin(fld2[:,oklon][oklat,:])
    meanval = np.nanmean(fld2[:,oklon][oklat,:])
    print(maxval)
    print(minval)
    print(meanval)

    ff = artofs_file.split('/')[-1].split('.f')[-1].split('.')[-0]
    fhour = ff
    if len(ff)==1:
        fhour = '00' + ff
    if len(ff)==2:
        fhour = '0' + ff
    if len(ff)==3:
        fhour = ff

    lonn = lon_hycom[:,oklon][oklat,:]
    latt = lat_hycom[:,oklon][oklat,:]
    field = fld2[:,oklon][oklat,:]
    field0 = Field0[:,oklon][oklat,:]

    #kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
    kw = dict(levels=np.arange(19,32+0.1,0.5))
    #fig,ax = plt.subplots(figsize = (8,6.5))
    #fig,ax = plt.subplots(figsize=(10, 5))
    fig,ax = plt.subplots()
    ax.set_facecolor("bisque")
    #cs = ax.contour(lonn-360,latt,field,np.arange(min_val,max_val,delta_contour),colors='grey')
    cs = ax.contour(lonn-360,latt,field,np.linspace(19, 32, 27),colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ctr = ax.contourf(lonn-360,latt,field,cmap=colormap,**kw,extend='both')
    #ax.axis('scaled')
    cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
    cbar.set_label(units,fontsize = 14)
    #plt.title(var_name + ' layer ' + klayer + ' ' + cycle + ' f' + f ,fontsize=18)
    ax.plot(lon_forec_track[::2], lat_forec_track[::2],'.-',markersize=10,color='k',alpha=0.5)
    ax.plot(lon_forec_track[::2][n], lat_forec_track[::2][n],'o-',color='k',markeredgecolor='green',markersize=8)
    ax.set_title('HYCOM SST ' + cycle + ' f' + fhour ,fontsize=18)
    #plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
    #plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
    #ax.set_xlim(lon_lim)
    #ax.set_ylim(lat_lim)
    ax.set_ylim([np.min(latt),np.max(latt)])
    ax.set_xlim([np.min(lonn)-360,np.max(lonn)-360])
    fname = cycle + '.rtofs.sst.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()

    kw = dict(levels=np.arange(-4,4+0.1,0.5))
    #fig,ax = plt.subplots(figsize = (8,6.5))
    #fig,ax = plt.subplots(figsize=(10, 5))
    fig,ax = plt.subplots()
    ax.set_facecolor("bisque")
    #cs = ax.contour(lonn-360,latt,field,np.arange(min_val,max_val,delta_contour),colors='grey')
    cs = ax.contour(lonn-360,latt,field-field0,np.arange(-4,4+0.1,0.5),colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ctr = ax.contourf(lonn-360,latt,field-field0,cmap='RdBu_r',**kw,extend='both')
    #ax.axis('scaled')
    cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
    cbar.set_label(units,fontsize = 14)
    ax.plot(lon_forec_track[::2], lat_forec_track[::2],'.-',markersize=10,color='k',alpha=0.5)
    ax.plot(lon_forec_track[::2][n], lat_forec_track[::2][n],'o-',color='k',markeredgecolor='green',markersize=8)
    #plt.title(var_name + ' layer ' + klayer + ' ' + cycle + ' f' + f ,fontsize=18)
    ax.set_title('HYCOM SST-SST0 ' + cycle + ' f' + fhour ,fontsize=18)
    ax.set_xlabel('dT min = '+ str(np.round(np.nanmin(field-field0),2)),fontsize=14)
    #plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
    #plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
    #ax.set_xlim(lon_lim)
    #ax.set_ylim(lat_lim)
    ax.set_ylim([np.min(latt),np.max(latt)])
    ax.set_xlim([np.min(lonn)-360,np.max(lonn)-360])
    fname = cycle + '.rtofs.sst_diff.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    #plt.close()




