#%% User input
cycle = '2023090706'
storm_id = '13l'

# folder ab files HYCOM
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
rtofs_folder = scratch_folder + 'HMON_2023/' + cycle + '/' + storm_id + '/'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/fix/'
rtofs_grid_file = rtofs_grid_depth_folder + 'hafs_hycom_hat10.basin.regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'hafs_hycom_hat10.basin.regional.depth'

var_name = 'temp'
zdepth = '1.00'
colormap='Spectral_r'
min_val = 18
max_val = 31
delta_val = 0.2  # delta in colorbar
delta_contour = 1 # delta in contour plot
units = '$^oC$'
#lon_lim = [-48,-38]
#lat_lim = [34,46]
lon_lim = [-47,-41]
lat_lim = [34,42]

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
def read_field_klayer(rtofs_file,var_name,zdepth):

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
            if line.split()[2] == var_name and line.split()[1] == zdepth:
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
folder_myutils= '/home/Maria.Aristizabal/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

sys.path.append(folder_myutils)
from my_models_utils import geo_coord_to_HYCOM_coord

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

##################################################################################
#%% finding list of hycom files
artofs_files = sorted(glob.glob(os.path.join(rtofs_folder,'*hmon_rtofs*3z*.a')))

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
rtofs_file = artofs_files[0][:-2]
fld = read_field_klayer(rtofs_file,var_name,zdepth)
fld = fld
if var_name == 'srfhgt':
    fld = fld*10

fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan
oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
field0 = fld2[:,oklon][oklat,:]

for artofs_file in artofs_files:
    fld = read_field_klayer(artofs_file[:-2],var_name,zdepth)
    if var_name == 'srfhgt':
        fld = fld*10
    fld2 = np.copy(fld)
    mask = fld2 > 10**5
    fld2[mask] = np.nan
    oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
    oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
    maxval = np.nanmax(fld2[:,oklon][oklat,:])
    minval = np.nanmin(fld2[:,oklon][oklat,:])
    meanval = np.nanmean(fld2[:,oklon][oklat,:])
    print(maxval)
    print(minval)
    print(meanval)

    ff = artofs_file.split('/')[-1].split('.f')[-1].split('.')[0]

    if len(ff)==1:
        fhour = '00' + ff
    if len(ff)==2:
        fhour = '0' + ff
    if len(ff)==3:
        fhour = ff

    lonn = lon_hycom[:,oklon][oklat,:]
    latt = lat_hycom[:,oklon][oklat,:]
    field = fld2[:,oklon][oklat,:]

    kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
    #fig,ax = plt.subplots(figsize = (8,6.5))
    fig,ax = plt.subplots()
    ax.set_facecolor("bisque")
    #cs = ax.contour(lonn-360,latt,field,np.arange(min_val,max_val,delta_contour),colors='grey')
    cs = ax.contour(lonn-360,latt,field,[26,27,28,29,30],colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ctr = ax.contourf(lonn-360,latt,field,cmap=colormap,**kw,extend='both')
    #ax.axis('scaled')
    cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
    cbar.set_label(units,fontsize = 14)
    #plt.title(var_name + ' layer ' + klayer + ' ' + cycle + ' f' + f ,fontsize=18)
    ax.set_title('HMON SST ' + cycle + ' f' + fhour ,fontsize=18)
    #plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
    #plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
    ax.set_xlim(lon_lim)
    ax.set_ylim(lat_lim)
    fname = cycle + '.rtofs.sst.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    plt.close()

    kw = dict(levels=np.arange(-4,4+0.1,0.2))
    #fig,ax = plt.subplots(figsize = (8,6.5))
    fig,ax = plt.subplots()
    ax.set_facecolor("bisque")
    #cs = ax.contour(lonn-360,latt,field,np.arange(min_val,max_val,delta_contour),colors='grey')
    cs = ax.contour(lonn-360,latt,field-field0,[0],colors='k',alpha=0.3)
    ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
    ctr = ax.contourf(lonn-360,latt,field-field0,cmap='seismic',**kw,extend='both')
    #ax.axis('scaled')
    cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
    cbar.set_label(units,fontsize = 14)
    #plt.title(var_name + ' layer ' + klayer + ' ' + cycle + ' f' + f ,fontsize=18)
    ax.set_title('HMON SST-SST0 ' + cycle + ' f' + fhour ,fontsize=18)
    ax.set_xlabel('dT min = '+ str(np.round(np.nanmin(field-field0),2)),fontsize=14)
    #plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
    #plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
    ax.set_xlim(lon_lim)
    ax.set_ylim(lat_lim)
    fname = cycle + '.rtofs.sst_diff.f' + fhour + '.png'
    fig.savefig(fname,bbox_inches='tight')
    plt.close()




