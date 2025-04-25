#%% User input
cycle = '20230908'

# folder ab files HYCOM
scratch1_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch2_folder = '/scratch2/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
artofs_file = scratch2_folder + 'RTOFS/rtofs.20230908/rtofs_glo.t00z.f120.archv.a'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch1_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

var_name = 'temp'
klayer = '1'
colormap='Spectral_r'
#colormap='jet'
#min_val = 18
#max_val = 31
min_val = -3.0
max_val = 33
delta_val = 2  # delta in colorbar
delta_contour = 1 # delta in contour plot
units = '$^oC$'
lon_lim = [-180,180]
lat_lim = [-90,90]

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
            if line.split()[0] == var_name and line.split()[4] == klayer:
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
fld = read_field_klayer(artofs_file[:-2],var_name,klayer)
fld = fld
if var_name == 'srfhgt':
    fld = fld*10

fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan
oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
field0 = fld2[:,oklon][oklat,:]

'''
print(artofs_file)
fld = read_field_klayer(artofs_file[:-2],var_name,klayer)
if var_name == 'srfhgt':
    fld = fld*10
fld2 = np.copy(fld)
mask = fld2 > 10**5
fld2[mask] = np.nan

oklon = np.logical_and(lon_hycom[0,:] >= lon_lim[0]+360,lon_hycom[0,:] <= lon_lim[1]+360)
oklat = np.logical_and(lat_hycom[:,0] >= lat_lim[0],lat_hycom[:,0] <= lat_lim[1])
'''

maxval = np.nanmax(fld2[:,oklon][oklat,:])
minval = np.nanmin(fld2[:,oklon][oklat,:])
meanval = np.nanmean(fld2[:,oklon][oklat,:])
print(maxval)
print(minval)
print(meanval)

ff = artofs_file.split('/')[-1].split('.f')[-1].split('.')[-3]
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

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
#kw = dict(levels=np.arange(19,32+0.1,0.5))
#fig,ax = plt.subplots(figsize = (8,6.5))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
#cs = ax.contour(lonn-360,latt,field,np.linspace(19, 32, 27),colors='k',alpha=0.3)
#cs = ax.contour(field,np.linspace(19, 32, 27),colors='k',alpha=0.3)
#ax.clabel(cs,cs.levels,inline=True,fmt='%1.1f',fontsize=10)
#ctr = ax.contourf(lonn-360,latt,field,cmap=colormap,**kw,extend='both')
#ctr = ax.contourf(field,cmap=colormap,**kw,extend='both')
ctr = ax.contourf(fld2,cmap=colormap,**kw,extend='both')
#ax.axis('scaled')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
#plt.title(var_name + ' layer ' + klayer + ' ' + cycle + ' f' + f ,fontsize=18)
ax.set_title('RTOFS SST ' + cycle + ' f' + fhour ,fontsize=18)
#plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
#plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
#ax.set_xlim(lon_lim)
#ax.set_ylim(lat_lim)
fname = cycle + '.rtofs.sst.f' + fhour + '.png'
#fig.savefig(fname,bbox_inches='tight')


