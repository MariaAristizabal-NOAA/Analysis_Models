#%% User input
cycle = '20240911'

# folder ab files HYCOM
scratch_folder = '/work/noaa/hwrf/noscrub/maristiz/'
artofs_file = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/rtofs_glo.t00z.f06.archv.a'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

ncoda_2dvar_file_temp_analinc = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_analinc'
ncoda_2dvar_file_salt_analinc = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/salint_sfc_1o4500x3298_2024091106_0000_analinc'
ncoda_2dvar_file_climerr = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_climerr'
ncoda_2dvar_file_climfld = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_climfld'
ncoda_2dvar_file_fcsterr = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_fcsterr'
ncoda_2dvar_file_flowfld = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_flowfld'
#ncoda_2dvar_file_jmindat = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_jmindat'
ncoda_2dvar_file_obsdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_obsdata'
ncoda_2dvar_file_statdat = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_statdat'
#ncoda_2dvar_file_voldata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_voldata'

##################################################################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import glob
import sys
import numpy.ma as ma
import struct

folder_utils4hycom= '/home/maristiz/Repos/NCEP_scripts/'
folder_myutils= '/home/maristiz/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

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
ijdm = idm*jdm
lon_hycom = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
lat_hycom = np.array(readgrids(rtofs_grid_file,'plat:',[0]))

depth_hycom = np.asarray(readdepth(rtofs_depth_file,'depth'))

#############################
# Read ncoda 2DVAR file
ncoda_2dvar_file = ncoda_2dvar_file_temp_analinc
#ncoda_2dvar_file = ncoda_2dvar_file_salt_analinc
#ncoda_2dvar_file = ncoda_2dvar_file_climerr
#ncoda_2dvar_file = ncoda_2dvar_file_climfld
#ncoda_2dvar_file = ncoda_2dvar_file_fcsterr
#ncoda_2dvar_file = ncoda_2dvar_file_flowfld
#ncoda_2dvar_file = ncoda_2dvar_file_jmindat
#ncoda_2dvar_file = ncoda_2dvar_file_obsdata
#ncoda_2dvar_file = ncoda_2dvar_file_statdat
#ncoda_2dvar_file = ncoda_2dvar_file_voldata

binary_file_ncoda = open(ncoda_2dvar_file,'rb')
binary_file_ncoda.seek(0)
fld_ncoda = binary_file_ncoda.read(ijdm*4)
fld_ncoda_unpack = struct.unpack('>'+str(ijdm)+'f',fld_ncoda)
fld_ncoda_array = np.array(fld_ncoda_unpack)
fld_ncoda_reshape = ma.reshape(fld_ncoda_array,(jdm,idm))
    
fld_ncoda_array_new_fill_value = np.copy(fld_ncoda_array)
mask = fld_ncoda_array_new_fill_value < -100
fld_ncoda_array_new_fill_value[mask] = 1.2676506002282294e+30
fld_ncoda_array_new_fill_value_reshape = ma.reshape(fld_ncoda_array_new_fill_value,(jdm,idm))
fld_ncoda_array_new_fill_value_byt = [struct.pack('>f',value) for value in fld_ncoda_array_new_fill_value]
fld_ncoda_array_new_fill_value_bytes = b''.join(fld_ncoda_array_new_fill_value_byt)
    
fld2_ncoda = np.copy(fld_ncoda_reshape)
mask = fld2_ncoda < -100
fld2_ncoda[mask] = np.nan
    
maxval_ncoda = np.nanmax(fld2_ncoda)
minval_ncoda = np.nanmin(fld2_ncoda)
print('from ncoda file ',maxval_ncoda)
print('from ncoda file ',minval_ncoda)

#############################################
colormap='Spectral_r'
min_val = -5
max_val = 5
delta_val = 1
units = '$^oC$'

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_ncoda,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title(ncoda_2dvar_file.split('/')[-1] ,fontsize=18)
fname = cycle + '_ncoda_sst.png'
fig.savefig(fname,bbox_inches='tight')

