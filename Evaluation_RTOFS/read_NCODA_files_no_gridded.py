#%% User input
cycle = '20240911'

# folder ab files HYCOM
scratch_folder = '/work/noaa/hwrf/noscrub/maristiz/'
artofs_file = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/rtofs_glo.t00z.f06.archv.a'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

ncoda_2dvar_file_jmindat = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_jmindat'
ncoda_2dvar_file_voldata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/seatmp_sfc_1o4500x3298_2024091106_0000_voldata'
ncoda_2dvar_file_ams_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/amssst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_goes_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/goesst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_him_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/himsst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_shp_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/shpsst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_jps_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/jpssst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_msg_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/msgsst_sfc_1o4500x3298_2024091106_0000_rawdata'
ncoda_2dvar_file_mtp_rawdata = scratch_folder + 'RTOFS_v2.4/rtofs.20240911/mtpsst_sfc_1o4500x3298_2024091106_0000_rawdata'


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
# Read ncoda 2DVAR file
#ncoda_2dvar_file = ncoda_2dvar_file_jmindat
#ncoda_2dvar_file = ncoda_2dvar_file_voldata
#ncoda_2dvar_file = ncoda_2dvar_file_goes_rawdata
#ncoda_2dvar_file = ncoda_2dvar_file_ams_rawdata
#ncoda_2dvar_file = ncoda_2dvar_file_him_rawdata
#ncoda_2dvar_file = ncoda_2dvar_file_shp_rawdata
#ncoda_2dvar_file = ncoda_2dvar_file_jps_rawdata
#ncoda_2dvar_file = ncoda_2dvar_file_msg_rawdata
ncoda_2dvar_file = ncoda_2dvar_file_mtp_rawdata

binary_file_ncoda = open(ncoda_2dvar_file,'rb')
binary_file_ncoda.seek(0)

#n = 10000000
n = 10000
fld_ncoda = binary_file_ncoda.read(n*4)
fld_ncoda_unpack = struct.unpack('>'+str(n)+'f',fld_ncoda)
fld_ncoda_array = np.array(fld_ncoda_unpack)

maxval_ncoda = np.nanmax(fld_ncoda_array)
minval_ncoda = np.nanmin(fld_ncoda_array)
print('from ncoda file ',maxval_ncoda)
print('from ncoda file ',minval_ncoda)

#############################################
#colormap='Spectral_r'
#min_val = 0
#max_val = 2085752.0
#delta_val = 208575
#units = '$^oC$'
#units = ' '

fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(fld_ncoda_array)
plt.ylim([-90,90])
plt.title(ncoda_2dvar_file.split('/')[-1] ,fontsize=18)

'''
kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_ncoda,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title(ncoda_2dvar_file.split('/')[-1] ,fontsize=18)
fname = cycle + '_ncoda_sst.png'
fig.savefig(fname,bbox_inches='tight')
'''
