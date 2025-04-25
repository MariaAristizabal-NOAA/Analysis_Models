#%% User input
cycle = '20240115'

# folder ab files HYCOM
scratch_folder = '/work/noaa/hwrf/noscrub/maristiz/'
artofs_file = scratch_folder + 'RTOFS/rtofs.20240115/rtofs_glo.t00z.f06.archv.test2.a'

# RTOFS grid file name
rtofs_grid_depth_folder = scratch_folder + 'RTOFS_fix/GRID_DEPTH_global/'
rtofs_grid_file = rtofs_grid_depth_folder + 'regional.grid'
rtofs_depth_file = rtofs_grid_depth_folder + 'regional.depth'

ncoda_2dvar_file_temp = scratch_folder + 'RTOFS/rtofs.20240115/seatmp_sfc_1o4500x3298_2024011506_0000_analfld'

ncoda_2dvar_file_salt = scratch_folder + 'RTOFS/rtofs.20240115/salint_sfc_1o4500x3298_2024011506_0000_analfld'

var_temp = 'temp'
var_salt = 'salin'
klayer = '1'
colormap='Spectral_r'

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

#sys.path.append(folder_myutils)
#from my_models_utils import geo_coord_to_HYCOM_coord

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
def replace_field_in_archv_from_2DVAR_NCODA(artofs_file,var_name,klayer,ncoda_2dvar_file):

    print('Reading ab files')
    rtofs_file = artofs_file[:-2] 
    
    lines = [line.rstrip() for line in open(rtofs_file+'.b')]
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)
    fld = ma.array([],fill_value=1.2676506002282294e+30)
    
    inFile_a = rtofs_file + '.a'
    
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
    
    fid = open(inFile_a,'rb')
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    fld = np.array(fld)
    fld_archv = ma.reshape(fld,(jdm,idm))
    
    maxval_raw_archv = np.nanmax(fld_archv)
    minval_raw_archv = np.nanmin(fld_archv)
    print('from raw archv file ',maxval_raw_archv)
    print('from raw archv file ',minval_raw_archv)
    
    fld2_archv = np.copy(fld_archv)
    mask = fld2_archv > 10**5
    fld2_archv[mask] = np.nan
    
    maxval_archv = np.nanmax(fld2_archv)
    minval_archv = np.nanmin(fld2_archv)
    print('from archv file ',maxval_archv)
    print('from archv file ',minval_archv)
    
    #############################
    # Read ncoda 2DVAR file
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
    
    #############################
    # Replace value in bynary a file
    binary_file = open(inFile_a,'r+b')
    binary_file.seek((nvar)*4*(npad+ijdm))
    binary_file.write(fld_ncoda_array_new_fill_value_bytes)
    binary_file.close()
    
    #############################
    # Replace max and min value for temp at surface layer in b file
    inFile_b = rtofs_file + '.b'
    bfile = open(inFile_b,'r+')
    
    linesb = open(inFile_b,'r+').readlines()
    
    for n,line in enumerate(linesb):
        if len(line) > 0:
            if line.split()[0] == var_name and line.split()[4] == klayer:
                print(n)
                linesb[n] =  line[0:44] + '{:.7E}'.format(minval_ncoda) + '   ' + '{:.7E}'.format(maxval_ncoda) + '\n'
    
    bfile.writelines(linesb)
    bfile.close()
    
    return fld2_archv, fld2_ncoda, nvar

#############################
# Replace SST on archv file
var_name = var_temp
fld2_temp_archv, fld2_temp_ncoda, nvar_temp = replace_field_in_archv_from_2DVAR_NCODA(artofs_file,var_name,klayer,ncoda_2dvar_file_temp)

# Replace SSS on archv file
var_name = var_salt
fld2_salt_archv, fld2_salt_ncoda, nvar_salt = replace_field_in_archv_from_2DVAR_NCODA(artofs_file,var_name,klayer,ncoda_2dvar_file_salt)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
'''
#var_name = var_salt
#var_name = var_temp
#ncoda_2dvar_file = ncoda_2dvar_file_temp
var_name = var_salt
ncoda_2dvar_file = ncoda_2dvar_file_salt
print('Reading ab files')
rtofs_file = artofs_file[:-2] 
    
lines = [line.rstrip() for line in open(rtofs_file+'.b')]
ijdm = idm*jdm
npad = 4096-(ijdm%4096)
fld = ma.array([],fill_value=1.2676506002282294e+30)

inFile_a = rtofs_file + '.a'

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

fid = open(inFile_a,'rb')
fid.seek((nvar)*4*(npad+ijdm),0)
fld = fid.read(ijdm*4)
fld = struct.unpack('>'+str(ijdm)+'f',fld)
fld = np.array(fld)
fld_archv = ma.reshape(fld,(jdm,idm))

maxval_raw_archv = np.nanmax(fld_archv)
minval_raw_archv = np.nanmin(fld_archv)
print('from raw archv file ',maxval_raw_archv)
print('from raw archv file ',minval_raw_archv)

fld2_archv = np.copy(fld_archv)
mask = fld2_archv > 10**5
fld2_archv[mask] = np.nan

maxval_archv = np.nanmax(fld2_archv)
minval_archv = np.nanmin(fld2_archv)
print('from archv file ',maxval_archv)
print('from archv file ',minval_archv)

#############################
# Read ncoda 2DVAR file
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

#############################
# Replace value in bynary a file
binary_file = open(inFile_a,'r+b')
binary_file.seek((nvar)*4*(npad+ijdm))
binary_file.write(fld_ncoda_array_new_fill_value_bytes)
binary_file.close()

#############################
# Replace max and min value for temp at surface layer in b file
inFile_b = rtofs_file + '.b'
bfile = open(inFile_b,'r+')

linesb = open(inFile_b,'r+').readlines()

for n,line in enumerate(linesb):
    if len(line) > 0:
        if line.split()[0] == var_name and line.split()[4] == klayer:
            print(n)
            linesb[n] =  line[0:44] + '{:.7E}'.format(minval_ncoda) + '   ' + '{:.7E}'.format(maxval_ncoda) + '\n'

bfile.writelines(linesb)
'''
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#############################
inFile_a = artofs_file 
rtofs_file = artofs_file[:-2] 
lines = [line.rstrip() for line in open(rtofs_file+'.b')]
ijdm = idm*jdm
npad = 4096-(ijdm%4096)

binary_file = open(inFile_a,'rb')

# Check SST in archv file was modified
binary_file.seek((nvar_temp)*4*(npad+ijdm))
fld = binary_file.read(ijdm*4)
fld_unpack = struct.unpack('>'+str(ijdm)+'f',fld)
fld_array = np.array(fld_unpack)
fld_modified_archv = ma.reshape(fld_array,(jdm,idm))

fld2_temp_modified_archv = np.copy(fld_modified_archv)
mask = fld2_temp_modified_archv > 10**5
fld2_temp_modified_archv[mask] = np.nan

print('From modified archv file temp ',np.nanmax(fld_modified_archv))
print('From modified archv file temp ',np.nanmin(fld_modified_archv))

# Check SSS in archv file was modified
binary_file.seek((nvar_salt)*4*(npad+ijdm))
fld = binary_file.read(ijdm*4)
fld_unpack = struct.unpack('>'+str(ijdm)+'f',fld)
fld_array = np.array(fld_unpack)
fld_modified_archv = ma.reshape(fld_array,(jdm,idm))

fld2_salt_modified_archv = np.copy(fld_modified_archv)
mask = fld2_salt_modified_archv > 10**5
fld2_salt_modified_archv[mask] = np.nan

print('From modified archv file salt ',np.nanmax(fld_modified_archv))
print('From modified archv file salt ',np.nanmin(fld_modified_archv))

#############################################
min_val = 5
max_val = 32
delta_val = 2
units = '$^oC$'

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_temp_archv,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SST original archv ' + cycle ,fontsize=18)
fname = cycle + '_rtofs_sst.png'
fig.savefig(fname,bbox_inches='tight')

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_temp_ncoda,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SST 2dVAR NCODA ' + cycle ,fontsize=18)
fname = cycle + '_ncoda_sst.png'
fig.savefig(fname,bbox_inches='tight')

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_temp_modified_archv,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SST modified archv ' + cycle ,fontsize=18)
fname = cycle + '_rtofs_sst_modified.png'
fig.savefig(fname,bbox_inches='tight')

#############################################
min_val = 32
max_val = 38
delta_val = 0.5
units = ' '

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_salt_archv,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SSS original archv ' + cycle ,fontsize=18)
fname = cycle + '_rtofs_sss_f.png'
fig.savefig(fname,bbox_inches='tight')

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_salt_ncoda,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SSS 2dVAR NCODA ' + cycle ,fontsize=18)
fname = cycle + '_ncoda_sss.png'
fig.savefig(fname,bbox_inches='tight')

kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax = plt.subplots(figsize=(10, 5))
ax.set_facecolor("bisque")
ctr = ax.contourf(fld2_salt_modified_archv,cmap=colormap,**kw,extend='both')
cbar = fig.colorbar(ctr,pad=0.04,extendrect=True)
cbar.set_label(units,fontsize = 14)
ax.set_title('RTOFS SSS modified archv ' + cycle ,fontsize=18)
fname = cycle + '_rtofs_sss_modified.png'
fig.savefig(fname,bbox_inches='tight')

