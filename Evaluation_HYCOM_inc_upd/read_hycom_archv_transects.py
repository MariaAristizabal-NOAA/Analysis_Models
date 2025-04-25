#%% User input

# RTOFS files 
#hycom_depth = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/hafs_hycom_hat10.basin.regional.depth_no_pacif'
#hycom_grid = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/regional.grid'
#afile = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2_HYCOM_increm_update/natl00l.2020072812.hafs_hycom_hat10.2020_210_18'
#afile = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/archv_mom6'
#afile = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/for_zulema/archv_rtofs'

hycom_depth = '/work/noaa/hwrf/noscrub/maristiz/RTOFS/fix/hafs_hycom_hat10.basin.regional.depth_99'
hycom_grid = '/work/noaa/hwrf/save/maristiz/hafs_develop_202112/fix/fix_hycom/hafs_hycom_hat10.basin.regional.grid'
afile = '/work/noaa/hwrf/scrub/maristiz/direct_increm_update_OBC3DF_mom6_minus_rtofs/com/2020072812/00L/natl00l.2020072812.hafs_hycom_hat10.2020_210_18'

#afile = 'archv_rtofs'
#afile = 'archv_inc_mom6'
#afile = 'incupd.2020_210_12'
#afile = 'diff_temp_mom6_minus_rtofs'
#afile = 'diff_saln_mom6_minus_rtofs'

'''
var_name = 'temp'
#colormap='nipy_spectral'
#colormap='seismic'
colormap='Spectral_r'
min_val = 0
max_val = 33
delta_val = 1  # delta in colorbar
delta_contour = 2 # delta in contour plot
units = '$^oC$'
lon_trans = True
lon_target = -60
lat_trans = True
lat_target = 24.5
ylim = [-4000,0]
'''


var_name = 'saln'
colormap='GnBu_r'
#colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = 30
max_val = 38.5
delta_val = 0.5
delta_contour = 0.5
units = ' '
lon_trans = True
lon_target = -50
lat_trans = True
lat_target = 26
ylim = [-500,0]


#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime
import os
import os.path
import glob
import numpy.ma as ma
import struct

import sys
sys.path.append('/home/Maria.Aristizabal/Repos/NCEP_scripts/')
from utils4HYCOM import readBinz, readgrids, parse_z
from utils4HYCOM import readdepth

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#### Copy Hycom-tools executables needed #################
os.system('cp /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HYCOM-tools/archive/src/archv2data3z .')
#os.system('cp /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/increm_mom6_minus_rtofs/archv2data3z_rtofs.in .')

#%% Link regional.depth and regional.grid files ##########
os.system('ln -sf ' + hycom_depth + '.a regional.depth.a')
os.system('ln -sf ' + hycom_depth + '.b regional.depth.b')
os.system('ln -sf ' + hycom_grid + '.a regional.grid.a')
os.system('ln -sf ' + hycom_grid + '.b regional.grid.b')

#%% Link archive file ##########
os.system('ln -sf ' + afile + '.a archv_rtofs.a')
os.system('ln -sf ' + afile + '.b archv_rtofs.b')

#%% Convert archive file to fixed z levels #########
#### and extracts different fields ##############
os.system('./archv2data3z < archv2data3z_rtofs_deeper.in')
os.system('mv fort.010a archv_zlevels_temp_rtofs.a')
os.system('mv fort.10 archv_zlevels_temp_rtofs.b')
os.system('mv fort.011a archv_zlevels_saln_rtofs.a')
os.system('mv fort.11 archv_zlevels_saln_rtofs.b')
os.system('mv fort.012a archv_zlevels_u_rtofs.a')
os.system('mv fort.12 archv_zlevels_u_rtofs.b')
os.system('mv fort.013a archv_zlevels_v_rtofs.a')
os.system('mv fort.13 archv_zlevels_v_rtofs.b')
os.system('mv fort.014a archv_zlevels_dens_rtofs.a')
os.system('mv fort.14 archv_zlevels_dens_rtofs.b')

#%% Reading hycom grid
print('Retrieving coordinates from RTOFS')
# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid+'.b')]
idm = int([line.split() for line in lines_grid if 'idm' in line][0][0])
jdm = int([line.split() for line in lines_grid if 'jdm' in line][0][0])
lon_hycom = np.array(readgrids(hycom_grid,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid,'plat:',[0]))
depth_hycom = np.asarray(readdepth(hycom_depth,'depth'))

#%% Getting grid points for transects
oklon = int(np.interp(lon_target +360,lon_hycom[0,:],np.arange(len(lon_hycom[0,:]))))
oklat = int(np.interp(lat_target,lat_hycom[:,0],np.arange(len(lat_hycom[:,0]))))

####################################################################
lines = [line.rstrip() for line in open(afile+'.b')]
ijdm = idm*jdm
npad = 4096-(ijdm%4096)
fld = ma.array([],fill_value=1.2676506002282294e+30)

inFile = afile + '.a'
fid = open(inFile,'rb')

for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == 'field':
            nheading = n + 1
            print(line.split()[0])

nvar = []
dens_layer = []
for n,line in enumerate(lines):
    if len(line) > 0:
        if line.split()[0] == var_name:
            nvar.append(n - nheading)
            dens_layer.append(line.split()[5])
nk = len(nvar)
dens_layer = np.asarray(dens_layer)

if var_name == 'temp':
    infile = 'archv_zlevels_temp_rtofs'
if var_name == 'saln':
    infile = 'archv_zlevels_saln_rtofs'
if var_name == 'u':
    infile = 'archv_zlevels_u_rtofs'
if var_name == 'v':
    infile = 'archv_zlevels_v_rtofs'

lines = [line.rstrip() for line in open(infile+'.b')]
depths = []
for n,line in enumerate(lines):
    if len(line) > 0:
        depths.append(float(lines[n].split('=')[1].split()[0]))
depths = np.asarray(depths)

fid = open(infile+'.a','rb')

if lon_trans:
    #%% Getting grid points for transects
    oklat = int(np.interp(lat_target,lat_hycom[:,0],np.arange(len(lat_hycom[:,0]))))
    trans_zonal = np.empty([len(depths),int(idm)])
    trans_zonal[:] = np.nan

if lat_trans:
    #%% Getting grid points for transects
    oklon = int(np.interp(lon_target +360,lon_hycom[0,:],np.arange(len(lon_hycom[0,:]))))
    trans_merid = np.empty([len(depths),int(jdm)])
    trans_merid[:] = np.nan

for d in np.arange(len(depths)):
    print(d)
    fid.seek((d)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    fld = np.array(fld)
    fld = ma.reshape(fld,(jdm,idm))
    fld2 = np.copy(fld)
    mask = fld2 > 10**5
    fld2[mask] = np.nan
    maxval = np.nanmax(fld2)
    minval = np.nanmin(fld2)
    print(maxval)
    print(minval)
    if lon_trans:
        trans_zonal[d,:] = fld2[oklat,:]
    if lat_trans:
        trans_merid[d,:] = fld2[:,oklon]

if lon_trans:
    kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
    fig,ax1 = plt.subplots(figsize = (11,6))
    ax1.set_facecolor("bisque")
    ctr = ax1.contour(lon_hycom[0,:]-360,-depths,trans_zonal,np.arange(min_val,max_val+0.1,delta_contour),colors='k')
    ax1.clabel(ctr,ctr.levels,inline=True,fontsize=12)
    cl = ax1.contourf(lon_hycom[0,:]-360,-depths,trans_zonal,cmap=colormap,**kw)
    cbar = fig.colorbar(cl, fraction=0.025, pad=0.04)
    cbar.set_label(units,fontsize = 14)
    ax1.set_title(var_name + ' Zonal Transect ',fontsize=18)
    ax1.set_ylim(ylim)

if lat_trans:
    kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
    fig,ax1 = plt.subplots(figsize = (11,6))
    ax1.set_facecolor("bisque")
    ctr = ax1.contour(lat_hycom[:,0],-depths,trans_merid,np.arange(min_val,max_val+0.1,delta_contour),colors='k')
    ax1.clabel(ctr,ctr.levels,inline=True,fontsize=12)
    cl = ax1.contourf(lat_hycom[:,0],-depths,trans_merid,cmap=colormap,**kw)
    cbar = fig.colorbar(cl, fraction=0.025, pad=0.04)
    cbar.set_label(units,fontsize = 14)
    ax1.set_title(var_name + ' Meridional Transect ',fontsize=18)
    ax1.set_ylim(ylim)


####################################################################

