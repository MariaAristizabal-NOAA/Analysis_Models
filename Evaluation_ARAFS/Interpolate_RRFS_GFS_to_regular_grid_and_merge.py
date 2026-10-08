#%% User input

gfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/COMGFSv16/gfs.20260902/00/atmos/gfs.t00z.pgrb2.0p25.f000'

rrfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.20260902/rrfs.t00z.prslev.13km.f000.na.nc'
#rrfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.t00z.prslev.13km.f000.na.grib2'

cartopyDataDir = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/local/share/cartopy'

lon_ini = -215
lon_end = -69
lat_ini = -21
lat_end = 80

# Resolution in degrees of regular grid
dx = 1/12  # Longitude resolution
dy = 1/12  # Latitude resolution

var_list = ['TMP','RH','UGRD','VGRD','HGT','PRES','SPFH','UGRD']
level_list = ['2 mb','5 mb','7 mb', '10 mb','20 mb','30 mb','50 mb','70 mb','100 mb','150 mb','200 mb','250 mb','300 mb','350 mb','400 mb', '450 mb','500 mb','550 mb','600 mb','650 mb','700 mb','750 mb','800 mb', '850 mb','900 mb','925 mb','950 mb','975 mb','1000 mb']

################################################################################
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import grib2io
import time

from collections import defaultdict

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.interpolate import griddata

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#####################################################################
# Read lon and lat

# GFS
gfs = grib2io.open(gfs_file,mode='r')
lat_gfs = gfs.select(shortName=var_list[0])[0].lats
lon_gfs = gfs.select(shortName=var_list[0])[0].lons

# RRFS
rrfs = nc.Dataset(rrfs_file)
lon_rrfs = np.asarray(rrfs['longitude'][:])
lat_rrfs = np.asarray(rrfs['latitude'][:])


####################################################################
# Create regular grid
# Calculate number of points needed for exact end bounds
num_lons = int(round((lon_end - (lon_ini)) / dx)) + 1
num_lats = int(round((lat_end - (lat_ini)) / dy)) + 1

lon_1d = np.linspace(lon_ini, lon_end, num_lons)
lat_1d = np.linspace(lat_ini, lat_end, num_lats)

# 2. Generate 2D Regular Grid
lon_reg, lat_reg = np.meshgrid(lon_1d, lat_1d)

#####################################################################
for var in var_list[0:1]:
    for level in level_list[0:1]:
        start_time = time.perf_counter()
        # Read GFS file
        var_gfs = gfs.select(shortName=var,level=level)[0].data
        
        #############################################################
        # Read RRFS grid
        var_rrfs = rrfs[var+'_'+"".join(level.split())][0,:,:]
        
        ######################################################################
        # Interpolate from RRFS to regular grid
        # Flatten Source Coordinates and Values to 1D
        points_src = (lon_rrfs.ravel()-360,lat_rrfs.ravel())
        values_src = var_rrfs.ravel()
        
        # Interpolate directly using griddata
        var_rrfs_interp = griddata(
            points=points_src,
            values=values_src,
            xi=(lon_reg, lat_reg),
            method='linear'  # Options: 'linear', 'nearest', 'cubic'
        )
        
        var_rrfs_interp[var_rrfs_interp>1000] = np.nan
        
        #####################################################################
        # Interpolate from GFS to regular grid
        # Flatten Source Coordinates and Values to 1D
        points_src = (lon_gfs.ravel()-360,lat_gfs.ravel())
        values_src = var_gfs.ravel()
        
        # Interpolate directly using griddata
        var_gfs_interp = griddata(
            points=points_src,
            values=values_src,
            xi=(lon_reg, lat_reg),
            method='linear'  # Options: 'linear', 'nearest', 'cubic'
        )
        
        var_gfs_interp[var_gfs_interp>1000] = np.nan
        
        #####################################################################
        # Merge GFS fields with RRFS fields on the same regular grid
        empty = np.isnan(var_rrfs_interp)
        var_rrfs_interp_merged = np.copy(var_rrfs_interp)
        var_rrfs_interp_merged[empty] = var_gfs_interp[empty]

        varr = gfs.select(shortName=var,level=level)[0]
        varr.data = var_rrfs_interp_merged
        varr.ny = num_lats
        varr.nx = num_lons

        rrfs_merged_file = 'rrfs.t00z.prslev.merged.gfs.9km.f000'
        with grib2io.open(rrfs_merged_file, mode='w') as fgrib2:
            print('Creating grib2 file '+rrfs_merged_file)
            fgrib2.write(varr)

        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"Execution time: {execution_time:.6f} seconds")

#####################################################################
'''
# GFS and RRFS fields merged on the same domain on regular grid
cartopy.config['data_dir'] = cartopyDataDir

lon_offset = 180
myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=myproj)
ax.axis('scaled')

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

cs = ax.contourf(lon_reg+180, lat_reg, tmp_rrfs_interp_merged,transform=transform)
fig.colorbar(cs)

ax.plot(lon_reg[:,-1]+180,lat_reg[:,-1],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[:,0]+180,lat_reg[:,0],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[-1,:]+180,lat_reg[-1,:],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[0,:]+180,lat_reg[0,:],color='cyan',linewidth=2,transform=transform)

ax.plot(lon_reg[::50,::50]+180, lat_reg[::50,::50],alpha=0.2,color='k',transform=transform)
ax.plot(lon_reg[::50,::50].T+180, lat_reg[::50,::50].T,alpha=0.2,color='k',transform=transform)
'''
