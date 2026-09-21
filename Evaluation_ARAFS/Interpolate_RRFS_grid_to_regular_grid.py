#%% User input

hgrid_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out_SUCCESFUL/ocean_hgrid_regional.nc'

topo_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out_SUCCESFUL/ocean_topog.nc'

fv3_file = '/gpfs/f6/drsa-hurr1/world-shared/scrub/Maria.Aristizabal/ARAFS_alaska_coupled_ocean_HYCOM_cutout_ar-cpu_120h_update_ocn_prep_a/com/2023010600/00E/arafs.2023010600.f003.grb2'

gfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/COMGFSv16/gfs.20260902/00/atmos/gfs.t00z.pgrb2.0p25.f000'

rrfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.20260902/rrfs.t00z.prslev.13km.f000.na.nc'
#rrfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.t00z.prslev.13km.f000.na.grib2'

cartopyDataDir = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/local/share/cartopy'

xlim = [-240,-40]
ylim = [-30,85]

################################################################################
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import grib2io

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from scipy.interpolate import griddata

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
# Read ocean grid, mask and depth
ncgrid = nc.Dataset(hgrid_file)
xgrid = np.asarray(ncgrid['x'][:])
ygrid = np.asarray(ncgrid['y'][:])

xgrid_360 = np.mod(xgrid, 360)

nctopo = nc.Dataset(topo_file)
depth = np.asarray(nctopo['depth'][:])
wet = np.asarray(nctopo['wet'][:])

################################################################################
# Read FV3 file
grb = grib2io.open(fv3_file,mode='r')
lat_arafs = grb.select(shortName='NLAT')[0].data
lon_arafs = grb.select(shortName='ELON')[0].data
tmp_arafs = grb.select(shortName='TMP',level='1000 mb')[0].data

lon_arafs_360 = np.mod(lon_arafs, 360)

'''
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
'''

################################################################################
# Read GFS file
grb = grib2io.open(gfs_file,mode='r')
lat_gfs = grb.select(shortName='TMP')[0].lats
lon_gfs = grb.select(shortName='TMP')[0].lons
tmp_gfs = grb.select(shortName='TMP',level='1000 mb')[0].data

lon_gfs_360 = np.mod(lon_gfs, 360)

'''
# The lon range in grib2 is typically between 0 and 360
# Cartopy's PlateCarree projection typically uses the lon range of -180 to 180
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
'''

################################################################################
# Read RRFS grid
rrfs_grid = nc.Dataset(rrfs_file)
lon_rrfs = np.asarray(rrfs_grid['longitude'][:])
lat_rrfs = np.asarray(rrfs_grid['latitude'][:])
tmp_rrfs = rrfs_grid['TMP_1000mb'][0,:,:]

'''
grb = grib2io.open(rrfs_file,mode='r')
lat_rrfs = grb.select(shortName='TMP')[0].lats
lon_rrfs = grb.select(shortName='TMP')[0].lons
tmp_rrfs = grb.select(shortName='TMP',level='1000 mb')[0].data

flat_indices = np.argsort(lon_rrfs, axis=None)

lon_rrfs_sort = np.reshape(np.ravel(lon_rrfs)[flat_indices],[lon_rrfs.shape[0],lon_rrfs.shape[1]])
lat_rrfs_sort = np.reshape(np.ravel(lat_rrfs)[flat_indices],[lon_rrfs.shape[0],lon_rrfs.shape[1]])
tmp_rrfs_sort = np.reshape(np.ravel(tmp_rrfs)[flat_indices],[lon_rrfs.shape[0],lon_rrfs.shape[1]])
'''

'''
plt.figure()
plt.contourf(lon_rrfs_sort, lat_rrfs_sort, tmp_rrfs_sort)
'''

################################################################################
# Create a regular grid 
# 1. Define 1D coordinate vectors
# Set resolution in degrees
dx = 1/12  # Longitude resolution
dy = 1/12  # Latitude resolution

# Calculate number of points needed for exact end bounds
num_lons = int(round((-69 - (-215)) / dx)) + 1
num_lats = int(round((80 - (-21)) / dy)) + 1

lon_1d = np.linspace(-215, -69, num_lons)
lat_1d = np.linspace(-21, 80, num_lats)

# 2. Generate 2D Regular Grid
lon_reg, lat_reg = np.meshgrid(lon_1d, lat_1d)

# --- Verification ---
print(f"1D Longitude shape: {lon_1d.shape} (Range: {lon_1d.min()} to {lon_1d.max()})")
print(f"1D Latitude shape:  {lat_1d.shape} (Range: {lat_1d.min()} to {lat_1d.max()})")
print(f"2D Grid shape (Ny, Nx): {lon_reg.shape}")

##############################################################################
# Interpolate from RRFS to regular grid
# Flatten Source Coordinates and Values to 1D
points_src = (lon_rrfs.ravel()-360,lat_rrfs.ravel())
values_src = tmp_rrfs.ravel()

# Interpolate directly using griddata
tmp_rrfs_interp = griddata(
    points=points_src,
    values=values_src,
    xi=(lon_reg, lat_reg),
    method='linear'  # Options: 'linear', 'nearest', 'cubic'
)

tmp_rrfs_interp[tmp_rrfs_interp>1000] = np.nan

##############################################################################
# Interpolate from GFS to regular grid
# Flatten Source Coordinates and Values to 1D
points_src = (lon_gfs.ravel()-360,lat_gfs.ravel())
values_src = tmp_gfs.ravel()

# Interpolate directly using griddata
tmp_gfs_interp = griddata(
    points=points_src,
    values=values_src,
    xi=(lon_reg, lat_reg),
    method='linear'  # Options: 'linear', 'nearest', 'cubic'
)

tmp_gfs_interp[tmp_gfs_interp>1000] = np.nan

##############################################################################
# Merge GFS fields with RRFS fields on the same regular grid
empty = np.isnan(tmp_rrfs_interp)
tmp_rrfs_interp_merged = np.copy(tmp_rrfs_interp)
tmp_rrfs_interp_merged[empty] = tmp_gfs_interp[empty]

###############################################################################
# RRFS domain
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

cs = ax.contourf(lon_rrfs-180, lat_rrfs, tmp_rrfs,transform=transform)
ax.plot(lon_rrfs[::30,::30][0,-1]+180,lat_rrfs[::30,::30][0,-1],color='green',label='RRFS Domain')
fig.colorbar(cs)

ax.plot(lon_rrfs[::10,::10]-180, lat_rrfs[::10,::10],alpha=0.5,transform=transform)
ax.plot(lon_rrfs[::10,::10].T-180, lat_rrfs[::10,::10].T,alpha=0.5,transform=transform)

###############################################################################
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

cs = ax.contourf(lon_rrfs-180, lat_rrfs, tmp_rrfs,alpha=0.8,transform=transform,cmap='Oranges')
ax.plot(lon_rrfs[::30,::30][0,-1]-180,lat_rrfs[::30,::30][0,-1],color='orange',label='RRFS Domain')

cs = ax.contourf(lon_gfs-180, lat_gfs, tmp_gfs, alpha=0.5,cmap='Greens',transform=transform)
ax.plot(lon_gfs[::30,::30][0,-1],lat_gfs[::30,::30][0,-1],color='green',label='GFS Domain')

ax.plot(xgrid[:,-1]+180,ygrid[:,-1],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[:,0]+180,ygrid[:,0],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[-1,:]+180,ygrid[-1,:],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[0,:]+180,ygrid[0,:],color='red',linewidth=2,transform=transform,label='Ocean domain')

cs = ax.contourf(lon_arafs-180, lat_arafs, tmp_arafs, alpha=0.8,transform=transform,cmap='Blues')
ax.plot(lon_arafs[::30,::30][0,-1]-180,lat_arafs[::30,::30][0,-1],color='blue',label='AR-AFS Domain')

ax.plot(lon_reg[:,-1]+180,lat_reg[:,-1],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[:,0]+180,lat_reg[:,0],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[-1,:]+180,lat_reg[-1,:],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[0,:]+180,lat_reg[0,:],color='cyan',linewidth=2,transform=transform)

plt.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.15),fontsize=14)

#ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

#plt.savefig("/ncrc/home1/Maria.Aristizabal/Figures_coupled_AR_AFS_paper/MOM6_fv3_domains.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.02)

###############################################################################
# RRFS domain onn original grid
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

cs = ax.contourf(lon_rrfs-180, lat_rrfs, tmp_rrfs,transform=transform)
ax.plot(lon_rrfs[::30,::30][0,-1]+180,lat_rrfs[::30,::30][0,-1],color='green',label='RRFS Domain')
fig.colorbar(cs)

ax.plot(lon_rrfs[::10,::10]-180, lat_rrfs[::10,::10],alpha=0.1,transform=transform)
ax.plot(lon_rrfs[::10,::10].T-180, lat_rrfs[::10,::10].T,alpha=0.1,transform=transform)

###############################################################################
# RRFS domain on regular grid
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

cs = ax.contourf(lon_reg+180, lat_reg, tmp_rrfs_interp,transform=transform)
fig.colorbar(cs)

ax.plot(lon_reg[:,-1]+180,lat_reg[:,-1],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[:,0]+180,lat_reg[:,0],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[-1,:]+180,lat_reg[-1,:],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[0,:]+180,lat_reg[0,:],color='cyan',linewidth=2,transform=transform)

ax.plot(lon_reg[::50,::50]+180, lat_reg[::50,::50],alpha=0.2,color='k',transform=transform)
ax.plot(lon_reg[::50,::50].T+180, lat_reg[::50,::50].T,alpha=0.2,color='k',transform=transform)

###############################################################################
# GFS field on original grid
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

cs = ax.contourf(lon_gfs-180, lat_gfs, tmp_gfs, alpha=0.5,transform=transform)

###############################################################################
# GFS domain on regular grid
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

levels = np.arange(290,318)
cs = ax.contourf(lon_reg+180, lat_reg, tmp_gfs_interp,transform=transform,levels=levels)
fig.colorbar(cs)

ax.plot(lon_reg[:,-1]+180,lat_reg[:,-1],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[:,0]+180,lat_reg[:,0],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[-1,:]+180,lat_reg[-1,:],color='cyan',linewidth=2,transform=transform)
ax.plot(lon_reg[0,:]+180,lat_reg[0,:],color='cyan',linewidth=2,transform=transform)

ax.plot(lon_reg[::50,::50]+180, lat_reg[::50,::50],alpha=0.2,color='k',transform=transform)
ax.plot(lon_reg[::50,::50].T+180, lat_reg[::50,::50].T,alpha=0.2,color='k',transform=transform)

###############################################################################
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

