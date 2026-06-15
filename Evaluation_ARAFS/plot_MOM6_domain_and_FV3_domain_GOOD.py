#%% User input

hgrid_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out/ocean_hgrid_regional.nc'

topo_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out/ocean_topog.nc'

fv3_file = '/gpfs/f6/drsa-hurr1/world-shared/scrub/Maria.Aristizabal/ARAFS_alaska_coupled_ocean_HYCOM_cutout_ar-cpu_120h_update_ocn_prep_a/com/2023010600/00E/arafs.2023010600.f003.grb2'

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

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
# Read ocean grid, mask and depth 
ncgrid = nc.Dataset(hgrid_file)
xgrid = np.asarray(ncgrid['x'][:])
ygrid = np.asarray(ncgrid['y'][:])

nctopo = nc.Dataset(topo_file)
depth = np.asarray(nctopo['depth'][:])
wet = np.asarray(nctopo['wet'][:])

#################################################################################
# Read FV3 file
grb = grib2io.open(fv3_file,mode='r')

lat = grb.select(shortName='NLAT')[0].data
lon = grb.select(shortName='ELON')[0].data

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

print('Extracting MSLET')
slp = grb.select(shortName='MSLET')[0].data
slp = slp * 0.01 # convert Pa to hPa

#################################################################################

fig, ax = plt.subplots(figsize=(12, 6))
plt.pcolor(xgrid[1::2,1::2],ygrid[1::2,1::2],depth,cmap=plt.cm.Spectral_r)
plt.colorbar()

fig, ax = plt.subplots(figsize=(12, 6))
#plt.axis('scaled')
plt.pcolor(xgrid[1::2,1::2],ygrid[1::2,1::2],wet,cmap=plt.cm.coolwarm)
plt.colorbar()
plt.plot(xgrid[1::2,1::2][::30,::30],ygrid[1::2,1::2][::30,::30],color='grey')
plt.plot(xgrid[1::2,1::2][::30,::30].T,ygrid[1::2,1::2][::30,::30].T,color='grey')
cslevels = np.arange(840,1040,4)
cs = ax.contourf(lon-180, lat, slp, levels=cslevels,alpha=0.5)

#################################################################################

cartopy.config['data_dir'] = cartopyDataDir

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=myproj)
ax.axis('scaled')

ax.pcolor(xgrid[1::2,1::2]+180,ygrid[1::2,1::2],wet,cmap=plt.cm.Reds,transform=transform,alpha=0.5)
ax.plot(xgrid[::60,::60]+180,ygrid[::60,::60],color='grey',transform=transform,alpha=0.7)
ax.plot((xgrid[::60,::60].T)+180,ygrid[::60,::60].T,color='grey',transform=transform,alpha=0.7)

ax.plot(xgrid[:,-1]+180,ygrid[:,-1],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[:,0]+180,ygrid[:,0],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[-1,:]+180,ygrid[-1,:],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[0,:]+180,ygrid[0,:],color='red',linewidth=2,transform=transform,label='Ocean domain')

cs = ax.contourf(lon, lat, slp, levels=cslevels,alpha=0.5,transform=transform)
ax.plot(lon[::30,::30][0,-1],lat[::30,::30][0,-1],color='yellow',label='Atm. Domain')

plt.legend(ncol=2, loc='lower center', bbox_to_anchor=(0.5, -0.15),fontsize=18)

#ax.plot(lon[::30,::30][:,-1],lat[::30,::30][:,-1],color='yellow')
#ax.plot(lon[::30,::30][:,0],lat[::30,::30][:,0],color='yellow')
#ax.plot(lon[::30,::30][-1,:],lat[::30,::30][-1,:],color='yellow')
#ax.plot(lon[::30,::30][0,:],lat[::30,::30][0,:],color='yellow')

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}

ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)


