#%% User input

gfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/COMGFSv16/gfs.20260401/00/atmos/gfs.t00z.pgrb2.0p25.f000'

rrfs_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.t00z.prslev.13km.f000.na.nc'

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
# Read GFS file
grb = grib2io.open(gfs_file,mode='r')

lat = grb.select(shortName='TMP')[0].lats
lon = grb.select(shortName='TMP')[0].lons

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

################################################################################
# Read RRFS grid
rrfs_grid = nc.Dataset(rrfs_file)
lon_rrfs = np.asarray(rrfs_grid['longitude'][:])
lat_rrfs = np.asarray(rrfs_grid['latitude'][:])
#field_rrfs = rrfs_grid['HGT_2mb'][0,:,:]
field_rrfs = rrfs_grid['TMP_1000mb'][0,:,:]

#################################################################################
# RRFS domain
cartopy.config['data_dir'] = cartopyDataDir

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=myproj)
ax.axis('scaled')

cs = ax.contourf(lon_rrfs-180, lat_rrfs, field_rrfs,transform=transform)
ax.plot(lon_rrfs[::30,::30][0,-1],lat_rrfs[::30,::30][0,-1],color='green',label='RRFS Domain')
fig.colorbar(cs)

ax.plot(lon_rrfs[::10,::10]-180, lat_rrfs[::10,::10], field_rrfs,alpha=0.5,transform=transform)
ax.plot(lon_rrfs[::10,::10].T-180, lat_rrfs[::10,::10].T, field_rrfs,alpha=0.5,transform=transform)

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.1', alpha=0.6, linestyle=(0, (5, 10)))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'black'}
gl.ylabel_style = {'size': 8, 'color': 'black'}


#################################################################################

cartopy.config['data_dir'] = cartopyDataDir

myproj = ccrs.PlateCarree(lon_offset)
transform = ccrs.PlateCarree(lon_offset)

# create figure and axes instances
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=myproj)
ax.axis('scaled')

'''
ax.pcolor(xgrid[1::2,1::2]+180,ygrid[1::2,1::2],wet,cmap=plt.cm.Reds,transform=transform,alpha=0.5)

ax.plot(xgrid[:,-1]+180,ygrid[:,-1],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[:,0]+180,ygrid[:,0],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[-1,:]+180,ygrid[-1,:],color='red',linewidth=2,transform=transform)
ax.plot(xgrid[0,:]+180,ygrid[0,:],color='red',linewidth=2,transform=transform,label='Ocean domain')
'''

#cs = ax.contourf(lon, lat, slp, levels=cslevels,alpha=0.5,transform=transform)
cs = ax.contourf(lon, lat, slp, alpha=0.5,transform=transform)
ax.plot(lon[::30,::30][0,-1],lat[::30,::30][0,-1],color='yellow',label='Atm. Domain')

cs = ax.contourf(lon_rrfs-180, lat_rrfs, field_rrfs,alpha=0.5,transform=transform)

ax.plot(lon_rrfs[::30,::30][0,-1],lat_rrfs[::30,::30][0,-1],color='green',label='RRFS Domain')
#fig.colorbar(cs)
#ax.plot(lon[::30,::30][0,-1],lat[::30,::30][0,-1],color='yellow',label='Atm. Domain')

plt.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, -0.15),fontsize=18)

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

#ax.set_extent([xlim[0]+lon_offset, xlim[1]+lon_offset, ylim[0], ylim[1]], crs=transform)

#plt.savefig("/ncrc/home1/Maria.Aristizabal/Figures_coupled_AR_AFS_paper/MOM6_fv3_domains.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.02)
