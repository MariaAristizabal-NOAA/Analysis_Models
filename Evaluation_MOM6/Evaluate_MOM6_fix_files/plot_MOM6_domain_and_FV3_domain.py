#%% User input

hgrid_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out/ocean_hgrid_regional.nc'

topo_file = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Scripts_to_create_MOM6_domain_fix_files_HYCOM_cut_out/ocean_topog.nc'

fv3_file = '/gpfs/f6/drsa-hurr1/world-shared/scrub/Maria.Aristizabal/ARAFS_alaska_coupled_ocean_ar-cpu_120h_a/com/2025121100/00E/arafs.2025121100.f003.grb2'

################################################################################
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import grib2io

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
'''
print('raw lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
if abs(np.max(lon) - 360.) < 10.:
    lon[lon>180] = lon[lon>180] - 360.
    lon_offset = 0.
else:
    lon_offset = 180.
lon = lon - lon_offset
print('new lonlat limit: ', np.min(lon), np.max(lon), np.min(lat), np.max(lat))
'''

lon_offset = 360
lon = lon - lon_offset

print('Extracting MSLET')
slp = grb.select(shortName='MSLET')[0].data
slp = slp * 0.01 # convert Pa to hPa

#################################################################################

fig, ax = plt.subplots(figsize=(12, 6))
plt.pcolor(xgrid[1::2,1::2],ygrid[1::2,1::2],depth,cmap=plt.cm.Spectral_r)
plt.colorbar()

fig, ax = plt.subplots(figsize=(12, 6))
plt.pcolor(xgrid[1::2,1::2],ygrid[1::2,1::2],wet,cmap=plt.cm.coolwarm)
plt.colorbar()

cslevels = np.arange(840,1040,4)
cs = ax.contourf(lon, lat, slp, levels=cslevels,alpha=0.5)
#cs = ax.contour(lon, lat, slp, levels=cslevels, colors='black', linewidths=0.6)

#################################################################################




