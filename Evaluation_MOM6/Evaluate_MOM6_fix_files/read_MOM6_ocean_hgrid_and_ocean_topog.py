#%% User input

hgrid_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/nhc/ocean_hgrid.nc'
topo_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/nhc/ocean_topog.nc'
#hgrid_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/jtnh/ocean_hgrid.nc'
#topo_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/jtnh/ocean_topog.nc'
#hgrid_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/jtsh/ocean_hgrid.nc'
#topo_file = '/work/noaa/hwrf/save/maristiz/hafsv2p1_final_20250409_Kd_output_MOM6/fix/fix_mom6/jtsh/ocean_topog.nc'

################################################################################
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Reading original file
#tidal_ampl_nc0 = xr.open_dataset(tidal_amplitude_file,decode_times=False)
ncgrid = nc.Dataset(hgrid_file)
xgrid = np.asarray(ncgrid['x'][:])
ygrid = np.asarray(ncgrid['y'][:])

nctopo = nc.Dataset(topo_file)
depth = np.asarray(nctopo['depth'][:])

#################################################################################
'''
plt.figure()
plt.plot(xgrid,ygrid,'.k')
plt.ylim([0,1])
plt.xlim([0,1])
'''

fig, ax = plt.subplots(figsize=(12, 6))
plt.contourf(xgrid[1::2,1::2],ygrid[1::2,1::2],depth,cmap=plt.cm.Spectral_r)
plt.colorbar()

plt.figure()
plt.pcolor(xgrid[1::2,1::2],ygrid[1::2,1::2],depth)
plt.colorbar()
plt.plot(xgrid,ygrid,'.k',markersize=0.5)
plt.ylim([7,8])
plt.xlim([-23,-21])

#################################################################################




