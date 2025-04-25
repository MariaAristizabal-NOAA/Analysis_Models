#%% User input

tidal_amplitude_file = '/work/noaa/hwrf/save/maristiz/MOM6-NOAA-EMC-032123/ocean_only/RTOFS_2DVAR_SST_baseline/INPUT/ocean_tidal_amplitude.nc'

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
tidal_ampl_nc0 = nc.Dataset(tidal_amplitude_file)
tideamp = np.asarray(tidal_ampl_nc0['tideamp'][:])
tidal_ampl_nc0.close()

print(np.nanmin(tideamp))
print(np.nanmax(tideamp))

#################################################################################
plt.figure()
plt.contourf(tideamp)
plt.colorbar()

#################################################################################
nan_values = np.isnan(tideamp)
tideamp[nan_values] = 0.0

print(np.min(tideamp))
print(np.max(tideamp))

plt.figure()
plt.contourf(tideamp)
plt.colorbar()

# Write new tideamp with no nans to file
tidal_ampl_nc1 = nc.Dataset(tidal_amplitude_file,'a')
tidal_ampl_nc1['tideamp'][:,:] = tideamp
tidal_ampl_nc1.close()

#%% Reading modified file
#tidal_ampl_nc2 = xr.open_dataset(tidal_amplitude_file,decode_times=False)
tidal_ampl_nc2 = nc.Dataset(tidal_amplitude_file)
tideamp = np.asarray(tidal_ampl_nc2['tideamp'])
tidal_ampl_nc2.close()

print(np.min(tideamp))
print(np.max(tideamp))

