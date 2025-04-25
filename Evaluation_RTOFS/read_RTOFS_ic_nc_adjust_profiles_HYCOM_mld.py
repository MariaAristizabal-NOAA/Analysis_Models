#! /usr/bin/env python3

#%% User input
'''
# forecasting cycle to be used
cycle = '2024011506'

rtofs_ic_file = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_IC/2024011506_SST_2DVAR_v2_NHC_test2/ocean_ts_ic.nc'
'''

ymdh = '2024011506'
rtofs_ic_netcdf_file = '/work2/noaa/hwrf/scrub/maristiz/hafsv2_SST_2DVAR_mld_archv_hnod/2024041506/00L/ocn_prep/mom6_init/ocean_ts_ic.nc'
rtofs_mld_netcdf_file = '/work2/noaa/hwrf/scrub/maristiz/hafsv2_SST_2DVAR_mld_archv_hnod/2024041506/00L/ocn_prep/mom6_init/ocean_ssh_ic.nc'
depth_var_name = 'depth'
temp_var_name = 'Temp'
salt_var_name = 'Salt'
mld_var_name = 'mld'

################################################################################
import argparse
import xarray as xr
import netCDF4 as nc
import numpy as np

################################################################################
'''
if __name__ == "__main__":
    # get command line args
    parser = argparse.ArgumentParser(
        description="Modify the temperature and salinity profile from the RTOFS initial condition file (netcdf format) so the SST and SSS from the 2DVAR NCODA analysis is merged with the original temperature profile")
    parser.add_argument('ymdh', type=str, help="HAFS 6-hourly cycle. Ex: 2024011506")
    parser.add_argument('rtofs_ic_netcdf_file', type=str, help="Name of the RTOFS initial condition file for temp. and salinity. Ex: ocean_ts_ic.nc")
    parser.add_argument('rtofs_mld_netcdf_file', type=str, help="Name of the RTOFS file that contains the mixed layer depth variable. Ex: ocean_ssh_ic.nc")
    parser.add_argument('depth_var_name', type=str, help="Name of the depth variable in the netcdf file. Ex: 'depth'")
    parser.add_argument('temp_var_name', type=str, help="Name of the temperature variable in the netcdf file. Ex: 'Temp'")
    parser.add_argument('salt_var_name', type=str, help="Name of the salinity variable in the netcdf file. Ex: 'Salt'")
    parser.add_argument('mld_var_name', type=str, help="Name of the mixed layer depth variable in the netcdf file. Ex: 'mld'")

    args = parser.parse_args()

    ymdh = args.ymdh
    rtofs_ic_netcdf_file = args.rtofs_ic_netcdf_file
    rtofs_mld_netcdf_file = args.rtofs_mld_netcdf_file
    depth_var_name = args.depth_var_name
    temp_var_name = args.temp_var_name
    salt_var_name = args.salt_var_name
    mld_var_name = args.mld_var_name
'''
################################################################################
#%% Reading grid
#rtofs_ic = xr.open_dataset(rtofs_ic_netcdf_file,decode_times=False)
rtofs_ic = nc.Dataset(rtofs_ic_netcdf_file,decode_times=False)
depth = np.asarray(rtofs_ic[depth_var_name][:])
temp = np.asarray(rtofs_ic[temp_var_name][0,:,:,:])
salt = np.asarray(rtofs_ic[salt_var_name][0,:,:,:])
rtofs_mld = nc.Dataset(rtofs_mld_netcdf_file,decode_times=False)
mld = np.asarray(rtofs_mld[mld_var_name][0,:,:])

rtofs_ic.close()
rtofs_mld.close()

#################################################################################
# Adjust profiles
fill_value = 1.2676506002282294e+30
adj_temp_prof = np.empty((len(depth),temp.shape[1],temp.shape[2]))
adj_temp_prof[:] = fill_value
adj_salt_prof = np.empty((len(depth),temp.shape[1],temp.shape[2]))
adj_salt_prof[:] = fill_value
for y in np.arange(temp.shape[1]):
    print(y)
    for x in np.arange(temp.shape[2]):
        temp_prof = temp[:,y,x]
        salt_prof = salt[:,y,x]
        mldd = mld[y,x] 
        if mldd != fill_value:
            inside_mld = np.where(depth <= mldd)[0]

            adj_temp_pr = np.copy(temp_prof)
            adj_temp_pr[inside_mld] = temp_prof[0]
            adj_temp_prof[:,y,x] = adj_temp_pr

            adj_salt_pr = np.copy(salt_prof)
            adj_salt_pr[inside_mld] = salt_prof[0]
            adj_salt_prof[:,y,x] = adj_salt_pr
        else:
            adj_temp_prof[:,y,x] = temp_prof
            adj_salt_prof[:,y,x] = salt_prof

#################################################################################
# Write new adjusted temp profiles to file
'''read_RTOFS_ic_nc_adjust_profiles_HYCOM_mld.py
nc_file = nc.Dataset(rtofs_ic_netcdf_file,'a')
nc_file[temp_var_name][0,:,:,:] = adj_temp_prof
nc_file[salt_var_name][0,:,:,:] = adj_salt_prof
nc_file.close()
'''
