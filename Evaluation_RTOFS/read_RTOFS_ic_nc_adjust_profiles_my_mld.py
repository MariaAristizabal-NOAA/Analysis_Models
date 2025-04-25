#! /usr/bin/env python3

#%% User input
'''
# forecasting cycle to be used
cycle = '2024011506'

rtofs_ic_file = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_IC/2024011506_SST_2DVAR_v2_NHC_test2/ocean_ts_ic.nc'
'''

#ymdh = 2024041506
#rtofs_ic_netcdf_file = '/work2/noaa/hwrf/scrub/maristiz/hafsv2_SST_2DVAR_hnod/2024041506/00L/ocn_prep/mom6_init/ocean_ts_ic.nc'
ymdh = 2020082506
rtofs_ic_netcdf_file = '/work2/noaa/hwrf/scrub/maristiz/hafsv2_h3a0/2020082506/13L/ocn_prep/mom6_init/ocean_ts_ic.nc'
depth_var_name = 'depth'
temp_var_name = 'Temp'
salt_var_name = 'Salt'

################################################################################
import argparse
import xarray as xr
import netCDF4 as nc
import numpy as np

from eos80 import dens

################################################################################
def MLD_dens_crit(drho,ref_depth,depth,temp,dens):
    # This function calculates the mixed layer depth and Mixed layer temperature
    # based on a density criteria: rho_at_ref_depth - rho <= drho
    # Inputs
    # drho: delta density from the mixed layer depth definition used
    # ref_depth: Reference depth from the mixed layer depth definition used
    # depth, temp and dens: 1D vectors depth, temperature and density
    # Output
    # MLD and MLT: mixed layer depth and Mixed layer temperature

    ok_ref_depth = np.where(depth >= ref_depth)[0][0]
    rho_ref_depth = dens[ok_ref_depth]
    delta_rho = -(rho_ref_depth - dens)
    ok_mld_rho = np.where(delta_rho <= drho)[0]

    if ok_mld_rho.size == 0:
        MLD = np.nan
        MLT = np.nan
    else:
        MLD = depth[ok_mld_rho[-1]]
        MLT = np.nanmean(temp[ok_mld_rho])

    return MLD, MLT

################################################################################
def MLD_temp_crit(dtemp,ref_depth,depth,temp):
    # This function calculates the mixed layer depth and Mixed layer temperature
    # based on a temperature criteria: T - T_at_ref_depth <= dtemp
    # Inputs
    # dtemp: delta temperature from the mixed layer depth definition used
    # ref_depth: Reference depth from the mixed layer depth definition used
    # depth and temp: 1D vectors depth and temperature
    # Output
    # MLD and MLT: mixed layer depth and Mixed layer temperature

    ok_ref_depth = np.where(depth >= ref_depth)[0][0]
    temp_ref_depth = temp[ok_ref_depth]
    delta_T = temp_ref_depth - temp
    ok_mld_temp = np.where(delta_T <= dtemp)[0]

    if ok_mld_temp.size == 0:
        MLD = np.nan
        MLT = np.nan
    else:
        MLD = depth[ok_mld_temp[-1]]
        MLT = np.nanmean(temp[ok_mld_temp])

    return MLD, MLT

################################################################################
'''
if __name__ == "__main__":
    # get command line args
    parser = argparse.ArgumentParser(
        description="Modify the temperature and salinity profile from the RTOFS initial condition file (netcdf format) so the SST and SSS from the 2DVAR NCODA analysis is merged with the original temperature profile")
    parser.add_argument('ymdh', type=str, help="HAFS 6-hourly cycle. Ex: 2024011506")
    parser.add_argument('rtofs_ic_netcdf_file', type=str, help="Name of the RTOFS initial condition file for temp. and salinity. Ex: ocean_ts_ic.nc")
    parser.add_argument('depth_var_name', type=str, help="Name of the depth variable in the netcdf file. Ex: 'depth'")
    parser.add_argument('temp_var_name', type=str, help="Name of the temperature variable in the netcdf file. Ex: 'Temp'")
    parser.add_argument('salt_var_name', type=str, help="Name of the salinity variable in the netcdf file. Ex: 'Salt'")

    args = parser.parse_args()

    ymdh = args.ymdh
    rtofs_ic_netcdf_file = args.rtofs_ic_netcdf_file
    depth_var_name = args.depth_var_name
    temp_var_name = args.temp_var_name
    salt_var_name = args.salt_var_name
'''

################################################################################
#%% Reading grid
rtofs_ic = xr.open_dataset(rtofs_ic_netcdf_file,decode_times=False)
depth = np.asarray(rtofs_ic[depth_var_name][:])
temp = np.asarray(rtofs_ic[temp_var_name][0,:,:,:])
salt = np.asarray(rtofs_ic[salt_var_name][0,:,:,:])

rtofs_ic.close()

#################################################################################
# Adjust profiles
drho = 0.125
ref_depth = 10 # meters
adj_temp_prof = np.empty((len(depth),temp.shape[1],temp.shape[2]))
adj_temp_prof[:] = np.nan
adj_salt_prof = np.empty((len(depth),temp.shape[1],temp.shape[2]))
adj_salt_prof[:] = np.nan
for y in np.arange(temp.shape[1]):
    print(y)
    for x in np.arange(temp.shape[2]):
        temp_prof = temp[:,y,x]
        salt_prof = salt[:,y,x]
        dens_prof = dens(salt_prof,temp_prof,depth)
        mld, _ = MLD_dens_crit(drho,ref_depth,depth,temp_prof,dens_prof)
        if np.isfinite(mld):
            inside_mld = np.where(depth <= mld)[0]

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
'''
nc_file = nc.Dataset(rtofs_ic_netcdf_file,'a')
nc_file[temp_var_name][0,:,:,:] = adj_temp_prof
nc_file[salt_var_name][0,:,:,:] = adj_salt_prof
nc_file.close()
'''

