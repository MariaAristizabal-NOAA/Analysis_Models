#%% User input

# forecasting cycle to be used
cycle = '2024011506'

rtofs_ic_file = '/work2/noaa/hwrf/save/maristiz/scripts_to_prep_MOM6/RTOFS_OBC_esmpy/2024011506_nhc_SST_2DVAR_test3/rtofs.f6_ocean_ts_obc.nc'

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
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Reading grid
rtofs_ic = xr.open_dataset(rtofs_ic_file,decode_times=False)
depth = np.asarray(rtofs_ic['Depth'][:])
temp = np.asarray(rtofs_ic['pot_temp'][0,:,:,:])

rtofs_ic.close()

#################################################################################
# Adjust profiles
max_depth = 150
dtemp = 0.2
ref_depth = 10 # meters
adj_temp_prof = np.empty((len(depth),temp.shape[1],temp.shape[2]))
adj_temp_prof[:] = np.nan
for y in np.arange(temp.shape[1]):
    print(y)
    for x in np.arange(temp.shape[2]):
        temp_prof = temp[:,y,x]
        mld, _ = MLD_temp_crit(dtemp,ref_depth,depth,temp_prof)
        if np.isfinite(mld):
            inside_mld = np.where(depth <= mld)[0]
            ok_mld = inside_mld[-1] 
            ok_max_depth = np.where(depth <= max_depth)[0][-1]

            adj_temp_pr = np.copy(temp_prof)
            adj_temp_pr[inside_mld] = temp_prof[0]

            for k in np.arange(1,ok_max_depth-ok_mld):
                adj_temp_pr[ok_mld+k] = (temp_prof[ok_mld+k] + adj_temp_pr[ok_mld+k-1])/2 

            adj_temp_prof[:,y,x] = adj_temp_pr
    
        else:
            adj_temp_prof[:,y,x] = temp_prof


y = 300
x = 2000
plt.figure()
plt.plot(temp[:,y,x],-depth,'.-')
plt.plot(adj_temp_prof[:,y,x],-depth,'.-')
plt.ylim([-250,0])
plt.xlim([10,30])

#################################################################################
# Write new adjusted temp profiles to file
nc_file = nc.Dataset(rtofs_ic_file,'a')
nc_file['pot_temp'][0,:,:,:] = adj_temp_prof
nc_file.close()

#%% Reading modified file
rtofs_ic_adj = xr.open_dataset(rtofs_ic_file,decode_times=False)
#rtofs_ic_adj = nc.Dataset(rtofs_ic_file)
temp_adj = np.asarray(rtofs_ic_adj['pot_temp'][0,:,:,:])

y = 300
x = 2000
plt.figure()
plt.plot(temp[:,y,x],-depth,'.-')
plt.plot(temp_adj[:,y,x],-depth,'.-')
plt.ylim([-250,0])
plt.xlim([10,30])

