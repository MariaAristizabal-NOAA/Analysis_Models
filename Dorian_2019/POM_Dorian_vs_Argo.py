
#%% User input

folder_bef = '/scratch2/NOS/nosofs/Maria.Aristizabal/dorian05l.2019082512_oper/'
folder_aft = '/scratch2/NOS/nosofs/Maria.Aristizabal/dorian05l.2019090412_oper/'

# POM grid file name
grid_file_bef = folder_bef + 'dorian05l.2019082512.pom.grid.nc'
grid_file_aft = folder_aft + 'dorian05l.2019090412.pom.grid.nc'

# POM files
file_bef = folder_bef + 'dorian05l.2019082512.pom.0000.nc'
file_aft = folder_aft + 'dorian05l.2019090412.pom.0000.nc'

# Name of 3D variable
var_name = 't'

#Argo lat and lon
argo_lat_bef = 27.09778
argo_lon_bef = -76.60314
argo_lat_aft = 26.77398
argo_lon_aft = -76.64989

#%% 
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import netCDF4
from datetime import datetime,timedelta
from matplotlib.dates import date2num, num2date
import matplotlib.dates as mdates
import os
import os.path
import glob
import cmocean
from mpl_toolkits.basemap import Basemap

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Reading POM before
pom_grid = xr.open_dataset(grid_file_bef)
lonc = np.asarray(pom_grid['east_e'][:])
latc = np.asarray( pom_grid['north_e'][:])
zlevc = np.asarray(pom_grid['zz'][:])
topoz = np.asarray(pom_grid['h'][:])
pom = xr.open_dataset(file_bef)
tpom_bef = np.asarray(pom['time'][:])[0]
oklon_bef = np.int(np.round(np.interp(argo_lon_bef,lonc[0,:],np.arange(len(lonc[0,:])))))
oklat_bef = np.int(np.round(np.interp(argo_lat_bef,latc[:,0],np.arange(len(latc[:,0])))))
temp_pom_bef = np.asarray(pom['t'][0,:,oklat_bef,oklon_bef])
salt_pom_bef = np.asarray(pom['s'][0,:,oklat_bef,oklon_bef])
topoz_pom_bef = np.asarray(topoz[oklat_bef,oklon_bef])
z_matrix_pom_bef = np.dot(topoz_pom_bef.reshape(-1,1),zlevc.reshape(1,-1))
temp_pom_bef[temp_pom_bef == 0.0] = np.nan
salt_pom_bef[salt_pom_bef == 0.0] = np.nan

#%% Reading POM after
pom_grid = xr.open_dataset(grid_file_aft)
lonc = np.asarray(pom_grid['east_e'][:])
latc = np.asarray( pom_grid['north_e'][:])
zlevc = np.asarray(pom_grid['zz'][:])
topoz = np.asarray(pom_grid['h'][:])
pom = xr.open_dataset(file_aft)
tpom_aft = np.asarray(pom['time'][:])[0]
oklon_aft = np.int(np.round(np.interp(argo_lon_aft,lonc[0,:],np.arange(len(lonc[0,:])))))
oklat_aft = np.int(np.round(np.interp(argo_lat_aft,latc[:,0],np.arange(len(latc[:,0])))))
temp_pom_aft = np.asarray(pom['t'][0,:,oklat_aft,oklon_aft])
salt_pom_aft = np.asarray(pom['s'][0,:,oklat_aft,oklon_aft])
topoz_pom_aft = np.asarray(topoz[oklat_aft,oklon_aft])
z_matrix_pom_aft = np.dot(topoz_pom_aft.reshape(-1,1),zlevc.reshape(1,-1))
temp_pom_aft[temp_pom_aft == 0.0] = np.nan
salt_pom_aft[salt_pom_aft == 0.0] = np.nan

