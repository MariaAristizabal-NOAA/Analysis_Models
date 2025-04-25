
#%% User input

# forecasting cycle to be used
cycle = '2018123112'

folder_mom6 = '/work/noaa/hwrf/noscrub/maristiz/MOM6_MLD_examp/'

folder_uom= '/work/noaa/hwrf/noscrub/maristiz/MOM6_MLD_examp/Upper_ocean_metrics/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import sys
import os
import os.path
import glob
import sys
import seawater as sw

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#################################################################################
#%% Get list files
files_mom6 = [folder_mom6 + 'ocn.ana.'+cycle+'.nc']
#files_mom6 = [folder_mom6 + 'ocn.bkg.'+cycle+'.nc']

################################################################################
#%% Reading MOM6 grid
mom6 = xr.open_dataset(files_mom6[0],decode_times=False)
lon_mom6 = np.asarray(mom6['lonh'][:])
lat_mom6 = np.asarray(mom6['lath'][:])
depth_mom6 = np.asarray(mom6['h'][:])[0,:,:,:]

#################################################################################
#%% Read MOM6 time
time_mom6 = []
for n,file in enumerate(files_mom6):
    print(file)
    mom6 = xr.open_dataset(file,decode_times=True)
    t = np.asarray(mom6.variables['Time'][:])[0]
    #timestamp = mdates.date2num(t)
    #timestamp = mdates.num2date(t)
    #time_mom6.append(mdates.num2date(timestamp))

#time_mom6 = np.asarray(time_mom6)

#################################################################################
#%% Get MLT around 2 degrees of storm eye HAFS/HYCOM
dtemp =  0.2
ref_depth = 2 #meters

MLT_mom6 = np.empty((len(lat_mom6),len(lon_mom6)))
MLT_mom6[:] = np.nan
MLD_mom6 = np.empty((len(lat_mom6),len(lon_mom6)))
MLD_mom6[:] = np.nan
OHC_mom6 = np.empty((len(lat_mom6),len(lon_mom6)))
OHC_mom6[:] = np.nan

file = files_mom6[0]
print(file)
mom6 = xr.open_dataset(file)
temp_mom6 = np.asarray(mom6['Temp'][0,:,:,:])
salt_mom6 = np.asarray(mom6['Salt'][0,:,:,:])

for x in np.arange(len(lon_mom6)):
    print(x)
    for y in np.arange(len(lat_mom6)):
        #print(y)
        if len(np.where(depth_mom6[:,y,x] > ref_depth)[0]) < 2:
            MLD_mom6[y,x] = np.nan
            MLT_mom6[y,x] = np.nan
        else:
           MLD_mom6[y,x] , MLT_mom6[y,x] = MLD_temp_crit(dtemp,ref_depth,depth_mom6[:,y,x],temp_mom6[:,y,x])

for x in np.arange(len(lon_mom6)):
    print(x)
    for y in np.arange(len(lat_mom6)):
        nok = depth_mom6[:,y,x] == 0.0
        depth_mom6[nok,y,x] = np.nan
        dens_prof = sw.dens(salt_mom6[:,y,x],temp_mom6[:,y,x],depth_mom6[:,y,x])
        OHC_mom6[y,x] = OHC_from_profile(depth_mom6[:,y,x],temp_mom6[:,y,x],dens_prof)

###########################################################################
#%% figures

kw = dict(levels=np.arange(0,281,10))
plt.figure()
plt.contourf(lon_mom6,lat_mom6,MLD_mom6,cmap='RdYlBu_r',**kw)
c = plt.colorbar()
c.set_label('m',fontsize=14)
plt.title('Mixed Layer Depth')
plt.axis('scaled')
plt.savefig('MLD_mom6')

kw = dict(levels=np.arange(-5,36,1))
plt.figure()
plt.contourf(lon_mom6,lat_mom6,MLT_mom6,cmap='Spectral_r',**kw)
c = plt.colorbar()
plt.contour(lon_mom6,lat_mom6,MLT_mom6,[26],colors='k')
c.set_label('$^oC$',fontsize=14)
plt.title('Mixed Layer Temperature')
plt.axis('scaled')
plt.savefig('MLT_mom6')

kw = dict(levels=np.arange(0,101,10))
plt.figure()
plt.contourf(lon_mom6,lat_mom6,OHC_mom6,cmap='gist_heat_r',**kw)
c = plt.colorbar()
c.set_label('$kJ/cm^2$',fontsize=14)
plt.title('Ocean Heat Content')
plt.axis('scaled')
plt.savefig('OHC_mom6')

