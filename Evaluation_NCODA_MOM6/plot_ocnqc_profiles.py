"""
  Plot ocnqc profiles
"""

import os
import sys
import yaml
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

###################################################################
# Parse the yaml config file
print('Parse the config file: plot_ocnqc_profiles.yml')
with open('plot_ocnqc_profiles.yml', 'rt') as f:
    conf = yaml.safe_load(f)
ocnqc_file = conf['ocnqc_file']
cartopyDataDir = conf['cartopyDataDir']

print(conf)

print('Reading '+ocnqc_file)
try:
    fga = open(ocnqc_file,'rb')
except:
    print('Could not open '+ocnqc_file)

fga.seek(0)
n_out  = np.fromfile(fga, dtype='>i4', count=2)[1]
max_lev  = np.fromfile(fga, dtype='>i4', count=1)[0]
vrsn  = np.fromfile(fga, dtype='>i4', count=1)[0]
btm  = np.fromfile(fga, dtype='>f', count=n_out+2)[2:]
lat  = np.fromfile(fga, dtype='>f', count=n_out+2)[2:]
lon  = np.fromfile(fga, dtype='>f', count=n_out+2)[2:]
ls  = np.fromfile(fga, dtype='>i4', count=n_out+2)[2:]
lt  = np.fromfile(fga, dtype='>i4', count=n_out+2)[2:]
sal_typ  = np.fromfile(fga, dtype='>i4', count=n_out+2)[2:]
sqc  = np.fromfile(fga, dtype='>f', count=n_out+2)[2:]
tmp_typ  = np.fromfile(fga, dtype='>i4', count=n_out+2)[2:]
tqc  = np.fromfile(fga, dtype='>f', count=n_out+2)[2:]

# Transform ocnqc longitude to geographic longitude
long = np.empty((len(lon)))
long[:] = np.nan
for xpos,lo in enumerate(lon):
    if lo > 180:
        long[xpos] = lo - 360
    else:
        long[xpos] = lo

lvl = np.empty((max_lev,n_out))
lvl[:] = np.nan
sal = np.empty((max_lev,n_out))
sal[:] = np.nan
sal_err = np.empty((max_lev,n_out))
sal_err[:] = np.nan
spbr = np.empty((max_lev,n_out))
spbr[:] = np.nan
tmp = np.empty((max_lev,n_out))
tmp[:] = np.nan
tmp_err = np.empty((max_lev,n_out))
tmp_err[:] = np.nan
tpbr = np.empty((max_lev,n_out))
tpbr[:] = np.nan
clm_sal = np.empty((max_lev,n_out))
clm_sal[:] = np.nan
cssd = np.empty((max_lev,n_out))
cssd[:] = np.nan
clm_tmp = np.empty((max_lev,n_out))
clm_tmp[:] = np.nan
ctsd = np.empty((max_lev,n_out))
ctsd[:] = np.nan
flg= np.empty((max_lev,n_out))
flg[:] = np.nan
for i in np.arange(n_out):
    print(i)
    lvl[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    sal[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    sal_err[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    spbr[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    tmp[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    tmp_err[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    tpbr[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    clm_sal[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    cssd[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    clm_tmp[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    ctsd[0:lt[i],i]  = np.fromfile(fga, dtype='>f', count=lt[i]+2)[2:]
    flg[0:lt[i],i]  = np.fromfile(fga, dtype='>i4', count=lt[i]+2)[2:]

dtg = np.fromfile(fga, dtype='str', count=n_out+2)[2:] 
rct = np.fromfile(fga, dtype='str', count=n_out+2)[2:] 
sgn = np.fromfile(fga, dtype='str', count=n_out+2)[2:] 

i=1000
myproj = ccrs.PlateCarree(0)
transform = ccrs.PlateCarree(0)
fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection=myproj)
ax.axis('equal')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
ax.plot(long,lat,'.k',markersize=1)
ax.plot(long[i],lat[i],'*r',markersize=5)
plt.title(ocnqc_file.split('/')[-1])

plt.figure(figsize=(4,7))
plt.plot(tmp[:,i],-lvl[:,i],'.')
plt.ylabel('Depth (m)')
plt.xlabel('Temperature ($^oC$)')
plt.title('Lat = '+str(lat[i])+' Lon ='+str(lon[i]))

plt.figure(figsize=(4,7))
plt.plot(sal[:,i],-lvl[:,i],'.')
plt.ylabel('Depth (m)')
plt.xlabel('Salinity ')
plt.title('Lat = '+str(lat[i])+' Lon ='+str(lon[i]))

myproj = ccrs.PlateCarree(0)
transform = ccrs.PlateCarree(0)
# create figure and axes instances
fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection=myproj)
ax.axis('equal')
ax.plot(long,lat,'.k',markersize=1)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')
plt.title(ocnqc_file.split('/')[-1]+'\n'+
          'Number of profiles = '+str(n_out))

##############################################################






