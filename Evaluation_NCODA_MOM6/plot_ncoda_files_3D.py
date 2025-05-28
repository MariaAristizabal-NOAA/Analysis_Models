"""
  Plot NCODA increment input fields
"""

import os
import sys
import yaml
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
#import cartopy
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature

sys.path.append('/home/Maria.Aristizabal/Analysis_Models/Evaluation_NCODA_MOM6/MyPython_Dmitry/ncoda_utils')
import mod_read_ncoda as rncoda

###################################################################
# Parse the yaml config file
print('Parse the config file: plot_ncoda.yml:')
with open('plot_ncoda_files_3D.yml', 'rt') as f:
    conf = yaml.safe_load(f)
rdate = conf['rdate']
mom6_file = conf['mom6_file']
ncoda_file = conf['ncoda_file']
rLayer = conf['rLayer']
min_value = conf['min_value']
max_value = float(conf['max_value'])
dvalue = float(conf['dvalue'])
cmap_name = conf['cmap_name']
#cartopy_dir = conf['cartopyDataDir']

print(conf)

# get bacground date:
rbgrd = rncoda.anls2bckground_date(rdate)

# Read MOM6 file
nc = xr.open_dataset(mom6_file)
xh = np.asarray(nc['xh'])
yh = np.asarray(nc['yh'])
zl = np.asarray(nc['zl'])

IDM = len(xh)
JDM = len(yh)
KDM = len(zl)

print('Reading '+ncoda_file)
field = rncoda.read_ncoda_inc(ncoda_file,IDM,JDM,KDM,rLayer=rLayer)

field[field == -999] = np.nan

##############################################################
fig,ax = plt.subplots(figsize=(8,4))

cflevels = np.arange(min_value,max_value,dvalue)
cmap = plt.get_cmap(cmap_name)
cf = ax.contourf(xh,yh,field,levels=cflevels,cmap=cmap,extend='both')
lb = ax.contour(xh,yh,field,levels=[26],colors='grey',alpha=0.7,linewidths=0.5)
cb = plt.colorbar(cf,orientation='vertical',pad=0.02,aspect=20,shrink=0.8,extendrect=True,ticks=cflevels[::4])
cb.ax.tick_params(labelsize=10)
ax.set_title(ncoda_file.split('/')[-1]+' Level = '+str(rLayer),fontsize=12)
ax.axis('scaled')

##############################################################
'''
central_longitude = 0

fig = plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
ax.axis('scaled')

cflevels = np.arange(-5,36,5)
cmap = plt.get_cmap('Spectral_r')
cf = ax.contourf(xh,yh,field,levels=cflevels,cmap=cmap,extend='both',transform=ccrs.PlateCarree())
lb = ax.contour(xh,yh,field,levels=[26],colors='grey',alpha=0.7,transform=ccrs.PlateCarree(),linewidths=0.5)
ax.clabel(lb, lb.levels,inline=True,fmt='%1.0f',fontsize=6,colors='grey')
cb = plt.colorbar(cf,orientation='vertical',pad=0.02,aspect=20,shrink=0.8,extendrect=True,ticks=cflevels[::4])
cb.ax.tick_params(labelsize=8)

# Add borders and coastlines
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3, facecolor='none', edgecolor='0.1')

ax.set_title(ncoda_file.split('/')[-1],fontsize=12)
'''




