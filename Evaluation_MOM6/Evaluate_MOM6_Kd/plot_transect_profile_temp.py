#!/usr/bin/env python3

"""This script plots the profile of Eddy Diffusivity Kd. """

import yaml
import xarray as xr
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_transect_profile_temp.yml:')
with open('plot_transect_profile_temp.yml', 'rt') as f:
    conf = yaml.safe_load(f)
ncfiles = conf['ncfiles']
labels = conf['labels']

#================================================================
nc = xr.open_dataset(ncfiles[0])
zl = np.asarray(nc.zl)

temp_prof = np.empty((len(ncfiles),len(zl)))
temp_prof[:] = np.nan
for i,ncfile in enumerate(ncfiles):
    nc = xr.open_dataset(ncfile)
    xh = np.asarray(nc.xh)
    yh = np.asarray(nc.yh)
    zl = np.asarray(nc.zl)
    temp = np.asarray(nc['temp'][0,:,:,:])

    #================================================================
    # Transect and profile
    xpos = int(len(xh)/2-1)
    ypos = int(len(yh)/2-1)
    temp_trans = temp[:,:,xpos]
    temp_prof[i,:] = temp[:,ypos,xpos]

    fig,ax = plt.subplots(figsize=(8,4))
    cmap = plt.colormaps['Spectral']
    ncolors = cmap.N
    levels = np.arange(20,31,1)
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N,extend='both')
    ctr = ax.pcolor(yh,-zl,temp_trans,cmap='Spectral',norm=norm)
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$m/s^2$',fontsize=14)
    ax.set_ylim([-200,0])
    plt.ylabel('Depth (m)')
    plt.xlabel('Distance (Km)')
    plt.title('Temperature '+labels[i])

#================================================================
fig,ax = plt.subplots(figsize=(6,8))
for i in np.arange(len(ncfiles)):
    plt.plot(temp_prof[i,:],-zl,'.-',label=labels[i])
plt.xlabel('Temperature ($^oC$)')
plt.ylabel('Depth')
plt.ylim([-200,0])
plt.xlim([20,32])
plt.legend()

