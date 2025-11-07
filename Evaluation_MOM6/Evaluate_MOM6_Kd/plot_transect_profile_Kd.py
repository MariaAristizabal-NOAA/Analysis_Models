#!/usr/bin/env python3

"""This script plots the profile of Eddy Diffusivity Kd. """

import yaml
import xarray as xr
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

#================================================================
# Parse the yaml config file
print('Parse the config file: plot_transect_profile_Kd.yml:')
with open('plot_transect_profile_Kd.yml', 'rt') as f:
    conf = yaml.safe_load(f)
ncfiles = conf['ncfiles']
labels = conf['labels']

#================================================================
nc = xr.open_dataset(ncfiles[0])
zi = np.asarray(nc.zi)

Kd_interf_prof = np.empty((len(ncfiles),len(zi)))
Kd_interf_prof[:] = np.nan
for i,ncfile in enumerate(ncfiles):
    nc = xr.open_dataset(ncfile)
    xh = np.asarray(nc.xh)
    yh = np.asarray(nc.yh)
    zi = np.asarray(nc.zi)
    kdd_interface = np.asarray(nc['Kd_interface'][0,:,:,:])

    #================================================================
    # Transect and profile
    xpos = int(len(xh)/2-1)
    ypos = int(len(yh)/2-1)
    xnorm = (xh - xh[xpos])/50
    ynorm = (yh - yh[ypos])/50
    Kd_trans = kdd_interface[:,ypos,:]
    Kd_interf_prof[i,:] = kdd_interface[:,ypos,xpos]

    fig,ax = plt.subplots(figsize=(8,4))
    cmap = plt.colormaps['YlOrRd']
    ncolors = cmap.N
    levels = np.arange(0,0.35,0.05)
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N,extend='max')
    ctr = ax.pcolor(xnorm,-zi,Kd_trans,cmap='YlOrRd',norm=norm)
    cbar = fig.colorbar(ctr,extendrect=True)
    cbar.set_label('$m/s^2$',fontsize=14)
    ax.set_ylim([-100,0])
    ax.set_xlim([-5,5])
    plt.ylabel('Depth (m)')
    plt.xlabel('Distance (RMW)')
    plt.title('Kd interface '+labels[i])

#================================================================
fig,ax = plt.subplots(figsize=(6,8))
for i in np.arange(len(ncfiles)):
    plt.plot(Kd_interf_prof[i,:],-zi,'.-',label=labels[i])
plt.xlabel('Kd')
plt.ylabel('Depth (m)')
plt.ylim([-100,0])
plt.title('Kd interface at ')
plt.legend()

