#User input

ncfile = '/work/noaa/hwrf/scrub/maristiz/hafs_restart_RTOFS_increm_update/com/2020072812/00L/natl00l.2020072812.hafs_hycom_hat10_3z.f006.nc'

'''
var_name = 'temp'
klayer = '1'
#colormap='nipy_spectral'
#colormap='seismic'
colormap='Spectral_r'
min_val = 12
max_val = 33
delta_val = 1  # delta in colorbar
delta_contours = 2 # delta in contour plot
units = '$^oC$'
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]
'''

'''
var_name = 'srfhgt'
klayer = '0'
colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = -110
max_val = 110
delta_val = 10
delta_contour = 30
units = 'Cm'
lon_lim = [-98.23,-7.52]
lat_lim = [1.03,45.78]
'''

var_name = 'salinity'
klayer = 0
colormap='GnBu_r'
#colormap='nipy_spectral'
#colormap='seismic'
#colormap='Spectral_r'
min_val = 30
max_val = 40
delta_val = 0.5
delta_contour = 2
units = ' '
lon_lim = [-98.22,-7.52]
lat_lim = [1.03,45.78]

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import os.path
import glob

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

# Reading HYCOM grid
hycom = xr.open_dataset(ncfile,decode_times=False)
lon_hycom = np.asarray(hycom['Longitude'][:])
lat_hycom = np.asarray(hycom['Latitude'][:])
depth_hycom = np.asarray(hycom['Z'][:])

thycom = np.asarray(hycom['MT'][:])
timestamp = mdates.date2num(thycom)[0]
time_hycom = mdates.num2date(timestamp)
#time_hycom = np.asarray(time_hycom)
#timestamp_hycom = np.asarray(timestamp_hycom)

lon_limh = [lon_lim[0] + 360, lon_lim[1] + 360]

if np.min(lon_hycom) < 0:
    oklon = np.where(np.logical_and(lon_hycom>lon_lim[0],lon_hycom<lon_lim[1]))[0]
else:
    oklon = np.where(np.logical_and(lon_hycom>lon_limh[0],lon_hycom<lon_limh[1]))[0]
oklat = np.where(np.logical_and(lat_hycom>lat_lim[0],lat_hycom<lat_lim[1]))[0]

var_hycom = np.asarray(hycom[var_name][0,klayer,oklat,:][:,oklon])
maxval = np.nanmax(var_hycom)
minval = np.nanmin(var_hycom)
meanval = np.nanmean(var_hycom)
print(maxval)
print(minval)
print(meanval)

####################################################################
kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
fig,ax1 = plt.subplots(figsize = (11,6))
ax1.set_facecolor("bisque")
plt.contour(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,np.arange(min_val,max_val,delta_contour),colors='grey')
plt.contour(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,[0],colors='k')
plt.contourf(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,cmap=colormap,**kw)
plt.axis('scaled')
cbar = plt.colorbar(fraction=0.025, pad=0.04)
cbar.set_label(units,fontsize = 14)
plt.title(var_name + ' Depth ' + str(depth_hycom[klayer]) +  '  ' + ncfile.split('/')[-1].split('.')[-2],fontsize=18)
plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
plt.xlim(lon_lim)
plt.ylim(lat_lim)

if var_name == 'srfhgt':
    var_mean = var_hycom - meanval
    kw = dict(levels=np.arange(min_val,max_val+0.1,delta_val))
    fig,ax1 = plt.subplots(figsize = (11,6))
    ax1.set_facecolor("bisque")
    plt.contour(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,np.arange(min_val,max_val,delta_contour),colors='grey')
    plt.contour(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,[0],colors='k')
    plt.contourf(lon_hycom[oklon]-360,lat_hycom[oklat],var_hycom,cmap=colormap,**kw)
    plt.axis('scaled')
    cbar = plt.colorbar(fraction=0.025, pad=0.04)
    cbar.set_label(units,fontsize = 14)
    plt.title(var_name + ' Depth ' + str(depth_hycom[klayer]),fontsize=18)
    plt.text(260-360,-5,'max val = ' + str(maxval),fontsize=14)
    plt.text(260-360,-8,'min val = ' + str(minval),fontsize=14)
    plt.xlim(lon_lim)
    plt.ylim(lat_lim)

####################################################################

