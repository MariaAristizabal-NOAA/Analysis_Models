
#%% User input
#cycle = '2019083112'
#cycle = '2019083100'
cycle = '2019082800'

Dir_HWRF_POM_oper = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.' + cycle + '_grb2_to_nc_oper/' 
Dir_HWRF_POM_exp = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.' + cycle +'_grb2_to_nc_exp/'
Dir_HWRF_HYCOM_exp = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.' + cycle +'_grb2_to_nc_exp/'
dir_figs = '/home/Maria.Aristizabal/Dorian_2019/Figures/'

#%%
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import netCDF4
from mpl_toolkits.basemap import Basemap
import cmocean
import glob
import os
from datetime import datetime,timedelta
import matplotlib.dates as mdates

#%% Get list HWRF files
HWRF_POM_oper = sorted(glob.glob(os.path.join(Dir_HWRF_POM_oper,'*.nc')))
HWRF_POM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_POM_exp,'*.nc')))
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(Dir_HWRF_HYCOM_exp,'*.nc')))

max_wind_10m_hwrf19_pom_oper = []
max_wind_10m_hwrf20_pom_exp = []
max_wind_10m_hwrf20_hycom_exp = []
time_hwrf = []
for fl in HWRF_POM_oper:
	HWRF = xr.open_dataset(fl)
	t_hwrf = np.asarray(HWRF.variables['time'][:])
	UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
	VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
	max_wind_10m_hwrf19_pom_oper.append(np.max(np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)))
	time_hwrf.append(t_hwrf)
time_hwrf = np.asarray(time_hwrf)

for fl in HWRF_POM_exp:
	HWRF = xr.open_dataset(fl)
	UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
	VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
	max_wind_10m_hwrf20_pom_exp.append(np.max(np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)))

for fl in HWRF_HYCOM_exp:
	HWRF = xr.open_dataset(fl)
	t_hwrf = np.asarray(HWRF.variables['time'][:])
	UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
	VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
	max_wind_10m_hwrf20_hycom_exp.append(np.max(np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)))
	#time_hwrf.append(t_hwrf)

#%% wind speed in knots
max_wind_10m_hwrf19_pom_oper = 1.94384 * np.asarray(max_wind_10m_hwrf19_pom_oper)
max_wind_10m_hwrf20_pom_exp = 1.94384 * np.asarray(max_wind_10m_hwrf20_pom_exp)
max_wind_10m_hwrf20_hycom_exp = 1.94384 * np.asarray(max_wind_10m_hwrf20_hycom_exp)

#%% Hurricane Dorian best track information
t0 = datetime(2019,8,24,12)
deltat= timedelta(hours=6) # every 6 hours
time_best_track = [t0+nstep*deltat for nstep in np.arange(63)] 
time_best_track = np.asarray(time_best_track)
wind_int_kt = np.array([ 30.,  35.,  35.,  40.,  40.,  45.,  45.,  50.,  50.,  45.,  45.,
        45.,  45.,  45.,  50.,  55.,  60.,  65.,  70.,  75.,  75.,  75.,
        80.,  90.,  95., 100., 115., 120., 125., 130., 130., 130., 150.,
       160., 155., 145., 135., 125., 120., 105., 100.,  95.,  95.,  95.,
        90.,  90., 100., 100., 100.,  95.,  85.,  80.,  80.,  80.,  80.,
        75.,  75.,  85.,  80.,  75.,  70.,  60.,  50.])

#%% Figure intensity
fig,ax1 = plt.subplots(figsize=(6, 4))
plt.ion()
plt.plot(time_best_track,wind_int_kt,'o-k',label='Best')
plt.plot(time_hwrf,max_wind_10m_hwrf19_pom_oper,'X-',color='mediumorchid',label='POM Oper',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_wind_10m_hwrf20_pom_exp,'^-',color='teal',label='POM Exp',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_wind_10m_hwrf20_hycom_exp,'H-',color='darkorange',label='HYCOM Exp',markeredgecolor='k',markersize=7)
plt.legend(loc='lower right')
#plt.legend()
plt.xlim([time_hwrf[0],time_hwrf[-1]])
xfmt = mdates.DateFormatter('%d \n %b')
ax1.xaxis.set_major_formatter(xfmt)
plt.ylim([20,165])
plt.title('Intensity Forecast Dorian '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([20,165])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

file_name = dir_figs + 'Dorian_intensity_cycle_'+cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 

#%% Intensity error

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

okt = np.logical_and(mdates.date2num(time_best_track) >= mdates.date2num(time_hwrf[0]),mdates.date2num(time_best_track) <= mdates.date2num(time_hwrf[-1]))
lead_time = np.arange(0,132,6)

int_err_hwrf19_pom_oper = (wind_int_kt[okt] - max_wind_10m_hwrf19_pom_oper[::2]) #*100/wind_int_kt[okt]
int_err_hwrf20_pom_exp = (wind_int_kt[okt] - max_wind_10m_hwrf20_pom_exp[::2]) #*100/wind_int_kt[okt]
int_err_hwrf20_hycom_exp = (wind_int_kt[okt] - max_wind_10m_hwrf20_hycom_exp[::2]) #*100/wind_int_kt[okt]

fig,ax1 = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time,int_err_hwrf19_pom_oper,'X-',color='mediumorchid',label='HWRF2019-POM Oper',markeredgecolor='k',markersize=7)
plt.plot(lead_time,int_err_hwrf20_pom_exp,'^-',color='teal',label='HRWF2020-POM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time,int_err_hwrf20_hycom_exp,'H-',color='darkorange',label='HWRF2020-HYCOM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time,np.tile(0,len(lead_time)),'--k')
plt.xlim([0,126])
ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
plt.title('Intensity Forecast Error Dorian '+ cycle,fontsize=18)
plt.ylabel('Forecast Error (Kt)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14)
plt.legend()

'''
#%% Read HWRF nc files
shtfl_maxwind = []
lhtfl_maxwind = []
for fl in HWRF_POM_oper:
	print(fl)
	HWRF = xr.open_dataset(fl)
	lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
	lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
	t_hwrf = np.asarray(HWRF.variables['time'][:])
	UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
	VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
	SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
	LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
	DSWRD_hwrf = np.asarray(HWRF.variables['DSWRF_surface'][0,:,:])
	USWRD_hwrf = np.asarray(HWRF.variables['USWRF_surface'][0,:,:])
	DLWRD_hwrf = np.asarray(HWRF.variables['DLWRF_surface'][0,:,:])
	ULWRD_hwrf = np.asarray(HWRF.variables['ULWRF_surface'][0,:,:])
	WTMP_hwrf = np.asarray(HWRF.variables['WTMP_surface'][0,:,:])
	SWRD_hwrf = DSWRD_hwrf - USWRD_hwrf
	LWRD_hwrf = DLWRD_hwrf - ULWRD_hwrf
	 
	wind_10m = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
	ok_maxwind = np.where(wind_10m == np.max(wind_10m))
	shtfl_maxwind.append(SHTFL_hwrf[ok_maxwind[0][0],ok_maxwind[1][0]])
	lhtfl_maxwind.append(LHTFL_hwrf[ok_maxwind[0][0],ok_maxwind[1][0]])
'''
'''
	#%% map wind vectors
	m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
	x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))

	plt.figure()
	plt.ion()
	m.drawcoastlines()
	m.fillcontinents()
	m.drawmapboundary()
	# draw parallels and meridians.
	#m.drawparallels(np.arange(-80,-60,5),labels=[True,True,True,True],dashes=[2,2])
	#m.drawmeridians(np.arange(15,35,5),labels=[False,False,False,True],dashes=[2,2])
	kw = dict(levels=np.linspace(0,70,8))
	plt.contourf(x,y,np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2),cmap=plt.cm.Spectral_r,**kw)
	c = plt.colorbar()
	q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
	xq,yq = m(-78,12.5)
	plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
	xc, yc = m(-77.4,27.0)
	plt.plot(xc,yc,'*r',markersize=7)
	plt.title('HWRF Wind speed on '+str(time_hwrf)[2:18],fontsize=16)
	c.set_label('m/s',rotation=90, labelpad=15, fontsize=16)
	c.ax.tick_params(labelsize=14)

	file_name = dir_figs + 'Dorian_wind_' + str(time_hwrf)[2:18]
	plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 
'''
'''
	#%% map sensible heat flux
	m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
	x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))
	plt.figure()
	plt.ion()
	m.drawcoastlines()
	m.fillcontinents()
	m.drawmapboundary()
	kw = dict(levels=np.linspace(-160,400,15))
	plt.contourf(x,y,SHTFL_hwrf,cmap=plt.cm.coolwarm,**kw)
	c = plt.colorbar()
	q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
	xq,yq = m(-78,12.5)
	plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
	xc, yc = m(-77.4,27.0)
	plt.plot(xc,yc,'*r',markersize=7)
	plt.title('HWRF Sensible Heat Flux on '+str(time_hwrf)[2:18],fontsize=16)
	c.set_label('$W/m^2$',rotation=90, labelpad=15, fontsize=16)
	c.ax.tick_params(labelsize=14) 
'''

'''
	#%% map latent heat flux
	m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
	x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))
	plt.figure()
	plt.ion()
	m.drawcoastlines()
	m.fillcontinents()
	m.drawmapboundary()
	kw = dict(levels=np.linspace(-200,1400,9))
	plt.contourf(x,y,LHTFL_hwrf,cmap=plt.cm.coolwarm,**kw)
	c = plt.colorbar()
	q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
	xq,yq = m(-78,12.5)
	plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
	xc, yc = m(-77.4,27.0)
	plt.plot(xc,yc,'*r',markersize=7)
	plt.title('HWRF Latent Heat Flux on '+str(time_hwrf)[2:18],fontsize=16)
	c.set_label('$W/m^2$',rotation=90, labelpad=15, fontsize=16)
	c.ax.tick_params(labelsize=14)
'''
'''
	#%% map Downward short wave radiation
	m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
	x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))
	plt.figure()
	plt.ion()
	m.drawcoastlines()
	m.fillcontinents()
	m.drawmapboundary()
	kw = dict(levels=np.linspace(0,800,9))
	plt.contourf(x,y,SWRD_hwrf,cmap=plt.cm.coolwarm,**kw)
	#plt.pause(0.1)
	c = plt.colorbar()
	q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
	xq,yq = m(-78,12.5)
	plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
	xc, yc = m(-77.4,27.0)
	plt.plot(xc,yc,'*r',markersize=7)
	plt.title('HWRF Shortwave Radiation on '+str(time_hwrf)[2:18],fontsize=16)
	c.set_label('$W/m^2$',rotation=90, labelpad=15, fontsize=16)
	c.ax.tick_params(labelsize=14)
'''

'''
	#%% map long wave radiation
	m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
	x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))
	plt.figure()
	plt.ion()
	m.drawcoastlines()
	m.fillcontinents()
	m.drawmapboundary()
	kw = dict(levels=np.linspace(0,800,9))
	plt.contourf(x,y,LWRD_hwrf,cmap=plt.cm.coolwarm)
	#plt.pause(0.1)
	c = plt.colorbar()	
	q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
	xq,yq = m(-78,12.5)
	plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})	
	xc, yc = m(-77.4,27.0)
	plt.plot(xc,yc,'*r',markersize=7)
	plt.title('HWRF Longwave Radiation on '+str(time_hwrf)[2:18],fontsize=16)
	c.set_label('$W/m^2$',rotation=90, labelpad=15, fontsize=16)
	c.ax.tick_params(labelsize=14)

'''
'''
#%% map SST 
m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=35,llcrnrlon=-80,urcrnrlon=-60,resolution='l')
x, y = m(*np.meshgrid(lon_hwrf,lat_hwrf))
plt.figure()
plt.ion()
m.drawcoastlines()
m.fillcontinents()
m.drawmapboundary()
kw = dict(levels=np.linspace(24,33,10))
plt.contourf(x,y,WTMP_hwrf-273.15,cmap=cmocean.cm.thermal,**kw)
c = plt.colorbar()
q = plt.quiver(x[::30,::30], y[::30,::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30])
xq,yq = m(-78,12.5)
plt.quiverkey(q,xq,yq,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
xc, yc = m(-77.4,27.0)
plt.plot(xc,yc,'*r',markersize=7)
plt.title('HWRF SST on '+str(time_hwrf)[2:18],fontsize=16)
c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=16)
c.ax.tick_params(labelsize=14)

'''
'''
shtfl_maxwind = np.asarray(shtfl_maxwind)
lhtfl_maxwind = np.asarray(lhtfl_maxwind)

fig,ax = plt.subplots()
plt.plot(time_hwrf,shtfl_maxwind,'.-',label='Sensible Heat Flux')
plt.plot(time_hwrf,lhtfl_maxwind,'.-',label='Latent Heat Flux')
plt.plot(time_hwrf,lhtfl_maxwind+shtfl_maxwind,'.-',label='Enthalpy')
plt.legend()
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
'''

