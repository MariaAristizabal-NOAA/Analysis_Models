
#%% User input

folder = '/scratch2/NOS/nosofs/Maria.Aristizabal/dorian05l.2019082800_exp/'

# POM grid file name
grid_file = folder + 'dorian05l.2019082800.pom.grid.nc'

# POM files
prefix = 'dorian05l.2019082800.pom.'

# Name of 3D variable
var_name = 't'

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

# for glider ng665 cycle 2019082800_exp
oklatpom_ng665 = np.tile(141,21)
oklonpom_ng665 = np.array([339, 339, 339, 339, 339, 339, 339, 339, 339, 339, 339, 339, 339,339, 339, 339, 339, 339, 339, 339, 339])

#%% Reading POM grid files
pom_grid = xr.open_dataset(grid_file)
lonc = np.asarray(pom_grid['east_e'][:])
latc = np.asarray( pom_grid['north_e'][:])
zlevc = np.asarray(pom_grid['zz'][:])
topoz = np.asarray(pom_grid['h'][:])

#%% Getting list of POM files
ncfiles = sorted(glob.glob(os.path.join(folder,prefix+'*0*.nc')))

# Reading POM time
time_pom = []
for i,file in enumerate(ncfiles):
    print(i)
    pom = xr.open_dataset(file)
    tpom = pom['time'][:]
    timestamp_pom = date2num(tpom)[0]
    time_pom.append(num2date(timestamp_pom))

time_pom = np.asarray(time_pom)
timestamp_pom = date2num(time_pom)

# Reading POM temperature
target_temp_pom = np.empty((len(zlevc),len(ncfiles),))
target_temp_pom[:] = np.nan
target_salt_pom = np.empty((len(zlevc),len(ncfiles),))
target_salt_pom[:] = np.nan
target_topoz_pom = np.empty((len(ncfiles),))
target_topoz_pom[:] = np.nan
target_dens_pom = np.empty((len(zlevc),len(ncfiles),))
target_dens_pom[:] = np.nan
for x,file in enumerate(ncfiles):
    print(x)
    pom = xr.open_dataset(file)

    # Interpolating latg and longlider into RTOFS grid
    #sublonpom = np.interp(timestamp_pom,timestampg,long)
    #sublatpom = np.interp(timestamp_pom,timestampg,latg)
    #oklonpom = np.int(np.round(np.interp(target_lon,lonc[0,:],np.arange(len(lonc[0,:])))))
    #oklatpom = np.int(np.round(np.interp(target_lat,latc[:,0],np.arange(len(latc[:,0])))))

    target_temp_pom[:,x] = np.asarray(pom['t'][0,:,oklatpom_ng665[x],oklonpom_ng665[x]])
    target_salt_pom[:,x] = np.asarray(pom['s'][0,:,oklatpom_ng665[x],oklonpom_ng665[x]])
    target_topoz_pom[x] = np.asarray(topoz[oklatpom_ng665[x],oklonpom_ng665[x]])
    target_dens_pom[:,x] = np.asarray(pom['rho'][0,:,oklatpom_ng665[x],oklonpom_ng665[x]])

z_matrix_pom = np.dot(target_topoz_pom.reshape(-1,1),zlevc.reshape(1,-1)).T
target_dens_pom = target_dens_pom * 1000 + 1000

target_temp_pom[target_temp_pom == 0.0] = np.nan
target_salt_pom[target_salt_pom == 0.0] = np.nan
target_dens_pom[target_dens_pom == 1000.0] = np.nan

#%%  Calculation of mixed layer depth based on dt, Tmean: mean temp within the 
# mixed layer and td: temp at 1 meter below the mixed layer
# for POM output            

dt = 0.2

MLD_dt_pom = np.empty(len(time_pom)) 
MLD_dt_pom[:] = np.nan
Tmean_dtemp_pom = np.empty(len(time_pom)) 
Tmean_dtemp_pom[:] = np.nan
Smean_dtemp_pom = np.empty(len(time_pom)) 
Smean_dtemp_pom[:] = np.nan
Td_pom = np.empty(len(time_pom)) 
Td_pom[:] = np.nan
for t,tt in enumerate(time_pom):
    d10 = np.where(z_matrix_pom[:,t] >= -10)[0][-1]
    T10 = target_temp_pom[d10,t]
    delta_T = T10 - target_temp_pom[:,t] 
    ok_mld = np.where(delta_T <= dt)[0]    
    if ok_mld.size == 0:
        MLD_dt_pom[t] = np.nan
        Tmean_dtemp_pom[t] = np.nan
        Smean_dtemp_pom[t] = np.nan
        Td_pom[t] = np.nan
    else:
        ok_mld_plus1m = np.where(z_matrix_pom >= z_matrix_pom[ok_mld[-1]] + 1)[0][0]
        MLD_dt_pom[t] = z_matrix_pom[ok_mld[-1],t]
        Tmean_dtemp_pom[t] = np.nanmean(target_temp_pom[ok_mld,t])
        Smean_dtemp_pom[t] = np.nanmean(target_salt_pom[ok_mld,t])
        Td_pom[t] = target_temp_pom[ok_mld_plus1m,t]

#%%  Calculation of mixed layer depth based on drho
# for POM output        
        
drho = 0.125

MLD_drho_pom = np.empty(len(time_pom)) 
MLD_drho_pom[:] = np.nan
Tmean_drho_pom = np.empty(len(time_pom)) 
Tmean_drho_pom[:] = np.nan
Smean_drho_pom = np.empty(len(time_pom)) 
Smean_drho_pom[:] = np.nan
for t,tt in enumerate(time_pom):
    d10 = np.where(z_matrix_pom[:,t] >= -10)[0][-1]
    rho10 = target_dens_pom[d10,t]
    delta_rho = -(rho10 - target_dens_pom[:,t]) 
    ok_mld = np.where(delta_rho <= drho)[0]
    if ok_mld.size == 0:
        MLD_drho_pom[t] = np.nan
        Tmean_drho_pom[t] = np.nan
        Smean_drho_pom[t] = np.nan
    else:
        MLD_drho_pom[t] = z_matrix_pom[ok_mld[-1],t] 
        Tmean_drho_pom[t] = np.nanmean(target_temp_pom[ok_mld,t]) 
        Smean_drho_pom[t] = np.nanmean(target_salt_pom[ok_mld,t]) 

#%% Temp profile during Dorian
tDorian = datetime(2019,8,28,18)
okpom = np.where(date2num(time_pom) == date2num(tDorian))[0][0]
temp_prof_pom = target_temp_pom[:,okpom]
depth_prof_pom = z_matrix_pom[:,okpom]

#%% Figure Temp time series profile
time_matrixpom = np.tile(date2num(time_pom),(z_matrix_pom.shape[0],1))
kw = dict(levels = np.linspace(19,31,13))

fig, ax = plt.subplots(figsize=(12, 2))
plt.ion()
#plt.contour(time_matrixpom,z_matrix_pom,target_temp_pom,colors = 'lightgrey',**kw)
plt.contour(time_matrixpom,z_matrix_pom,target_temp_pom,[26],colors = 'k')
plt.contourf(time_matrixpom,z_matrix_pom,target_temp_pom,cmap=cmocean.cm.thermal,**kw)
plt.plot(time_matrixpom[0,:],MLD_dt_pom,'^-',label='MLD dt',color='indianred',linewidth=2 )
plt.plot(time_matrixpom[0,:],MLD_drho_pom,'^-',label='MLD drho',color='seagreen',linewidth=2 )
ax.set_ylim(-200,0)
yl = ax.set_ylabel('Depth (m)',fontsize=14) #,labelpad=20)
cbar = plt.colorbar()
cbar.ax.set_ylabel('($^\circ$C)',fontsize=14)
ax.set_title('Along Track Temperature Profile HWRF-POM',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(-200,0))) #ng665
plt.plot(tDorian,np.arange(-200,0),'--k')
plt.legend()

file = "/home/Maria.Aristizabal/Dorian_2019/Figures/HWRF_pom_temp_Dorian_2019082818_exp.png"
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Figure salinity time series profile
time_matrixpom = np.tile(date2num(time_pom),(z_matrix_pom.shape[0],1))
#kw = dict(levels = np.linspace(35.5,37.3,19))
kw = dict(levels = np.linspace(35.5,37.3,19))

fig, ax = plt.subplots(figsize=(12, 2))
plt.ion()
#plt.contour(time_matrixpom,z_matrix_pom,target_salt_pom,colors = 'lightgrey',**kw)
plt.contourf(time_matrixpom,z_matrix_pom,target_salt_pom,cmap=cmocean.cm.haline ,**kw)
plt.plot(time_matrixpom[0,:],MLD_dt_pom,'^-',label='MLD dt',color='indianred',linewidth=2 )
plt.plot(time_matrixpom[0,:],MLD_drho_pom,'^-',label='MLD drho',color='seagreen',linewidth=2 )
ax.set_ylim(-200,0)
yl = ax.set_ylabel('Depth (m)',fontsize=14) #,labelpad=20)
cbar = plt.colorbar()
#cbar.ax.set_ylabel('($^\circ$C)',fontsize=16)
ax.set_title('Along Track Salinity Profile HWRF-POM',fontsize=14)
xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(-200,0))) #ng665
plt.plot(tDorian,np.arange(-200,0),'--k')
plt.legend()

file = "/home/Maria.Aristizabal/Dorian_2019/Figures/HWRF_pom_salt_Dorian_2019082818_exp.png"
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)


np.savez('POM_Dorian_2019082800_exp.npz',time_pom_exp=time_pom,timestamp_pom_exp=timestamp_pom,target_temp_pom_exp=target_temp_pom,target_salt_pom_exp=target_salt_pom,z_matrix_pom_exp=z_matrix_pom,MLD_dt_pom_exp=MLD_dt_pom,MLD_drho_pom_exp=MLD_drho_pom,Tmean_dtemp_pom_exp=Tmean_dtemp_pom,Tmean_drho_pom_exp=Tmean_drho_pom,temp_prof_pom_exp=temp_prof_pom,depth_prof_pom_exp=depth_prof_pom,target_dens_pom_exp=target_dens_pom,Smean_dtemp_pom_exp=Smean_dtemp_pom,Smean_drho_pom_exp=Smean_drho_pom)
         

