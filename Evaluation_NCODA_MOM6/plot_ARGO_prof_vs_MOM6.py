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

#sys.path.append('/home/Maria.Aristizabal/RTOFS_utilities_Dmitry/MyPython')
#import mod_read_ncoda as rncoda

sys.path.append('/home/Maria.Aristizabal/RTOFS_utilities_Dmitry/anls_mom6')
import mod_mom6

###################################################################
# Parse the yaml config file
print('Parse the config file: plot_ARGO_prof_vs_MOM6.yml:')
with open('plot_ARGO_prof_vs_MOM6.yml', 'rt') as f:
    conf = yaml.safe_load(f)
rdate = conf['rdate']
ncoda_file = conf['ncoda_file']
mom6_backg = conf['mom6_backg']
mom6_anal = conf['mom6_anal']
mom6_f24 = conf['mom6_f24']

print(conf)

# get bacground date:
#rbgrd = rncoda.anls2bckground_date(rdate)

# Read MOM6 files
print('Reading ' + mom6_anal)
nc = xr.open_dataset(mom6_anal)
xh_anal = np.asarray(nc['xh'])
yh_anal = np.asarray(nc['yh'])
zl_anal = np.asarray(nc['zl'])
time_anal = nc['Time']
temp_anal = np.asarray(nc['potT'])[0,:,:,:]

_,zm_anal = mod_mom6.get_zz_zm(mom6_anal, f_btm=True)


print('Reading ' + mom6_f24)
nc = xr.open_dataset(mom6_f24)
xh_f24 = np.asarray(nc['xh'])
yh_f24 = np.asarray(nc['yh'])
zl_f24 = np.asarray(nc['zl'])
time_f24 = nc['Time']
temp_f24 = np.asarray(nc['potT'])[0,:,:,:]

_,zm_f24 = mod_mom6.get_zz_zm(mom6_f24, f_btm=True)

nc = xr.open_dataset(mom6_backg)
xh_backg = np.asarray(nc['xh'])
yh_backg = np.asarray(nc['yh'])
zl_backg = np.asarray(nc['zl'])
time_backg = nc['Time']
temp_backg = np.asarray(nc['potT'])[0,:,:,:]

_,zm_backg = mod_mom6.get_zz_zm(mom6_backg, f_btm=True)

print('Reading '+ncoda_file)
try:
    fga = open(ncoda_file,'rb')
except:
    print('Could not open '+ncoda_file)

fga.seek(0)
n_data  = np.fromfile(fga, dtype='>i4', count=2)[1]
time  = np.fromfile(fga, dtype='>f', count=1)[0]
lat  = np.fromfile(fga, dtype='>f', count=n_data+2)[2:]
lon = np.fromfile(fga, dtype='>f', count=n_data+2)[2:]
lvl = np.fromfile(fga, dtype='>f', count=n_data+2)[2:]
ndx = np.fromfile(fga, dtype='>i4', count=n_data+2)[2:]
tmp_typ = np.fromfile(fga, dtype='>i4', count=n_data+2)[2:]
tmp = np.fromfile(fga, dtype='>f', count=n_data+2)[2:]
#sgn = np.fromfile(fga, dtype='>str', count=n_data+2)[2:]

xhg = np.empty((len(xh_backg)))
xhg[:] = np.nan
# Transform MOM6 longitude to geographic longitude
for xpos,x in enumerate(xh_backg):
    if x < -180:
        xhg[xpos] = x + 360
    else:
        xhg[xpos] = x

oksort = np.argsort(xhg)
xhg_sort = xhg[oksort]
temp_backg_sort = temp_backg[:,:,oksort]
temp_anal_sort = temp_anal[:,:,oksort]
temp_f24_sort = temp_f24[:,:,oksort]
zm_backg_sort = zm_backg[:,:,oksort]
zm_anal_sort = zm_anal[:,:,oksort]
zm_f24_sort = zm_f24[:,:,oksort]

# Transform NCODA longitude to geographic longitude
long = np.empty((len(lon)))
long[:] = np.nan
for xpos,lo in enumerate(lon):
    if lo > 180:
        long[xpos] = lo - 360 
    else:
        long[xpos] = lo

'''
plt.figure()
plt.contourf(xh_backg,yh_backg,temp_f24[0,:,:],cmap='Spectral_r')
plt.plot(long,lat,'.k',markersize=1)
plt.colorbar()
plt.title(mom6_f24.split('/')[-1])
'''

plt.figure(figsize=(8,5))
plt.contourf(xhg_sort,yh_backg,temp_f24_sort[0,:,:],cmap='Spectral_r')
plt.plot(long,lat,'.k',markersize=1)
plt.colorbar(shrink=0.6)
plt.axis('scaled')
plt.title('NCODA: '+ncoda_file.split('/')[7]+'/'+ncoda_file.split('/')[-1]+'\n'+
          'Number Argo profiles = '+str(np.max(ndx))+'\n'+
          'Total number of data points = '+str(n_data))

temp_backg_interp_NCODA_lvl = []
temp_f24_interp_NCODA_lvl = []
nprof = np.unique(ndx)
for n in nprof[0:3]:
    print('Reading profile '+str(n))
    lat_prof = lat[ndx==n][0] 
    lon_prof = long[ndx==n][0] 
    lvl_prof = lvl[ndx==n] 
    tmp_prof = tmp[ndx==n]
    
    okx = int(np.round(np.interp(lon_prof,xhg_sort,np.arange(len(xh_backg)))))
    oky = int(np.round(np.interp(lat_prof,yh_backg,np.arange(len(yh_backg)))))
    temp_backg_prof = temp_backg_sort[:,oky,okx]
    zm_backg_prof = zm_backg_sort[:,oky,okx]
    temp_f24_prof = temp_f24_sort[:,oky,okx]
    zm_f24_prof = zm_f24_sort[:,oky,okx]
    temp_anal_prof = temp_anal_sort[:,oky,okx]
    zm_anal_prof = zm_anal_sort[:,oky,okx]

    tmp_interp_zm_backg = np.interp(-zm_backg_prof,lvl_prof,tmp_prof)
    tmp_interp_zm_anal = np.interp(-zm_anal_prof,lvl_prof,tmp_prof)
    tmp_interp_zm_f24 = np.interp(-zm_f24_prof,lvl_prof,tmp_prof)

    plt.figure()
    plt.plot(tmp_prof,-lvl_prof,'.-r')
    plt.plot(tmp_interp_zm_backg,zm_backg_prof,'.-g')

##############################################################
    tback = np.asarray(time_backg)[0].strftime()[0:13]
    tanal = np.asarray(time_anal)[0].strftime()[0:13]
    tf24 = np.asarray(time_f24)[0].strftime()[0:13]

    plt.figure(figsize=(8,5))
    plt.contourf(xhg_sort,yh_backg,temp_f24_sort[0,:,:],cmap='Spectral_r')
    plt.plot(lon_prof,lat_prof,'.k',markersize=5)
    plt.colorbar(shrink=0.6)
    plt.axis('scaled')
    plt.title('Argo Float Number '+str(n)+
              ' Lat = '+str(np.round(lat_prof,2))+', Lon = '+str(np.round(lon_prof,2)))

    plt.subplots(figsize=(7,8))
    plt.plot(tmp_prof,-lvl_prof,'.-b',label='NCODA')
    plt.plot(temp_backg_prof,zm_backg_prof,'.-',color='orange',label='Background '+tback)
    plt.plot(temp_anal_prof,zm_anal_prof,'.-',color='green',label='Analysis '+tanal)
    plt.plot(temp_f24_prof,zm_f24_prof,'.-',color='red',label='Forecast '+tf24)
    plt.title('Argo Float Temperature'+'\n'+
              'NCODA: '+ncoda_file.split('/')[7]+'/'+ncoda_file.split('/')[-1]+'\n'+
              'Background: '+mom6_backg.split('MOM6')[-1]+'\n'+
              'Analysis: '+mom6_anal.split('MOM6')[-1]+'\n'+
              'Forecast: '+mom6_f24.split('MOM6')[-1]+'\n'+
              'Lat = '+str(np.round(lat_prof,2))+', Lon = '+str(np.round(lon_prof,2)))
    plt.legend()
    plt.ylim(-1000)

    plt.figure(figsize=(7,8))
    plt.plot(tmp_interp_zm_backg-temp_backg_prof,zm_backg_prof,'.-',color='orange',label='Background')
    plt.plot(tmp_interp_zm_anal-temp_anal_prof,zm_anal_prof,'.-',color='green',label='Analysis')
    plt.plot(tmp_interp_zm_f24-temp_f24_prof,zm_f24_prof,'.-',color='red',label='Forecast')
    plt.plot(np.zeros(1000),np.arange(-1000,0),'.-k')
    plt.xlim([-4,4])
    plt.ylim([-1000,0])
    plt.title('Observation Anomaly')
    plt.xlabel('Anomaly $(^oC)$')
    plt.legend()


##############################################################
    '''
    oklat = np.where(np.logical_and(lat>26.40,lat<26.41))
    lat_prof = lat[oklat][0]
    lon_prof = long[oklat][0]
    lvl_prof = lvl[oklat]
    tmp_prof = tmp[oklat]
    ndx_prof = ndx[oklat]
    okx = int(np.round(np.interp(lon_prof,xhg_sort,np.arange(len(xh_backg)))))
    oky = int(np.round(np.interp(lat_prof,yh_backg,np.arange(len(yh_backg)))))
    temp_backg_prof = temp_backg_sort[:,oky,okx]
    zm_backg_prof = zm_backg_sort[:,oky,okx]
    temp_f24_prof = temp_f24_sort[:,oky,okx]
    zm_f24_prof = zm_f24_sort[:,oky,okx]
    '''






