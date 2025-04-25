
#%% User input

#storm_id = '06l'
#cycle ='2021063018'
cycle ='2021090206'

Dir_input = '/scratch2/AOML/aoml-hafs1/Lew.Gramer/ocean/HAFS_A/uncoupled/grb2/' + cycle + '/'
Dir_output = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2a_uncoupled/' + cycle + '/'
#Dir_input = '/scratch2/AOML/aoml-hafs1/Lew.Gramer/ocean/HAFS_A/coupled/grb2/' + cycle + '/'
#Dir_output = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2a_coupled/' + cycle + '/'

import os 
import glob

#%% Load neccesary modules to use grib2 
#os.system('module load intel/19.0.5.281 wgrib2/2.0.8')
#%% Find files
#grb2_files = sorted(glob.glob(os.path.join(Dir_input,'*hwrfprs*storm*grb2')))
grb2_files = sorted(glob.glob(os.path.join(Dir_input,'*hafsprs*synoptic*grb2')))
#grb2_files = sorted(glob.glob(os.path.join(Dir_input,'*hafs*grb2')))


#%% convert from  grib2 to nc files
#for fl in grb2_files[7:]:
#for fl in [grb2_files[25],grb2_files[31]]:
for fl in grb2_files:
	print(fl)
	#ll = 'wgrib2 -v ' + fl + ' -match "TMP" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'TMP.nc'
	#os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "VVEL" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'VVEL.nc'
	#os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "DZDT" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'DZDT.nc'
	#os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "TMP" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'TMP.nc'
	#os.system(ll) 

	ll = 'wgrib2 -v ' + fl + ' -match "RH" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'RH.nc'
	os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "UGRD" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'UGRD.nc'
	#os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "VGRD" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'VGRD.nc'
	#os.system(ll) 

	#ll = 'wgrib2 -v ' + fl + ' -match "TMP|RH|UGRD|VGRD" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'nc'
	#ll = 'wgrib2 -v ' + fl + ' -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP|PRES|TMP):(surface)|:(UGRD|VGRD):(10 m above ground):" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'nc'
	#os.system(ll) 

'''
#%% convert from  grib2 to nc files
for fl in grb2_files:
	print(fl)
	ll = 'wgrib2 -v ' + fl + ' -match ":(UGRD|VGRD):850 mb|:(UGRD|VGRD):200 mb:" -netcdf ' + Dir_output + fl.split('/')[-1][0:-5] + '_shear' + '.nc'

	os.system(ll) 
'''
'''
os.system('wgrib2 -v /scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/dorian05l.2019083018.hwrfprs.storm.0p015.f003.grb2 -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf /home/Maria.Aristizabal/Dorian_2019/HWRF_nc_files/dorian05l.2019083018.hwrfprs.storm.0p015.f003.nc')
'''
