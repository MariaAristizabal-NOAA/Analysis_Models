#!/usr/bin/env python

# Usage
# python read_ncoda_prof_txt_file.py ncoda_file min_lon max_lon min_lat max_lat var_to_read rtofs_ab_file1 rtofs_ab_file2 rtofs_ab_file3 rtofs_grid_file rtofs_depth_file

# Example:

# From an ipython console
# run RTOFS_vs_Argo.py /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/profiles_qc/prof_argo_rpt.2021101600.txt -160 -150 20 30 Temp /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/rtofs_run4/rtofs.20211017/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/operational/20211017/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/GRID_DEPTH/regional.grid /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/GRID_DEPTH/regional.depth  

# run RTOFS_vs_Argo.py /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/profiles_qc/prof_argo_rpt.2021101600.txt -100 -90 -5 0 Temp /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/rtofs_run4/rtofs.20211017/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/operational/20211017/rtofs_glo.t00z.n00.archv /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/GRID_DEPTH/regional.grid /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/GRID_DEPTH/regional.depth  

####################################################################################
def get_lat_lot_ncoda_txt_file(ncoda_file):

    file = ncoda_file
    ff = open(file,'r')
    f = ff.readlines()

    lat = []
    lon = []
    for l in f:
        if len(l.split()) > 3 and l.split()[2]=='Lat=':
            lat.append(l.split()[3])
            lon.append(l.split()[5])

    return lon, lat

####################################################################################
def get_prof_ncoda_txt_file(ncoda_file,target_lon,target_lat,var_to_read):

    file = ncoda_file
    ff = open(file,'r')
    f = ff.readlines()

    depth_vector = []
    temp_vector = []
    for n,l in enumerate(f):
        if len(l.split()) > 3 and l.split()[3]==target_lat and l.split()[5]==target_lon:
            lat = l.split()[3]
            lon = l.split()[5]
            DTG = l.split()[9]
            Rcpt =  l.split()[11]
            Sign =  l.split()[13]
            print('lon =',lon)
            print('lat =',lat)
            print('DTG =',DTG)
            print('Sign =',Sign)
            print(' ')
            if f[n+5].split()[1] == var_to_read:
                for lines in f[n+6:]:
                    if len(lines.split()) != 0:
                        depth_vector.append(float(lines.split()[0]))
                        temp_vector.append(float(lines.split()[1]))
                    else:
                        break

    return lon, lat, DTG, Rcpt, Sign, depth_vector, temp_vector

####################################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
def find_grid_position_hycom(lon,lat,target_lon,target_lat):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),\
                                             left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklon,oklat

##################################################################################
def get_profile_rtofs_ab_file_desn_layers(nz,oklon,oklat,ab_file):

    layers = np.arange(0,nz)
    target_temp_rtofs = np.empty((nz))
    target_temp_rtofs[:] = np.nan
    target_z_rtofs = np.empty((nz))
    target_z_rtofs[:] = np.nan
    timeRTOFS = []

    lines = [line.rstrip() for line in open(ab_file + '.b')]
    time_stamp = lines[-1].split()[2]
    hycom_days = lines[-1].split()[3]
    tzero = datetime(1901,1,1,0,0)
    timeRT = tzero+timedelta(float(hycom_days))
    timeRTOFS.append(timeRT)
    timestampRTOFS = mdates.date2num(timeRT)

    ztmp = readVar(ab_file,'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
    target_ztmp = ztmp[oklat,oklon]
    for lyr in tuple(layers):
        print(lyr)
        temp_RTOFS = readVar(ab_file,'archive','temp',[lyr+1])
        target_temp_rtofs[lyr] = temp_RTOFS[oklat,oklon]
        dp = readVar(ab_file,'archive','thknss',[lyr+1])/9806
        target_ztmp = np.append(target_ztmp,dp[oklat,oklon])

    target_z_rtofs = np.cumsum(target_ztmp[0:-1]) + np.diff(np.cumsum(target_ztmp))/2
    

    time_rtofs = np.asarray(timeRTOFS)

    return target_temp_rtofs, target_z_rtofs, time_rtofs

##################################################################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import glob
import sys

folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readBinz, readgrids, readdepth, readVar

sys.path.append(folder_myutils)
from my_models_utils import geo_coord_to_HYCOM_coord

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

##################################################################################
import sys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    print(f"Arguments count: {len(sys.argv)}")
    for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>6}: {arg}")

    ncoda_file = sys.argv[1]
    min_lon = float(sys.argv[2])
    max_lon = float(sys.argv[3])
    min_lat = float(sys.argv[4])
    max_lat = float(sys.argv[5])
    var_to_read = sys.argv[6]
    rtofs_ab_file1 = sys.argv[7]
    rtofs_ab_file2 = sys.argv[8]
    rtofs_ab_file3 = sys.argv[9]
    rtofs_grid_file = sys.argv[10]
    rtofs_depth_file = sys.argv[11]

    lon, lat = get_lat_lot_ncoda_txt_file(ncoda_file)

    oklon = [n for n in np.arange(len(lon)) if float(lon[n]) >= min_lon and float(lon[n]) <= max_lon]
    oklonlat = [n for n in oklon if float(lat[n]) >= min_lat and float(lat[n]) <= max_lat]

    for n in oklonlat:
        target_lon = lon[n]
        target_lat = lat[n]
        lon_argo, lat_argo, DTG, Rcpt, Sign, depth_argo, temp_argo = get_prof_ncoda_txt_file(ncoda_file,target_lon,target_lat,var_to_read)

        temp_argo = np.asarray(temp_argo)
        depth_argo = np.asarray(depth_argo)

        ##################################################################################
        #%% Reading HYCOM grid from ab files
        print('Retrieving coordinates from RTOFS')
        lines_grid = [line.rstrip() for line in open(rtofs_grid_file +'.b')]
        lon_rtofs = np.array(readgrids(rtofs_grid_file,'plon:',[0]))
        lat_rtofs = np.array(readgrids(rtofs_grid_file,'plat:',[0]))
        depth_rtofs_ab = np.asarray(readdepth(rtofs_depth_file,'depth'))

        ##################################################################################
        #%% find grid point
        target_lonH,target_latH = geo_coord_to_HYCOM_coord(float(target_lon),float(target_lat))

        oklon_ab, oklat_ab = find_grid_position_hycom(lon_rtofs,lat_rtofs,target_lonH,target_latH)

        ##################################################################################
        #%% get target temp prof from ab files
        nz = 41
        target_temp_rtofs_par3, target_z_rtofs_par3, time_rtofs_par3 = get_profile_rtofs_ab_file_desn_layers(nz,oklon_ab,oklat_ab,rtofs_ab_file1)

        target_temp_rtofs_run4, target_z_rtofs_run4, time_rtofs_run4 = get_profile_rtofs_ab_file_desn_layers(nz,oklon_ab,oklat_ab,rtofs_ab_file2)

        target_temp_rtofs_oper, target_z_rtofs_oper, time_rtofs_oper = get_profile_rtofs_ab_file_desn_layers(nz,oklon_ab,oklat_ab,rtofs_ab_file3)

        ##################################################################################
        #%% Plot profiles

        plt.figure(figsize=(5,7))
        plt.plot(temp_argo,-depth_argo,'.-',color='b',label = 'Argo')
        plt.plot(target_temp_rtofs_oper,-target_z_rtofs_oper,'.-',color='magenta',label = 'Operational '+str(time_rtofs_run4[0])[0:13])
        plt.plot(target_temp_rtofs_par3,-target_z_rtofs_par3,'.-',color='orange',label = 'Parallel3 '+str(time_rtofs_par3[0])[0:13])
        plt.plot(target_temp_rtofs_run4,-target_z_rtofs_run4,'.-',color='green',label = 'Run4 '+str(time_rtofs_run4[0])[0:13])
        plt.ylim([-2000,0])
        plt.legend()
        plt.title('Temperature Argo Sign = ' + Sign + '\n DTG = ' + DTG + ' Rcpt = ' + Rcpt + '\n Lon_argo = ' + lon_argo + ' ,Lat_argo = ' + lat_argo + '\n Lon_rtofs = ' + str(np.round(lon_rtofs[oklat_ab,oklon_ab]-360,4)) + ' Lat_rtofs = ' + str(np.round(lat_rtofs[oklat_ab,oklon_ab],4)))



