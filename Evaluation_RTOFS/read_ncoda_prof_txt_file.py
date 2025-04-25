#!/usr/bin/env python

# Usage
# python read_ncoda_prof_txt_file.py ncoda_file min_lon max_lon min_lat max_lat var_to_read 

# Example:
# python read_ncoda_prof_txt_file.py '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/profiles_qc/prof_argo_rpt.2021101600.txt' -160 -150 20 30 'Temp'   

# Or from an ipython console
# run read_ncoda_prof_txt_file.py '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_DA/parallel3/rtofs.20211017_ab/profiles_qc/prof_argo_rpt.2021101600.txt' -160 -150 20 30 'Temp'

 
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

    lon, lat = get_lat_lot_ncoda_txt_file(ncoda_file)
 
    oklon = [n for n in np.arange(len(lon)) if float(lon[n]) >= min_lon and float(lon[n]) <= max_lon]
    oklonlat = [n for n in oklon if float(lat[n]) >= min_lat and float(lat[n]) <= max_lat]

    for n in oklonlat:
        target_lon = lon[n]
        target_lat = lat[n]
        lon_argo, lat_argo, DTG, Rcpt, Sign, depth_argo, temp_argo = get_prof_ncoda_txt_file(ncoda_file,target_lon,target_lat,var_to_read)

        temp_argo = np.asarray(temp_argo)
        depth_argo = np.asarray(depth_argo)

        plt.figure(figsize=(5,7))
        plt.plot(temp_argo,-depth_argo,'.-',label = 'Argo')
        plt.ylim([-1800,0])
        plt.title('Temperature Argo Sign = '+ Sign + '\n Lon = ' + lon_argo + ' ,Lat = ' + lat_argo + '\n DTG = ' + 
DTG + ' Rcpt = ' + Rcpt ) 
        plt.legend()

