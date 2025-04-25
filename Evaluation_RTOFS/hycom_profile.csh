#!/bin/csh
set echo
setenv DH /scratch2/NCEPDEV/marine/Zulema.Garraffo/HYCOM-tools/bin
setenv Dh /scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_dell/hycom/data/data_2020042400
setenv Df /scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs/hycom/fix

#change directory and link the archive from which you want the profile
#setenv D  ... your directory
#cd $D

setenv D ./data
mkdir -p $D
cd $D


ln -sf $Dh/archv.2020_115_00.[a,b] .
ln -s $Df/regional.grid.? . 
setenv lon -90
setenv lat 25
setenv i `$DH/hycom_lonlat2ij $lon $lat  | awk '{print $1}'`
setenv j `$DH/hycom_lonlat2ij $lon $lat  | awk '{print $2}'`
$DH/hycom_profile archv.2020_115_00.a $i $j > profile_${lon}_${lat}.out
