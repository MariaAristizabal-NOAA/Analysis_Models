"""
  Plot SST using SST from 2Dvar analysis fields

"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime

sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/MyPython/ncoda_utils')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/MyPython/draw_map')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/MyPython')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/MyPython/hycom_utils')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/MyPython/hausdorff')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/RTOFS_utils_Dmitry/validation_rtofs')

from mod_utils_fig import bottom_text
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
import mod_gulfstream as mgulf
import mod_time as mtime
import mod_read_ncoda as mncoda
importlib.reload(mncoda)
importlib.reload(mmisc)
importlib.reload(mgulf)
importlib.reload(mtime)

rdate0  = '20240402'
sfx     = 'n-24'

# Input directories with RTOFS and 2DVAR ssh:
#pthbase  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/'
#pthrtofs = f'{pthbase}rtofs/hycom/rtofs.{rdate0}/'
#pthanls  = f'{pthbase}rtofs.prod/rtofs.{rdate0}/glbl_var/restart/'
#pthgrid  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

pthrtofs = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS/rtofs.20240402/'
pthanls = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS/rtofs.20240402/ncoda/glbl_var/restart/'
pthgrid = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/RTOFS_fix/GRID_DEPTH_global/'

fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'
fina   = pthrtofs + fhcm + '.a'
finb   = pthrtofs + fhcm + '.b'

yrday  = mtime.rdate2jday(rdate0)
dnmb0  = mtime.rdate2datenum(rdate0)
# Date of plotted fields:
if sfx == 'n-24':
  dnmbP = dnmb0-1
  hr = 0
elif sfx[0] == 'f':
  hr = int(sfx[1:])
  dnmbP = dnmb0+float(hr)/24.
rdateP = mtime.dnumb2rdate(dnmbP, ihours=False)

dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]


huge = 1.e20
rg   = 9806.

# RTOFS grid topo
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LONr, LATr, HHr = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# 2D analysis grid:
#/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraD5b/rtofs.20231215/glbl_var/restart/
# 3241x2441, of .11107 degrees - these are anomalies
# The mean SSH is binary ieee, dimensions 2881x1441, i
# regular grid, size 0.125, starting at long=0 and lat=-90
# 32bit ieee unformatted, values -999 over land. 
#flon  = 'grdlon_sfc_1o3241x2441_2023121400_0000_datafld'
#flat  = 'grdlat_sfc_1o3241x2441_2023121400_0000_datafld'
#ftopo = 'depths_sfc_1o3241x2441_2023123100_0000_datafld'
import mod_valid_utils as mvalid
importlib.reload(mvalid)
fdlon  = mvalid.find_fileindir(pthanls, 'grdlon_sfc_1o3241x2441_*datafld')
fdlat  = mvalid.find_fileindir(pthanls, 'grdlat_sfc_1o3241x2441_*datafld')
fdtopo = mvalid.find_fileindir(pthanls, 'depths_sfc_1o3241x2441_*datafld')

jgrd  = 2441
igrd  = 3241
#jgrd  = 3298
#igrd  = 4500
LON   = mncoda.read_2Dbin(fdlon, igrd, jgrd)
LAT   = mncoda.read_2Dbin(fdlat, igrd, jgrd) 
HH    = mncoda.read_2Dbin(fdtopo, igrd, jgrd)
if np.min(HH) > -1.e-30:
  HH = -HH

HH    = np.where(HH>=-1.e-30, 100., HH)
IDM   = HH.shape[1]
JDM   = HH.shape[0]

LON = np.where(LON > 180, LON-360, LON)

# Read SSH 2D analysis:
#fanls  = f'seahgt_sfc_1o3241x2441_{rdateP}00_0000_analfld'
fanls  = f'seatmp_sfc_1o3241x2441_{rdateP}00_0000_analfld'
fdanls = os.path.join(pthanls,fanls)
SST    = mncoda.read_2Dbin(fdanls, igrd, jgrd)
SST    = np.where(SST <= -998, np.nan, SST)

# Plot fields:
kw = dict(levels=np.arange(-3,33+0.1,2))
plt.figure(figsize=(10,5))
plt.contourf(SST,cmap='Spectral_r',**kw,extend='both')
plt.colorbar(extendrect=True)
plt.title(fanls ,fontsize=18)




