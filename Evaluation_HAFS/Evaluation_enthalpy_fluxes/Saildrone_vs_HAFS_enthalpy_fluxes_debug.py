#%% User input
# forecasting cycle to be used

# Lee
'''
cycle = '2023090706'
storm_num = '13'
basin = 'al'
storm_id = '13l'
storm_name = 'lee'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_update26feb2024.nc'
'''

cycle = '2024081612'
storm_num = '05'
basin = 'al'
storm_id = '05l'
storm_name = 'Ernesto'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2024/sd1031_hurricane_2024.nc'

#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/sd1064_hurricane_2023.nc'
#url_saildrone_adcp = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2023/From_Greg_Foltz/sd1069_hurricane_2023_adcp_cbe4_c409_047b_U1696433734196.nc'

folder_myutils= '/home/Maria.Aristizabal/Utils/'

# folder CORE-algorithm
folder_COARE = '/home/Maria.Aristizabal/Analysis/COARE-algorithm/Python/COARE3.6'

#############################################################################
import xarray as xr
import cfgrib
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob


sys.path.append(folder_COARE)
#from coare36vn_zrf_et import coare36vn_zrf_et
from coare36vn_zrf_et_debug import coare36vn_zrf_et

#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#------------------------------------------------------------------------------

def psit_26(zeta = None):
    # computes temperature structure function
    dzeta = np.minimum(50,0.35 * zeta)
    psi = - ((1 + 0.6667 * zeta) ** 1.5 + np.multiply(0.6667 * (zeta - 14.28),np.exp(- dzeta)) + 8.525)
    k = np.array(np.where(zeta < 0))
    x = (1 - 15 * zeta[k]) ** 0.5
    psik = 2 * np.log((1 + x) / 2)
    x = (1 - 34.15 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi

#------------------------------------------------------------------------------

def psiu_26(zeta = None):
    # computes velocity structure function
    dzeta = np.minimum(50,0.35 * zeta)
    a = 0.7
    b = 3 / 4
    c = 5
    d = 0.35
    psi = - (a * zeta + np.multiply(b * (zeta - c / d),np.exp(- dzeta)) + b * c / d)
    k = np.array(np.where(zeta < 0))
    x = (1 - 15 * zeta[k]) ** 0.25
    psik = 2 * np.log((1 + x) / 2) + np.log((1 + np.multiply(x,x)) / 2) - 2 * np.arctan(x) + 2 * np.arctan(1)
    x = (1 - 10.15 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi
#------------------------------------------------------------------------------

def psiu_40(zeta = None):
    # computes velocity structure function
    dzeta = np.minimum(50,0.35 * zeta)
    a = 1
    b = 3 / 4
    c = 5
    d = 0.35
    psi = - (a * zeta + np.multiply(b * (zeta - c / d),np.exp(- dzeta)) + b * c / d)
    k = np.array(np.where(zeta < 0))
    x = (1 - 18 * zeta[k]) ** 0.25
    psik = 2 * np.log((1 + x) / 2) + np.log((1 + np.multiply(x,x)) / 2) - 2 * np.arctan(x) + 2 * np.arctan(1)
    x = (1 - 10 * zeta[k]) ** 0.3333
    psic = 1.5 * np.log((1 + x + x ** 2) / 3) - np.sqrt(3) * np.arctan((1 + 2 * x) / np.sqrt(3)) + 4 * np.arctan(1) / np.sqrt(3)
    f = zeta[k] ** 2.0 / (1 + zeta[k] ** 2)
    psi[k] = np.multiply((1 - f),psik) + np.multiply(f,psic)
    return psi
#------------------------------------------------------------------------------

def bucksat(T = None,P = None,Tf = None):
    # computes saturation vapor pressure [mb]
    # given T [degC] and P [mb] Tf is freezing pt
    exx = np.multiply(np.multiply(6.1121,np.exp(np.multiply(17.502,T) / (T + 240.97))),(1.0007 + np.multiply(3.46e-06,P)))
    ii = np.array(np.where(T < Tf))
    if np.size(ii) != 0:
        exx[ii] = np.multiply(np.multiply((1.0003 + 4.18e-06 * P[ii]),6.1115),np.exp(np.multiply(22.452,T[ii]) / (T[ii] + 272.55)))

    return exx
#------------------------------------------------------------------------------

def qsat26sea(T = None,P = None,Ss = None,Tf = None):
    # computes surface saturation specific humidity [g/kg]
    # given T [degC] and P [mb]
    ex = bucksat(T,P,Tf)
    fs = 1 - 0.02 * Ss / 35
    es = np.multiply(fs,ex)
    qs = 622 * es / (P - 0.378 * es)
    return qs

#------------------------------------------------------------------------------

def qsat26air(T = None,P = None,rh = None):
    # computes saturation specific humidity [g/kg]
    # given T [degC] and P [mb]
    Tf = 0
    es = bucksat(T,P,Tf)
    em = np.multiply(0.01 * rh,es)
    q = 622 * em / (P - 0.378 * em)
    return q,em

#------------------------------------------------------------------------------

def grv(lat = None):
    # computes g [m/sec^2] given lat in deg
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 2.32718e-05
    c3 = 1.262e-07
    c4 = 7e-10
    phi = lat * np.pi / 180
    x = np.sin(phi)
    g = gamma * (1 + c1 * x ** 2 + c2 * x ** 4 + c3 * x ** 6 + c4 * x ** 8)
    return g
#------------------------------------------------------------------------------

def RHcalc(T = None,P = None,Q = None,Tf = None):
    # computes relative humidity given T,P, & Q
    es = np.multiply(np.multiply(6.1121,np.exp(np.multiply(17.502,T) / (T + 240.97))),(1.0007 + np.multiply(3.46e-06,P)))
    ii = np.array(np.where(T < Tf))
    if np.size(ii) != 0:
        es[ii] = np.multiply(np.multiply(6.1115,np.exp(np.multiply(22.452,T[ii]) / (T[ii] + 272.55))),(1.0003 + 4.18e-06 * P[ii]))
    em = np.multiply(Q,P) / (np.multiply(0.378,Q) + 0.622)
    RHrf = 100 * em / es
    return RHrf

#------------------------------------------------------------------------------

def albedo_vector(sw_dn = None,jd = None,lon = None,lat = None,eorw = None):
    #  Computes transmission and albedo from downwelling sw_dn using
    #  lat   : latitude in degrees (positive to the north)
    #  lon   : longitude in degrees (positive to the west)
    #  jd    : yearday
    #  sw_dn : downwelling solar radiation measured at surface
    #  eorw  : 'E' if longitude is positive to the east, or 'W' if otherwise

    # updates:
    #   20-10-2021: ET vectorized function

    if eorw == 'E':
        #     disp('lon is positive to east so negate for albedo calculation');
        lon = - lon
    elif eorw == 'W':
        #     disp('lon is already positive to west so go ahead with albedo calculation');
        pass
    else:
        print('please provide sign information on whether lon is deg E or deg W')

    alb = np.full([np.size(sw_dn)],np.nan)
    lat = lat * np.pi / 180
    lon = lon * np.pi / 180
    SC = 1380
    utc = (jd - np.fix(jd)) * 24
    h = np.pi * utc / 12 - lon
    declination = 23.45 * np.cos(2 * np.pi * (jd - 173) / 365.25)
    solarzenithnoon = (lat * 180 / np.pi - declination)
    solaraltitudenoon = 90 - solarzenithnoon
    sd = declination * np.pi / 180
    gamma = 1
    gamma2 = gamma * gamma

    sinpsi = np.multiply(np.sin(lat),np.sin(sd)) - np.multiply(np.multiply(np.cos(lat),np.cos(sd)),np.cos(h))
    psi = np.multiply(np.arcsin(sinpsi),180) / np.pi
    solarmax = np.multiply(SC,sinpsi) / gamma2
    #solarmax=1380*sinpsi*(0.61+0.20*sinpsi);

    T = np.minimum(2,sw_dn / solarmax)

    Ts = np.arange(0,1+0.05,0.05)
    As = np.arange(0,90+2,2)

    a = np.array([[0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.062,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.072,0.07,0.068,0.065,0.065,0.063,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.06,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.087,0.083,0.079,0.073,0.07,0.068,0.066,0.065,0.064,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.115,0.108,0.098,0.086,0.082,0.077,0.072,0.071,0.067,0.067,0.065,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.061,0.061,0.061,0.061,0.06,0.059,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059],[0.163,0.145,0.13,0.11,0.101,0.092,0.084,0.079,0.072,0.072,0.068,0.067,0.064,0.063,0.062,0.061,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058],[0.235,0.198,0.174,0.15,0.131,0.114,0.103,0.094,0.083,0.08,0.074,0.074,0.07,0.067,0.065,0.064,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058],[0.318,0.263,0.228,0.192,0.168,0.143,0.127,0.113,0.099,0.092,0.084,0.082,0.076,0.072,0.07,0.067,0.065,0.064,0.062,0.062,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.057,0.058,0.058,0.058,0.058,0.057,0.057,0.057],[0.395,0.336,0.29,0.248,0.208,0.176,0.151,0.134,0.117,0.107,0.097,0.091,0.085,0.079,0.075,0.071,0.068,0.067,0.065,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.058,0.058,0.058,0.057,0.057,0.057,0.057,0.057,0.057,0.057,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.055],[0.472,0.415,0.357,0.306,0.252,0.21,0.176,0.154,0.135,0.125,0.111,0.102,0.094,0.086,0.081,0.076,0.072,0.071,0.068,0.066,0.065,0.063,0.062,0.061,0.06,0.059,0.058,0.057,0.057,0.057,0.056,0.055,0.055,0.055,0.055,0.055,0.055,0.054,0.053,0.054,0.053,0.053,0.054,0.054,0.053,0.053],[0.542,0.487,0.424,0.36,0.295,0.242,0.198,0.173,0.15,0.136,0.121,0.11,0.101,0.093,0.086,0.081,0.076,0.073,0.069,0.067,0.065,0.064,0.062,0.06,0.059,0.058,0.057,0.056,0.055,0.055,0.054,0.053,0.053,0.052,0.052,0.052,0.051,0.051,0.05,0.05,0.05,0.05,0.051,0.05,0.05,0.05],[0.604,0.547,0.498,0.407,0.331,0.272,0.219,0.185,0.16,0.141,0.127,0.116,0.105,0.097,0.089,0.083,0.077,0.074,0.069,0.066,0.063,0.061,0.059,0.057,0.056,0.055,0.054,0.053,0.053,0.052,0.051,0.05,0.05,0.049,0.049,0.049,0.048,0.047,0.047,0.047,0.046,0.046,0.047,0.047,0.046,0.046],[0.655,0.595,0.556,0.444,0.358,0.288,0.236,0.19,0.164,0.145,0.13,0.119,0.107,0.098,0.09,0.084,0.076,0.073,0.068,0.064,0.06,0.058,0.056,0.054,0.053,0.051,0.05,0.049,0.048,0.048,0.047,0.046,0.046,0.045,0.045,0.045,0.044,0.043,0.043,0.043,0.042,0.042,0.043,0.042,0.042,0.042],[0.693,0.631,0.588,0.469,0.375,0.296,0.245,0.193,0.165,0.145,0.131,0.118,0.106,0.097,0.088,0.081,0.074,0.069,0.065,0.061,0.057,0.055,0.052,0.05,0.049,0.047,0.046,0.046,0.044,0.044,0.043,0.042,0.042,0.041,0.041,0.04,0.04,0.039,0.039,0.039,0.038,0.038,0.038,0.038,0.038,0.038],[0.719,0.656,0.603,0.48,0.385,0.3,0.25,0.193,0.164,0.145,0.131,0.116,0.103,0.092,0.084,0.076,0.071,0.065,0.061,0.057,0.054,0.051,0.049,0.047,0.045,0.043,0.043,0.042,0.041,0.04,0.039,0.039,0.038,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.034,0.034,0.034,0.034,0.034,0.034],[0.732,0.67,0.592,0.474,0.377,0.291,0.246,0.19,0.162,0.144,0.13,0.114,0.1,0.088,0.08,0.072,0.067,0.062,0.058,0.054,0.05,0.047,0.045,0.043,0.041,0.039,0.039,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.033,0.032,0.032,0.032,0.031,0.031,0.031,0.03,0.03,0.03,0.03,0.03],[0.73,0.652,0.556,0.444,0.356,0.273,0.235,0.188,0.16,0.143,0.129,0.113,0.097,0.086,0.077,0.069,0.064,0.06,0.055,0.051,0.047,0.044,0.042,0.039,0.037,0.035,0.035,0.035,0.034,0.033,0.033,0.032,0.032,0.032,0.029,0.029,0.029,0.029,0.028,0.028,0.028,0.028,0.027,0.027,0.028,0.028],[0.681,0.602,0.488,0.386,0.32,0.252,0.222,0.185,0.159,0.142,0.127,0.111,0.096,0.084,0.075,0.067,0.062,0.058,0.054,0.05,0.046,0.042,0.04,0.036,0.035,0.033,0.032,0.032,0.031,0.03,0.03,0.03,0.03,0.029,0.027,0.027,0.027,0.027,0.026,0.026,0.026,0.026,0.026,0.026,0.026,0.026],[0.581,0.494,0.393,0.333,0.288,0.237,0.211,0.182,0.158,0.141,0.126,0.11,0.095,0.083,0.074,0.066,0.061,0.057,0.053,0.049,0.045,0.041,0.039,0.034,0.033,0.032,0.031,0.03,0.029,0.028,0.028,0.028,0.028,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.453,0.398,0.342,0.301,0.266,0.226,0.205,0.18,0.157,0.14,0.125,0.109,0.095,0.083,0.074,0.065,0.061,0.057,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.425,0.37,0.325,0.29,0.255,0.22,0.2,0.178,0.157,0.14,0.122,0.108,0.095,0.083,0.074,0.065,0.061,0.056,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.026,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025]])

    if T.size==1:   ### for single value function
        Tchk = np.abs(Ts - T)
        i = np.array(np.where(Tchk == Tchk.min()))
        if psi < 0:
            alb = np.array([0])
            solarmax = np.array([0])
            T = np.array([0])
            j = np.array([0])
            psi = np.array([0])
        else:
            Achk = np.abs(As - psi)
            j = np.array(np.where(Achk == Achk.min()))
            szj = j.shape
            if szj[0] > 0:
                alb = a[i,j].flatten()
            else:
                #       print('no j found, not assigning alb to anything');
                pass
    else:  ### for vectorized function
        for k in np.arange(0,np.size(sinpsi)).reshape(-1):
            Tchk = np.abs(Ts - T[k])
            i = np.array(np.where(Tchk == Tchk.min()))
            if psi[k] < 0:
                alb[k] = 0
                solarmax[k] = 0
                T[k] = 0
                j = 0
                psi[k] = 0
            else:
                Achk = np.abs(As - psi[k])
                j = np.array(np.where(Achk == Achk.min()))
                szj = j.shape
                if szj[0] > 0:
                    alb[k] = a[i,j]
                else:
                    #       disp('no j found, not assigning alb to anything');
                    pass

    #disp([num2str(jd) '  ' num2str(sw_dn) '  ' num2str(alb) '  ' num2str(T) '  ' num2str(i) '  ' num2str(j)])
    return alb,T,solarmax,psi

#############################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#############################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

#############################################################################

url = url_saildrone
gdata = xr.open_dataset(url)#,decode_times=False)

# Variables in the same order as required by CORE

# u = water-relative wind speed magnitude (m/s) at height zu (m)
wind_from_mean = np.asarray(gdata.WIND_FROM_MEAN)
wind_speed_mean = np.asarray(gdata.WIND_SPEED_MEAN)
water_current_speed_mean = np.asarray(gdata.WATER_CURRENT_SPEED_MEAN)
water_current_direccion_mean = np.asarray(gdata.WATER_CURRENT_DIRECTION_MEAN)

# zu height of wind speed observations
zu = gdata.WIND_SPEED_MEAN.installed_height

# t = air temperature (degC) at height zt (m)
temp_air_mean = np.asarray(gdata.TEMP_AIR_MEAN)

# zt height of air temp. observations
zt = gdata.TEMP_AIR_MEAN.installed_height

# rh = relative humidity (#) at height zq (m)
rh_mean = np.asarray(gdata.RH_MEAN)

# zq height of relative humidity observations
zq = gdata.RH_MEAN.installed_height

# P = sea level air pressure (mb)
baro_pres_mean = np.asarray(gdata.BARO_PRES_MEAN)

# ts = seawater temperature (degC)
temp_sb37_mean = np.asarray(gdata.TEMP_SBE37_MEAN)

# sw_dn = downward (positive) shortwave radiation (W/m^2)
Sw_dn = np.tile(200,len(temp_sb37_mean)) # default value inside COARE

# lat = latitude defined positive to north
latitude = np.asarray(gdata.latitude)

# lw_dn = downward (positive) longwave radiation (W/m^2)
Lw_dn = 400-1.6*np.abs(latitude) # default value inside COARE

# lon = longitude defined positive to east
longitude = np.asarray(gdata.longitude)

# jd = year day or julian day, where day Jan 1 00:00 UTC = 0
time = np.asarray(gdata.time)
timestamps = mdates.date2num(time)
times = np.asarray(mdates.num2date(timestamps))
JD = np.asarray([t.timetuple().tm_yday for t in times])

# zi = PBL height (m) (default or typical value = 600m)
zi = 600

# rain = rain rate (mm/hr)
rain = 0

#  Ss = sea surface salinity (PSU)
sal_sb37_mean = np.asarray(gdata.SAL_SBE37_MEAN)

# cp = phase speed of dominant waves (m/s) computed from peak period
wave_dominant_period = np.asarray(gdata.WAVE_DOMINANT_PERIOD)

# sigH = significant wave height (m)
wave_significant_height = np.asarray(gdata.WAVE_SIGNIFICANT_HEIGHT)

dataset_id = gdata.drone_id

##############################################################################
# Obtain enthalpy fluxes using COARE algorithm

oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

u = wind_speed_mean[oktimeg]
ok = np.isfinite(u)

u = u[ok]
zu = zu
t = temp_air_mean[oktimeg][ok] 
zt = zt
rh = rh_mean[oktimeg][ok] 
zq = zq
P = baro_pres_mean[oktimeg][ok] 
ts = temp_sb37_mean[oktimeg][ok] 
sw_dn = Sw_dn[oktimeg][ok]
lw_dn = Lw_dn[oktimeg][ok]
lat = latitude[oktimeg][ok]
lon = longitude[oktimeg][ok]
jd = JD[oktimeg][ok]
zi = zi
rain = rain
Ss = sal_sb37_mean[oktimeg][ok] 
cp = wave_dominant_period[oktimeg][ok] 
sigH = wave_significant_height[oktimeg][ok] 

###########################################################
A = coare36vn_zrf_et(u, zu, t, zt, rh, zq, P, ts, sw_dn, lw_dn, lat, lon, jd, zi, rain, Ss, cp, sigH, zrf_u=10.0, zrf_t=10.0, zrf_q=10.0)

jcoolx = 1

if u.size ==1 and t.size ==1:
    u = np.copy(np.asarray([u], dtype=float)).flatten()
    zu = np.copy(np.asarray([zu], dtype=float)).flatten()
    t = np.copy(np.asarray([t], dtype=float)).flatten()
    zt = np.copy(np.asarray([zt], dtype=float)).flatten()
    rh = np.copy(np.asarray([rh], dtype=float)).flatten()
    zq = np.copy(np.asarray([zq], dtype=float)).flatten()
    P = np.copy(np.asarray([P], dtype=float)).flatten()
    ts = np.copy(np.asarray([ts], dtype=float)).flatten()
    sw_dn = np.copy(np.asarray([sw_dn], dtype=float)).flatten()
    lw_dn = np.copy(np.asarray([lw_dn], dtype=float)).flatten()
    lat = np.copy(np.asarray([lat], dtype=float)).flatten()
    lon = np.copy(np.asarray([lon], dtype=float)).flatten()
    jd = np.copy(np.asarray([jd], dtype=float)).flatten()
    zi = np.copy(np.asarray([zi], dtype=float)).flatten()
    rain = np.copy(np.asarray([rain], dtype=float)).flatten()
    Ss = np.copy(np.asarray([Ss], dtype=float)).flatten()
    zrf_u = np.copy(np.asarray([zrf_u], dtype=float)).flatten()
    zrf_t = np.copy(np.asarray([zrf_t], dtype=float)).flatten()
    zrf_q = np.copy(np.asarray([zrf_q], dtype=float)).flatten()

N = np.size(u)
jcool = jcoolx * np.ones(N)

if cp is not None and cp.size==1:
    cp = np.copy(np.asarray([cp], dtype=float)).flatten()
elif cp is None:
    cp = np.nan * np.ones(N)

if sigH is not None and sigH.size==1:
    sigH = np.copy(np.asarray([sigH], dtype=float)).flatten()
elif sigH is None:
    sigH = np.nan * np.ones(N)

us = 0 * u
# convert rh to specific humidity after accounting for salt effect on freezing
# point of water
Tf = - 0.0575 * Ss + 0.00171052 * Ss ** 1.5 - np.multiply(0.0002154996 * Ss,Ss)
Qs = qsat26sea(ts,P,Ss,Tf) / 1000
P_tq = P - (0.125 * zt)
Q,Pv = qsat26air(t,P_tq,rh)

Q = Q / 1000
ice = np.zeros(N)
iice = np.array(np.where(ts < Tf))
ice[iice] = 1
jcool[iice] = 0
zos = 0.0005
#***********  set constants ***********************************************
zref = 10
Beta = 1.2
von = 0.4
fdg = 1.0
T2K = 273.16
grav = grv(lat)
#***********  air constants ***********************************************
Rgas = 287.1
Le = (2.501 - 0.00237 * ts) * 1000000.0
cpa = 1004.67
cpv = cpa * (1 + 0.84 * Q)
rhoa = P_tq * 100.0 / (np.multiply(Rgas * (t + T2K),(1 + 0.61 * Q)))
# Pv is the partial pressure due to wate vapor in mb
rhodry = (P_tq - Pv) * 100.0 / (Rgas * (t + T2K))
visa = 1.326e-05 * (1 + np.multiply(0.006542,t) + 8.301e-06 * t ** 2 - 4.84e-09 * t ** 3)
lapse = grav / cpa

#***********  cool skin constants  ***************************************
### includes salinity dependent thermal expansion coeff for water
tsw = ts
ii = np.array(np.where(ts < Tf))
if np.size(ii) != 0:
    tsw[ii] = Tf[ii]
Al35 = 2.1e-05 * (tsw + 3.2) ** 0.79
# Al0 = (2.2 * real((tsw - 1) ** 0.82) - 5) * 1e-05
Al0_i=(tsw - 1) ** 0.82
Al0 = (2.2 * Al0_i.real - 5) * 1e-05
Al = Al0 + np.multiply((Al35 - Al0),Ss) / 35
###################
bets = 0.00075
be = bets * Ss
####  see "Computing the seater expansion coefficients directly from the
####  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
cpw = 4000
rhow = 1022
visw = 1e-06
tcw = 0.6
bigc = 16 * grav * cpw * (rhow * visw) ** 3.0 / (tcw ** 2 * rhoa ** 2)
wetc = np.multiply(0.622 * Le,Qs) / (Rgas * (ts + T2K) ** 2)

#alb,T_sw,solarmax_sw,psi_sw = albedo_vector(sw_dn,jd,lon,lat,eorw='E')
eorw='E'

####################################################
if eorw == 'E':
    #     disp('lon is positive to east so negate for albedo calculation');
    lon = - lon
elif eorw == 'W':
    #     disp('lon is already positive to west so go ahead with albedo calculation');
    pass
else:
    print('please provide sign information on whether lon is deg E or deg W')

alb = np.full([np.size(sw_dn)],np.nan)
lat = lat * np.pi / 180
lon = lon * np.pi / 180
SC = 1380
utc = (jd - np.fix(jd)) * 24
h = np.pi * utc / 12 - lon
declination = 23.45 * np.cos(2 * np.pi * (jd - 173) / 365.25)
solarzenithnoon = (lat * 180 / np.pi - declination)
solaraltitudenoon = 90 - solarzenithnoon
sd = declination * np.pi / 180
gamma = 1
gamma2 = gamma * gamma
sinpsi = np.multiply(np.sin(lat),np.sin(sd)) - np.multiply(np.multiply(np.cos(lat),np.cos(sd)),np.cos(h))
psi = np.multiply(np.arcsin(sinpsi),180) / np.pi
solarmax = np.multiply(SC,sinpsi) / gamma2
#solarmax=1380*sinpsi*(0.61+0.20*sinpsi);
T = np.minimum(2,sw_dn / solarmax)
Ts = np.arange(0,1+0.05,0.05)
As = np.arange(0,90+2,2)

a = np.array([[0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.062,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061],[0.072,0.07,0.068,0.065,0.065,0.063,0.062,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.061,0.06,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.087,0.083,0.079,0.073,0.07,0.068,0.066,0.065,0.064,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06],[0.115,0.108,0.098,0.086,0.082,0.077,0.072,0.071,0.067,0.067,0.065,0.063,0.062,0.061,0.061,0.06,0.06,0.06,0.06,0.061,0.061,0.061,0.061,0.06,0.059,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059],[0.163,0.145,0.13,0.11,0.101,0.092,0.084,0.079,0.072,0.072,0.068,0.067,0.064,0.063,0.062,0.061,0.061,0.061,0.06,0.06,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058],[0.235,0.198,0.174,0.15,0.131,0.114,0.103,0.094,0.083,0.08,0.074,0.074,0.07,0.067,0.065,0.064,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058],[0.318,0.263,0.228,0.192,0.168,0.143,0.127,0.113,0.099,0.092,0.084,0.082,0.076,0.072,0.07,0.067,0.065,0.064,0.062,0.062,0.06,0.06,0.06,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.058,0.057,0.058,0.058,0.058,0.058,0.057,0.057,0.057],[0.395,0.336,0.29,0.248,0.208,0.176,0.151,0.134,0.117,0.107,0.097,0.091,0.085,0.079,0.075,0.071,0.068,0.067,0.065,0.063,0.062,0.061,0.06,0.06,0.06,0.059,0.059,0.058,0.058,0.058,0.057,0.057,0.057,0.057,0.057,0.057,0.057,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.056,0.055],[0.472,0.415,0.357,0.306,0.252,0.21,0.176,0.154,0.135,0.125,0.111,0.102,0.094,0.086,0.081,0.076,0.072,0.071,0.068,0.066,0.065,0.063,0.062,0.061,0.06,0.059,0.058,0.057,0.057,0.057,0.056,0.055,0.055,0.055,0.055,0.055,0.055,0.054,0.053,0.054,0.053,0.053,0.054,0.054,0.053,0.053],[0.542,0.487,0.424,0.36,0.295,0.242,0.198,0.173,0.15,0.136,0.121,0.11,0.101,0.093,0.086,0.081,0.076,0.073,0.069,0.067,0.065,0.064,0.062,0.06,0.059,0.058,0.057,0.056,0.055,0.055,0.054,0.053,0.053,0.052,0.052,0.052,0.051,0.051,0.05,0.05,0.05,0.05,0.051,0.05,0.05,0.05],[0.604,0.547,0.498,0.407,0.331,0.272,0.219,0.185,0.16,0.141,0.127,0.116,0.105,0.097,0.089,0.083,0.077,0.074,0.069,0.066,0.063,0.061,0.059,0.057,0.056,0.055,0.054,0.053,0.053,0.052,0.051,0.05,0.05,0.049,0.049,0.049,0.048,0.047,0.047,0.047,0.046,0.046,0.047,0.047,0.046,0.046],[0.655,0.595,0.556,0.444,0.358,0.288,0.236,0.19,0.164,0.145,0.13,0.119,0.107,0.098,0.09,0.084,0.076,0.073,0.068,0.064,0.06,0.058,0.056,0.054,0.053,0.051,0.05,0.049,0.048,0.048,0.047,0.046,0.046,0.045,0.045,0.045,0.044,0.043,0.043,0.043,0.042,0.042,0.043,0.042,0.042,0.042],[0.693,0.631,0.588,0.469,0.375,0.296,0.245,0.193,0.165,0.145,0.131,0.118,0.106,0.097,0.088,0.081,0.074,0.069,0.065,0.061,0.057,0.055,0.052,0.05,0.049,0.047,0.046,0.046,0.044,0.044,0.043,0.042,0.042,0.041,0.041,0.04,0.04,0.039,0.039,0.039,0.038,0.038,0.038,0.038,0.038,0.038],[0.719,0.656,0.603,0.48,0.385,0.3,0.25,0.193,0.164,0.145,0.131,0.116,0.103,0.092,0.084,0.076,0.071,0.065,0.061,0.057,0.054,0.051,0.049,0.047,0.045,0.043,0.043,0.042,0.041,0.04,0.039,0.039,0.038,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.034,0.034,0.034,0.034,0.034,0.034],[0.732,0.67,0.592,0.474,0.377,0.291,0.246,0.19,0.162,0.144,0.13,0.114,0.1,0.088,0.08,0.072,0.067,0.062,0.058,0.054,0.05,0.047,0.045,0.043,0.041,0.039,0.039,0.038,0.037,0.036,0.036,0.035,0.035,0.034,0.033,0.032,0.032,0.032,0.031,0.031,0.031,0.03,0.03,0.03,0.03,0.03],[0.73,0.652,0.556,0.444,0.356,0.273,0.235,0.188,0.16,0.143,0.129,0.113,0.097,0.086,0.077,0.069,0.064,0.06,0.055,0.051,0.047,0.044,0.042,0.039,0.037,0.035,0.035,0.035,0.034,0.033,0.033,0.032,0.032,0.032,0.029,0.029,0.029,0.029,0.028,0.028,0.028,0.028,0.027,0.027,0.028,0.028],[0.681,0.602,0.488,0.386,0.32,0.252,0.222,0.185,0.159,0.142,0.127,0.111,0.096,0.084,0.075,0.067,0.062,0.058,0.054,0.05,0.046,0.042,0.04,0.036,0.035,0.033,0.032,0.032,0.031,0.03,0.03,0.03,0.03,0.029,0.027,0.027,0.027,0.027,0.026,0.026,0.026,0.026,0.026,0.026,0.026,0.026],[0.581,0.494,0.393,0.333,0.288,0.237,0.211,0.182,0.158,0.141,0.126,0.11,0.095,0.083,0.074,0.066,0.061,0.057,0.053,0.049,0.045,0.041,0.039,0.034,0.033,0.032,0.031,0.03,0.029,0.028,0.028,0.028,0.028,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.453,0.398,0.342,0.301,0.266,0.226,0.205,0.18,0.157,0.14,0.125,0.109,0.095,0.083,0.074,0.065,0.061,0.057,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.027,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025],[0.425,0.37,0.325,0.29,0.255,0.22,0.2,0.178,0.157,0.14,0.122,0.108,0.095,0.083,0.074,0.065,0.061,0.056,0.052,0.048,0.044,0.04,0.038,0.033,0.032,0.031,0.03,0.029,0.028,0.027,0.026,0.026,0.026,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025]])

if T.size==1:   ### for single value function
    Tchk = np.abs(Ts - T)
    i = np.array(np.where(Tchk == Tchk.min()))
    if psi < 0:
        alb = np.array([0])
        solarmax = np.array([0])
        T = np.array([0])
        j = np.array([0])
        psi = np.array([0])
    else:
        Achk = np.abs(As - psi)
        j = np.array(np.where(Achk == Achk.min()))
        szj = j.shape
        if szj[0] > 0:
            alb = a[i,j].flatten()
        else:
            #       print('no j found, not assigning alb to anything');
            pass
else:  ### for vectorized function
    for k in np.arange(0,np.size(sinpsi)).reshape(-1):
        Tchk = np.abs(Ts - T[k])
        i = np.array(np.where(Tchk == Tchk.min()))
        if psi[k] < 0:
            alb[k] = 0
            solarmax[k] = 0
            T[k] = 0
            j = 0
            psi[k] = 0
        else:
            Achk = np.abs(As - psi[k])
            j = np.array(np.where(Achk == Achk.min()))
            szj = j.shape
            if szj[0] > 0:
                alb[k] = a[i,j]
            else:
                #       disp('no j found, not assigning alb to anything');
                pass

#disp([num2str(jd) '  ' num2str(sw_dn) '  ' num2str(alb) '  ' num2str(T) '  ' num2str(i) '  ' num2str(j)])
#return alb,T,solarmax,psi

###########################################################
sw_net = np.multiply((1 - alb),sw_dn)


'''
timeS = times[oktimeg]
timestampS = timestamps[oktimeg]

tauS = A[:,1]
shtfluxS = A[:,2]
lhtfluxS = A[:,3]
CdS = A[:,12]
ChS = A[:,13]
CeS = A[:,14]

fig,ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(timeS,shtfluxS,'.-',markersize=1,color='lightcoral',label='Sensible heat flux')
ax1.plot(timeS,lhtfluxS,'.-',markersize=1,color='lightseagreen',label='latent heat flux')
ax1.legend(loc='center right')
plt.ylabel('($W/m^2$)',fontsize=14)
xfmt = mdates.DateFormatter('%b-%d')
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Enthalpy Fluxes using COARE3.6. Saildrone '+dataset_id)
ax2 = ax1.twinx()
ax2.plot(timeS,P,'.-',markersize=1,label='Sea Level Pressure')
plt.ylabel('(hPa)',fontsize=14)
ax2.legend()
'''

##############################################################################
# Read adcp saildrone data
'''
url_adcp = url_saildrone_adcp
gdata_adcp = xr.open_dataset(url_adcp)#,decode_times=False)

latitude_adcp = np.asarray(gdata_adcp.latitude)
longitude_adcp = np.asarray(gdata_adcp.longitude)
time_adcp = np.asarray(gdata_adcp.time)
depth_adcp = np.asarray(gdata_adcp.depth)
dataset_id_adcp = gdata_adcp.drone_id
vel_east = np.asarray(gdata_adcp.vel_east)
vel_north = np.asarray(gdata_adcp.vel_north)
'''
#########################################################################
'''
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,CdS,'.-',color='blue',label='Saildrone sd'+dataset_id)
#for i in np.arange(len(exp_names)):
#    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Cd '+ cycle,fontsize=18)
#plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

#########################################################################
fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,ChS,'.-',color='blue',label='Saildrone sd'+dataset_id)
#for i in np.arange(len(exp_names)):
#    plt.plot(target_timeSfv3[i],target_lhtflux[i,0:len(target_timeSfv3[i])],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.legend()
plt.legend(loc='lower right')
plt.title('Ch '+ cycle,fontsize=18)
#plt.ylabel('($W/m^2$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)
'''
