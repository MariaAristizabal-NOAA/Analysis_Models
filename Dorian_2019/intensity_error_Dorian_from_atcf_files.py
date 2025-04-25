
#%% User input
scratch_dir = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/'
figs_dir = '/home/Maria.Aristizabal/Dorian_2019/Figures/'
storm_name = 'Dorian'
storm_id = 'dorian05l'
model1 = 'HWRF2019_POM'
mtype1 = 'oper'
model2 = 'HWRF2020_POM'
mtype2 = 'exp'
model3 = 'HWRF2020_HYCOM'
mtype3 = 'exp'
atcf_suffix = '.trak.hwrf.atcfunix'

cycles = ['2019082800','2019082806','2019082812','2019082818','2019082900','2019082906','2019082912','2019082918','2019083000','2019083006','2019083012','2019083018','2019083100','2019083106','2019083112','2019083118','2019090100','2019090106']
#cycles = ['2019082800','2019082806','2019082812']
#cycles = ['2019083000','2019083006','2019083012']
#cycles = ['2019083118']

RI_time = [2019,9,1,6] #year, month, day, hour

#%%
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
#import netCDF4
#from mpl_toolkits.basemap import Basemap
import cmocean
import glob
import os
from datetime import datetime,timedelta
import matplotlib.dates as mdates

#plt.switch_backend('agg')

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

Int_err_model1 = np.empty((len(cycles),22))
Int_err_model1[:] = np.nan
Int_err_model2 = np.empty((len(cycles),22))
Int_err_model2[:] = np.nan
Int_err_model3 = np.empty((len(cycles),22))
Int_err_model3[:] = np.nan
Int_err_model1_noRI = np.empty((len(cycles),22))
Int_err_model1_noRI[:] = np.nan
Int_err_model2_noRI = np.empty((len(cycles),22))
Int_err_model2_noRI[:] = np.nan
Int_err_model3_noRI = np.empty((len(cycles),22))
Int_err_model3_noRI[:] = np.nan
for c,cycle in enumerate(cycles):
    print(cycle)

    year = int(cycle[0:4])
    month = int(cycle[4:6])
    day = int(cycle[6:8])
    hour = int(cycle[8:10]) 
    t00 = datetime(year,month,day,hour)
    deltat= timedelta(hours=3) # every 6 hours
    time_hwrf1 = [t00+nstep*deltat for nstep in np.arange(43)]
    time_hwrf1 = np.asarray(time_hwrf1)
    time_hwrf2 = np.asarray(time_hwrf1)
    time_hwrf3 = np.asarray(time_hwrf1)

    Dir_model1 = scratch_dir + model1 + '_' + storm_name + '/' +model1 + '_' + storm_id + '.' + cycle + '_' + mtype1 + '/'
    Dir_model2 = scratch_dir + model2 + '_' + storm_name + '/' +model2 + '_' + storm_id + '.' + cycle + '_' + mtype2 + '/'
    Dir_model3 = scratch_dir + model3 + '_' + storm_name + '/' +model3 + '_' + storm_id + '.' + cycle + '_' + mtype3 + '/'

#%% Get list HWRF files
    HWRF1 = sorted(glob.glob(os.path.join(Dir_model1,'*.trak.hwrf.atcfunix')))
    HWRF2 = sorted(glob.glob(os.path.join(Dir_model2,'*.trak.hwrf.atcfunix')))
    HWRF3 = sorted(glob.glob(os.path.join(Dir_model3,'*.trak.hwrf.atcfunix')))

    if len(HWRF1) != 0:
        ff1 = open(HWRF1[0],'r')
        f1 = ff1.readlines()
        intt1 = []
        lead_time1 = []
        for l in f1:
            intt1.append(float(l.split(',')[8]))
            lead_time1.append(int(l.split(',')[5][1:4]))
        intt1 = np.asarray(intt1)
        lead_time_track1, ind = np.unique(lead_time1,return_index=True)
        int_track1 = intt1[ind]
        max_wind_10m_model1 = int_track1
    else:
        lead_time_track1 = np.empty(len(time_hwrf1))
        lead_time_track1[:] = np.nan
        max_wind_10m_model1 = np.empty(len(time_hwrf1))
        max_wind_10m_model1[:] = np.nan

    if len(HWRF2) != 0:
        ff2 = open(HWRF2[0],'r')
        f2 = ff2.readlines()
        intt2 = []
        lead_time2 = []
        for l in f2:
            intt2.append(float(l.split(',')[8]))
            lead_time2.append(int(l.split(',')[5][1:4]))
        intt2 = np.asarray(intt2)
        lead_time_track2, ind = np.unique(lead_time2,return_index=True)
        int_track2 = intt2[ind]
        max_wind_10m_model2 = int_track2
    else:
        lead_time_track2 = np.empty(len(time_hwrf2))
        lead_time_track2[:] = np.nan
        max_wind_10m_model2 = np.empty(len(time_hwrf2))
        max_wind_10m_model2[:] = np.nan

    if len(HWRF3) != 0:
        ff3 = open(HWRF3[0],'r')
        f3 = ff3.readlines()
        intt3 = []
        lead_time3= []
        for l in f3:
            intt3.append(float(l.split(',')[8]))
            lead_time3.append(int(l.split(',')[5][1:4]))
        intt3 = np.asarray(intt3)
        lead_time_track3, ind = np.unique(lead_time3,return_index=True)
        int_track3 = intt3[ind]
        max_wind_10m_model3 = int_track3
    else:
        lead_time_track3 = np.empty(len(time_hwrf3))
        lead_time_track3[:] = np.nan
        max_wind_10m_model3 = np.empty(len(time_hwrf3))
        max_wind_10m_model3[:] = np.nan

    year = int(cycle[0:4])
    month = int(cycle[4:6])
    day = int(cycle[6:8])
    hour = int(cycle[8:10]) 
    t00 = datetime(year,month,day,hour)
    deltat= timedelta(hours=3) # every 6 hours
    time_hwrf1 = [t00+nstep*deltat for nstep in np.arange(len(lead_time_track1))]
    time_hwrf1 = np.asarray(time_hwrf1)
    time_hwrf2 = np.asarray(time_hwrf1)
    time_hwrf3 = np.asarray(time_hwrf1)

    okt1 = np.logical_and(mdates.date2num(time_best_track) >= mdates.date2num(time_hwrf1[0]),mdates.date2num(time_best_track) <= mdates.date2num(time_hwrf1[-1]))
    okt2 = np.logical_and(mdates.date2num(time_best_track) >= mdates.date2num(time_hwrf2[0]),mdates.date2num(time_best_track) <= mdates.date2num(time_hwrf2[-1]))
    okt3 = np.logical_and(mdates.date2num(time_best_track) >= mdates.date2num(time_hwrf3[0]),mdates.date2num(time_best_track) <= mdates.date2num(time_hwrf3[-1]))

    t01 = time_best_track[okt1][0] 
    t02 = time_best_track[okt2][0] 
    t03 = time_best_track[okt3][0] 
    lead_time1 = np.asarray([np.int((t-t01).total_seconds()/3600) for t in time_best_track[okt1]])
    lead_time2 = np.asarray([np.int((t-t02).total_seconds()/3600) for t in time_best_track[okt2]])
    lead_time3 = np.asarray([np.int((t-t03).total_seconds()/3600) for t in time_best_track[okt3]])

    int_err_model1 = (wind_int_kt[okt1] - max_wind_10m_model1[::2])
    int_err_model2 = (wind_int_kt[okt2] - max_wind_10m_model2[::2])
    int_err_model3 = (wind_int_kt[okt3] - max_wind_10m_model3[::2])

    lead_time = np.arange(0,132,6)

    oki1 = [np.where(l == lead_time)[0][0] for l in lead_time1]
    oki2 = [np.where(l == lead_time)[0][0] for l in lead_time2]
    oki3 = [np.where(l == lead_time)[0][0] for l in lead_time3]

    Int_err_model1[c,oki1] = int_err_model1
    Int_err_model2[c,oki2] = int_err_model2
    Int_err_model3[c,oki3] = int_err_model3

    # Rapid intensification
    no_RI1 = time_best_track[okt1] <= datetime(RI_time[0],RI_time[1],RI_time[2],RI_time[3])
    no_RI2 = time_best_track[okt2] <= datetime(RI_time[0],RI_time[1],RI_time[2],RI_time[3])
    no_RI3 = time_best_track[okt3] <= datetime(RI_time[0],RI_time[1],RI_time[2],RI_time[3])
    okn1 = [np.where(l == lead_time)[0][0] for l in lead_time1[no_RI1]]
    okn2 = [np.where(l == lead_time)[0][0] for l in lead_time2[no_RI2]]
    okn3 = [np.where(l == lead_time)[0][0] for l in lead_time3[no_RI3]]
    lead_time_RI = lead_time1[no_RI1][-1]

    Int_err_model1_noRI[c,okn1] = int_err_model1[no_RI1]
    Int_err_model2_noRI[c,okn2] = int_err_model2[no_RI2]
    Int_err_model3_noRI[c,okn3] = int_err_model3[no_RI3]

    '''
    #%% Figure intensity
    fig,ax1 = plt.subplots(figsize=(6, 4))
    plt.ion()
    plt.plot(time_best_track,wind_int_kt,'o-k',label='Best')
    plt.plot(time_hwrf,max_wind_10m_model1,'X-',color='mediumorchid',label='POM Oper',markeredgecolor='k',markersize=7)
    plt.plot(time_hwrf,max_wind_10m_model2,'^-',color='teal',label='POM Exp',markeredgecolor='k',markersize=7)
    plt.plot(time_hwrf,max_wind_10m_model3,'H-',color='darkorange',label='HYCOM Exp',markeredgecolor='k',markersize=7)
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

    file_name = figs_dir + 'Dorian_intensity_cycle_'+cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 
    '''

    #%% Intensity error
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
   
    fig,ax1 = plt.subplots(figsize=(10, 5))
    plt.ion()
    plt.plot(lead_time1,int_err_model1,'X-',color='mediumorchid',label=model1+' (IC clim.)',markeredgecolor='k',markersize=7)
    plt.plot(lead_time2,int_err_model2,'^-',color='teal',label=model2+' (IC RTOFS)',markeredgecolor='k',markersize=7)
    plt.plot(lead_time3,int_err_model3,'H-',color='darkorange',label=model3+' (IC RTOFS)',markeredgecolor='k',markersize=7)
    plt.plot(lead_time1,np.tile(0,len(lead_time1)),'--k')
    plt.plot(np.tile(lead_time_RI,len(np.arange(50))),np.arange(50),'--k')
    plt.xlim([0,126])
    ax1.xaxis.set_major_locator(MultipleLocator(12))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax1.xaxis.set_minor_locator(MultipleLocator(3))
    plt.title('Intensity Forecast Error Dorian '+ cycle,fontsize=18)
    plt.ylabel('Forecast Error (Kt)',fontsize=14)
    plt.xlabel('Forecast Lead Time (Hr)',fontsize=14)
    plt.legend(loc='upper left',fontsize=14)
    file_name = figs_dir + 'Dorian_intensity_error_cycle_'+cycle
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 

#%%
Int_err_model1_noRI_mean = np.nanmean(Int_err_model1_noRI,0)
Int_err_model1_noRI_min = np.nanmin(Int_err_model1_noRI,0)
Int_err_model1_noRI_max = np.nanmax(Int_err_model1_noRI,0)

Int_err_model2_noRI_mean = np.nanmean(Int_err_model2_noRI,0)
Int_err_model2_noRI_min = np.nanmin(Int_err_model2_noRI,0)
Int_err_model2_noRI_max = np.nanmax(Int_err_model2_noRI,0)

Int_err_model3_noRI_mean = np.nanmean(Int_err_model3_noRI,0)
Int_err_model3_noRI_min = np.nanmin(Int_err_model3_noRI,0)
Int_err_model3_noRI_max = np.nanmax(Int_err_model3_noRI,0)

#%%
lead_time = np.arange(0,132,6)

fig,ax = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time,Int_err_model1_noRI_mean,'X-',color='mediumorchid',label=model1+' (IC clim.)',markeredgecolor='k',markersize=7)
plt.plot(lead_time,Int_err_model2_noRI_mean,'^-',color='teal',label=model2+' (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time,Int_err_model3_noRI_mean,'H-',color='darkorange',label=model3+' (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time,np.tile(0,len(lead_time)),'--k')
plt.plot(np.tile(lead_time_RI,len(np.arange(50))),np.arange(50),'--k')
plt.xlim([0,126])
#plt.xlim([0,84])
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
plt.title('Intensity Forecast Mean Error Dorian (no RI) '+ cycles[0] +'-'+cycles[-1],fontsize=18)
plt.ylabel('Forecast Error (Kt)',fontsize=16)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=16)
plt.legend(loc='upper left',fontsize=14)
file_name = figs_dir + 'Dorian_mean_intensity_error-'+cycle[0] + '-'+cycle[-1]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 

#%%
lead_time = np.arange(0,132,6)

fig,ax = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time,Int_err_model1_noRI_mean,'X-',color='mediumorchid',label=model1+' (IC clim.)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time,Int_err_model1_noRI_min,Int_err_model1_noRI_max,color='mediumorchid',alpha=0.1)
plt.plot(lead_time,Int_err_model1_noRI_min,'-',color='mediumorchid',alpha=0.5)
plt.plot(lead_time,Int_err_model1_noRI_max,'-',color='mediumorchid',alpha=0.5)

plt.plot(lead_time,Int_err_model2_noRI_mean,'^-',color='teal',label=model2+' (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time,Int_err_model2_noRI_min,Int_err_model2_noRI_max,color='teal',alpha=0.1)
plt.plot(lead_time,Int_err_model2_noRI_min,'-',color='teal',alpha=0.5)
plt.plot(lead_time,Int_err_model2_noRI_max,'-',color='teal',alpha=0.5)

plt.plot(lead_time,Int_err_model3_noRI_mean,'H-',color='darkorange',label=model3+' (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time,Int_err_model3_noRI_min,Int_err_model3_noRI_max,color='darkorange',alpha=0.1)
plt.plot(lead_time,Int_err_model3_noRI_min,'-',color='darkorange',alpha=0.5)
plt.plot(lead_time,Int_err_model3_noRI_max,'-',color='darkorange',alpha=0.5)

#plt.plot(np.tile(lead_time_RI,len(np.arange(50))),np.arange(50),'--k')
plt.plot(lead_time,np.tile(0,len(lead_time)),'--k')
plt.xlim([0,126])
#plt.xlim([0,84])
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
plt.title('Intensity Forecast Mean Error Dorian (no RI) '+ cycles[0] +'-'+cycles[-1],fontsize=18)
plt.ylabel('Forecast Error (Kt)',fontsize=16)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=16)
plt.legend(loc='lower right',fontsize=14)
file_name = figs_dir + 'Dorian_mean_intensity_error2-'+cycle[0] + '-'+cycle[-1]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 

#%%
lead_time = np.arange(0,132,6)

fig,ax = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time[0:15],Int_err_model1_noRI_mean[0:15],'X-',color='mediumorchid',label=model1+' (IC clim.)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time[0:15],Int_err_model1_noRI_min[0:15],Int_err_model1_noRI_max[0:15],color='mediumorchid',alpha=0.1)
plt.plot(lead_time[0:15],Int_err_model1_noRI_min[0:15],'-',color='mediumorchid',alpha=0.5)
plt.plot(lead_time[0:15],Int_err_model1_noRI_max[0:15],'-',color='mediumorchid',alpha=0.5)

plt.plot(lead_time[0:15],Int_err_model2_noRI_mean[0:15],'^-',color='teal',label=model2+' (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time[0:15],Int_err_model2_noRI_min[0:15],Int_err_model2_noRI_max[0:15],color='teal',alpha=0.1)
plt.plot(lead_time[0:15],Int_err_model2_noRI_min[0:15],'-',color='teal',alpha=0.5)
plt.plot(lead_time[0:15],Int_err_model2_noRI_max[0:15],'-',color='teal',alpha=0.5)

plt.plot(lead_time[0:15],Int_err_model3_noRI_mean[0:15],'H-',color='darkorange',label=model3+' (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time[0:15],Int_err_model3_noRI_min[0:15],Int_err_model3_noRI_max[0:15],color='darkorange',alpha=0.1)
plt.plot(lead_time[0:15],Int_err_model3_noRI_min[0:15],'-',color='darkorange',alpha=0.5)
plt.plot(lead_time[0:15],Int_err_model3_noRI_max[0:15],'-',color='darkorange',alpha=0.5)

#plt.plot(np.tile(lead_time_RI,len(np.arange(50))),np.arange(50),'--k')
plt.plot(lead_time[0:15],np.tile(0,len(lead_time[0:15])),'--k')
#plt.xlim([0,126])
plt.xlim([0,84])
ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
plt.title('Intensity Forecast Mean Error Dorian (no RI) '+ cycles[0] +'-'+cycles[-1],fontsize=18)
plt.ylabel('Forecast Error (Kt)',fontsize=16)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=16)
plt.legend(loc='upper left',fontsize=13)
file_name = figs_dir + 'Dorian_mean_intensity_error3_'+cycles[0] + '_'+cycles[-1]
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
