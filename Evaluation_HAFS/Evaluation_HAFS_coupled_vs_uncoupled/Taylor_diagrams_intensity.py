
Mettc_file = '/scratch1/NCEPDEV/hwrf/save/John.Steffen/data/metTC/tc_stat_lew_full2.tcst'

#===================================================================================================
def taylor_template(angle_lim,std_lim):

    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure(figsize=(7,7))
    tr = PolarAxes.PolarTransform()

    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))

    STDgrid=np.arange(0,std_lim,.5)
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))

    ra0, ra1 =0, angle_lim
    cz0, cz1 = 0, std_lim
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)

    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)

    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["left"].label.set_text("Normalized Standard Deviation")
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")

    ax1.axis["bottom"].set_visible(False)
    ax1 = ax1.get_aux_axes(tr)

    plt.grid(linestyle=':',alpha=0.5)

    return fig,ax1

#===================================================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#===================================================================================================
# Read in adeckFile as pandas dataFrame
print('Read in METtc file')
deck = pd.read_csv(Mettc_file,sep='\s+')

coupled = deck['AMODEL'] == 'HAFSc'
uncoupled = deck['AMODEL'] == 'HAFSu'

GoM_coupled = np.logical_and(deck['BLAT'][coupled] >= 20,deck['BLAT'][coupled] <= 30)


plt.figure()
plt.plot(deck['BLAT'][coupled],deck['ALAT'][coupled],'.b')
plt.plot(deck['BLAT'][uncoupled],deck['ALAT'][uncoupled],'.r')

plt.figure()
plt.plot(deck['BLON'][coupled],deck['ALON'][coupled],'.b')
plt.plot(deck['BLON'][uncoupled],deck['ALON'][uncoupled],'.r')

plt.figure()
plt.plot(deck['BMAX_WIND'][uncoupled],deck['AMAX_WIND'][uncoupled],'.r',label='uncupled')
plt.plot(deck['BMAX_WIND'][coupled],deck['AMAX_WIND'][coupled],'.b',label='coupled')
plt.plot(np.arange(150),np.arange(150),'--k')
plt.legend()
plt.title('Max Wind')
plt.xlabel('Bdeck')
plt.ylabel('Adeck')

corr_intensity_uncoup = np.corrcoef(deck['BMAX_WIND'][uncoupled],deck['AMAX_WIND'][uncoupled])[0,1]
corr_intensity_coup = np.corrcoef(deck['BMAX_WIND'][coupled],deck['AMAX_WIND'][coupled])[0,1]

std_intensity_bdeck_uncoup = np.std(deck['BMAX_WIND'][uncoupled])
std_intensity_bdeck_coup = np.std(deck['BMAX_WIND'][coupled])
std_intensity_adeck_uncoup = np.std(deck['AMAX_WIND'][uncoupled])
std_intensity_adeck_coup = np.std(deck['AMAX_WIND'][coupled])

#===================================================================================================
angle_lim = np.pi/2
std_lim = 1.5
theta_uncoup = np.arccos(corr_intensity_uncoup)
rr_uncoup = std_intensity_adeck_uncoup/std_intensity_bdeck_uncoup
theta_coup = np.arccos(corr_intensity_coup)
rr_coup = std_intensity_adeck_coup/std_intensity_bdeck_coup

fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^']
ax1.plot(0,1,'o',label='Bdeck',color='green',markersize=10,markeredgecolor='k')
ax1.plot(theta_uncoup,rr_uncoup,markers[0],color = 'r',markersize=8,markeredgecolor='k',label='Adeck uncoupled')
ax1.plot(theta_coup,rr_coup,markers[1],color = 'b',markersize=8,markeredgecolor='k',label='Adeck coupled')
plt.legend()


