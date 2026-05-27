#!/usr/bin/env python3

"""This script is to create a grib2 file with the 6h accumulated precipitation from the 3h accumulated precipitation."""

import os

import yaml
import numpy as np
import datetime

import grib2io

# Parse the yaml config file
print('Parse the config file: plot_atmos.yml:')
with open('plot_atmos.yml', 'rt') as f:
    conf = yaml.safe_load(f)

ff3 = np.arange(6,121,3)

for f in ff3:
    fm3 = f-3
    if f < 10:
        fhour = 'f00'+str(f)
    if f >= 10 and f < 100:
        fhour = 'f0'+str(f)
    if f > 100:
        fhour = 'f'+str(f)
    if fm3 < 10:
        fhourm3 = 'f00'+str(fm3)
    if fm3 >= 10 and fm3 < 100:
        fhourm3 = 'f0'+str(fm3)
    if fm3 > 100:
        fhourm3 = 'f'+str(fm3)

    print(fhour)
    print(fhourm3)

    fname = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+fhour+'.grb2'
    grib2file = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fname)
    print(f'grib2file: {grib2file}')
    grb = grib2io.open(grib2file,mode='r')
    print('Extracting the 3h accumulated precipitation')
    levstr='surface'
    apcp_msg = grb.select(shortName='APCP', level=levstr)[0]
    apcp_data = grb.select(shortName='APCP', level=levstr)[0].data

    fnamem3 = conf['stormModel'].lower()+'.'+conf['ymdh']+'.'+fhourm3+'.grb2'
    grib2filem3 = os.path.join(conf['COMarafs']+conf['ymdh']+'/00E/', fnamem3)
    print(f'grib2file: {grib2filem3}')
    grbm3 = grib2io.open(grib2filem3,mode='r')
    print('Extracting the 3h accumulated precipitation')
    levstr='surface'
    apcpm3_data = grbm3.select(shortName='APCP', level=levstr)[0].data

    apcp6h = apcp_msg
    apcp6h.timeRangeOfStatisticalProcess = 6
    apcp6h.duration = datetime.timedelta(seconds=21600)
    apcp6h.leadTime = datetime.timedelta(seconds=int(f*3600))
    apcp6h.data = apcpm3_data + apcp_data

    with grib2io.open(grib2file, mode='a') as fgrib2:
        print('Adding 6h acc to '+grib2file)
        fgrib2.write(apcp6h)



