from netCDF4 import Dataset
import numpy as np
from pyproj import Geod
import math

#-----------------------------------------------
# coordinate rotation                           
#-----------------------------------------------
def rotate(u,v,phi):
    # phi is in radians
    u2 =  u*math.cos(phi) + v*math.sin(phi)
    v2 = -u*math.sin(phi) + v*math.cos(phi)
    return u2,v2

#-----------------------------------------------
# full cross-section transport calculation      
#-----------------------------------------------
def calc_transport(dates,fcst):
    """
    Calculate the transport of water across the Florida Straits
    This extracts the section and integrates the flow through it.
    """
    transport=[]
    if fcst==0:
        fcst_str='n024'
    else:
        fcst_str='f{:03d}'.format(fcst)
    cable_loc=np.loadtxt(refDir+'/eightmilecable.dat',dtype='int',usecols=(0,1))
    eightmile_lat = 26.5167
    eightmile_lon = -78.7833%360
    wpb_lat = 26.7153425
    wpb_lon = -80.0533746%360
    cable_angle = math.atan((eightmile_lat-wpb_lat)/(eightmile_lon-wpb_lon))
    g=Geod(ellps='WGS84')
    
    for date in dates:
        print('processing',date.strftime('%Y%m%d'),'fcst',fcst)
        rundate=date-timedelta(fcst/24.)  # calc rundate from fcst and date
        ufile=archDir+'/'+rundate.strftime('%Y%m%d')+'/rtofs_glo_3dz_'+fcst_str+'_daily_3zuio.nc'
        vfile=archDir+'/'+rundate.strftime('%Y%m%d')+'/rtofs_glo_3dz_'+fcst_str+'_daily_3zvio.nc'
        
        try:
            udata=Dataset(ufile)
            vdata=Dataset(vfile)
        except:
            print(rundate,fcst,'not found -- continuing')
            transport.append(np.nan)
            continue
        
        lon=udata['Longitude'][:]
        lat=udata['Latitude'][:]
        depth=udata['Depth'][:]
                                        
        #ssh.set_fill_value(value=0.0)

        usection=np.zeros((depth.shape[0],cable_loc.shape[0]))
        vsection=np.zeros((depth.shape[0],cable_loc.shape[0]))
        
        #rows=np.array(cable_loc[:,0], dtype=np.intp)
        #cols=np.array(cable_loc[:,1], dtype=np.intp)
        #levs=np.arange(33, dtype=np.intp)
        
        #x[np.ix_(rows, columns)]
        
        udata=udata['u'][:].squeeze()
        vdata=vdata['v'][:].squeeze()
        
        for ncol,(row,col) in enumerate(cable_loc):
            usection[:,ncol]=udata[:,row,col].filled(fill_value=0.0)
            vsection[:,ncol]=vdata[:,row,col].filled(fill_value=0.0)
            
        lon=lon[cable_loc[:,0],cable_loc[:,1]]
        lat=lat[cable_loc[:,0],cable_loc[:,1]]

        # compute the distances along the track
        _,_,dist=g.inv(lon[0:-1],lat[0:-1],lon[1:],lat[1:])
        depth=np.diff(depth)
        usection=usection[:-1,:-1]
        vsection=vsection[:-1,:-1]

        dist,depth=np.meshgrid(dist,depth)        
        u,v=rotate(usection,vsection,cable_angle)        
        trans1=(v*dist*depth).sum()/1e6        
        #print(date.strftime('%Y-%m-%d'),' transport:',transport,'Sv')
        transport.append(trans1)
    
    return transport

