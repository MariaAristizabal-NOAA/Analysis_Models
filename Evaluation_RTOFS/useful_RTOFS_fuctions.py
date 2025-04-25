##################################################################################
#%% Find the grid point positions for rtofs output given a longitude and latitude
def find_grid_position_hycom(lat,lon,target_lat,target_lon):

    # search in xi_rho direction
    oklatmm = []
    oklonmm = []

    for pos_xi in np.arange(lat.shape[1]):
        pos_eta = np.round(np.interp(target_lat,lat[:,pos_xi],np.arange(len(lat[:,pos_xi])),left=np.nan,right=np.nan))
        if np.isfinite(pos_eta):
            oklatmm.append((pos_eta).astype(int))
            oklonmm.append(pos_xi)

    pos = np.round(np.interp(target_lon,lon[oklatmm,oklonmm],np.arange(len(lon[oklatmm,oklonmm])))).astype(int)
    oklat = oklatmm[pos]
    oklon = oklonmm[pos]

    return oklat,oklon

##################################################################################
def get_glider_transect_from_RTOFS_ncfiles(ncfiles,lon,lat,depth,time_name,temp_name,salt_name,long,latg,tstamp_glider):

    target_temp = np.empty((len(depth),len(ncfiles)))
    target_temp[:] = np.nan
    target_salt = np.empty((len(depth),len(ncfiles)))
    target_salt[:] = np.nan
    target_depth = np.empty((len(depth),len(ncfiles)))
    target_depth[:] = np.nan
    target_time = []

    for x,file in enumerate(ncfiles):
        print(x)
        model = xr.open_dataset(file)
        t = model[time_name][:]
        tstamp_model = mdates.date2num(t)[0]
        target_time.append(mdates.num2date(tstamp_model))

        # Interpolating latg and longlider onto RTOFS grid
        sublon = np.interp(tstamp_model,tstamp_glider,long)
        sublat = np.interp(tstamp_model,tstamp_glider,latg)

        oklat, oklon = find_grid_position_hycom(lat,lon,sublat,sublon)

        target_temp[:,x] = np.asarray(model[temp_name][0,:,oklat,oklon])
        target_salt[:,x] = np.asarray(model[salt_name][0,:,oklat,oklon])

    return target_time, target_temp, target_salt

##################################################################################
def get_glider_transect_from_RTOFS_archv(afiles,lat,lon,nz,long,latg,tstamp_glider,jdm,idm):

    #target_t, target_temp_hafs_oc, target_salt_hafs_oc, depth = \
    #get_glider_transect_from_RTOFS_archv(afiles,lat,lon,nz,lon_glid,lat_glid,tstamp_glider,jdm,idm)
    #long = lon_glid
    #latg = lat_glid 
    
    layers = np.arange(0,nz)

    target_temp = np.empty((nz,len(afiles)))
    target_temp[:] = np.nan
    target_salt = np.empty((nz,len(afiles)))
    target_salt[:] = np.nan
    target_depth = np.empty((nz,len(afiles)))
    target_depth[:] = np.nan
    target_time = []

    for x,file in enumerate(afiles):
        print(file)
        lines = [line.rstrip() for line in open(file[:-1]+'b')]
        time_stamp = lines[-1].split()[2]
        hycom_days = lines[-1].split()[3]
        tzero = datetime(1901,1,1,0,0)
        timeRT = tzero+timedelta(float(hycom_days)-1)
        tstamp_model = mdates.date2num(timeRT)
        target_time.append(timeRT)

        # Interpolating latg and longlider onto RTOFS grid
        sublon = np.interp(tstamp_model,tstamp_glider,long)
        sublat = np.interp(tstamp_model,tstamp_glider,latg)

        oklat, oklon = find_grid_position_hycom(lat,lon,sublat,sublon)

        #ztmp = readVar(file[:-2],'archive','srfhgt',[0])*0.01 # converts [cm] to [m]
        #target_ztmp = ztmp[oklat,oklon]
        target_ztmp = reatrive_one_value_rtofs_archv(file[:-2],'srfhgt',str(0),jdm,idm,oklat,oklon)*0.01
        for lyr in tuple(layers):
            print(lyr)
            target_temp[lyr,x] = reatrive_one_value_rtofs_archv(file[:-2],'temp',str(lyr+1),jdm,idm,oklat,oklon)
            target_salt[lyr,x] = reatrive_one_value_rtofs_archv(file[:-2],'salin',str(lyr+1),jdm,idm,oklat,oklon)
            #dp = readVar(file[:-2],'archive','thknss',[lyr+1])/2/9806
            target_dp = reatrive_one_value_rtofs_archv(file[:-2],'thknss',str(lyr+1),jdm,idm,oklat,oklon)/9806
            target_ztmp = np.append(target_ztmp,target_dp)

        target_depth[:,x] = np.cumsum(target_ztmp[0:-1]) + np.diff(np.cumsum(target_ztmp))/2

    return target_time, target_temp, target_salt, target_depth

#####################################################################
def reatrive_one_value_rtofs_archv(rtofs_file,var_name,klayer,jdm,idm,oklat,oklon):

    lines = [line.rstrip() for line in open(rtofs_file+'.b')]
    ijdm = idm*jdm
    npad = 4096-(ijdm%4096)
    fld = ma.array([],fill_value=1.2676506002282294e+30)

    inFile = rtofs_file + '.a'

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == 'field':
                nheading = n + 1
                #print(line.split()[0])

    for n,line in enumerate(lines):
        if len(line) > 0:
            if line.split()[0] == var_name and line.split()[4] == klayer:
                nvar = n - nheading
                #print(nvar)
                #print(n)

    fid = open(inFile,'rb')
    fid.seek((nvar)*4*(npad+ijdm),0)
    fld = fid.read(ijdm*4)
    fld = struct.unpack('>'+str(ijdm)+'f',fld)
    #fld = np.array(fld)
    #fld = ma.reshape(fld,(jdm,idm))
    target_fld = fld[idm*(oklat)+oklon]

    #fid.seek((nvar)*4*(npad+idm*(oklat)+oklon),0)
    #fld = fid.read(4*(idm*(oklat)+oklon))
    #fld = struct.unpack('>'+str(idm*(oklat)+oklon)+'f',fld)
    #target_fld = fld[idm*oklat+oklon-50:idm*oklat+oklon+50]

    return target_fld

#####################################################################

