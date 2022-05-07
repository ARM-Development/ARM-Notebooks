"""
# prepare E3SM vertical profiles at ARM sites
# input data is E3SM regional output
# output is variables at the nearest column
"""

import os
import numpy as np
from ..subroutines.time_format_change import timeunit2cday, yyyymmdd2cday, cday2mmdd
from ..subroutines.read_netcdf import read_E3SM
from netCDF4 import Dataset

def find_nearest(xall, yall, x, y):
    distance = np.square(xall-x) + np.square(yall-y)
    idx = distance.argmin()
    return(idx)

def run_prep(settings):
    #%% variables from settings
    campaign = settings['campaign']
    lat0 = settings['lat0']
    lon0 = settings['lon0']
    site = settings['site']
    start_date = settings['start_date']
    end_date = settings['end_date']
    Model_List = settings['Model_List']
    E3SM_hourly_path = settings['E3SM_hourly_path']
    E3SM_hourly_filehead = settings['E3SM_hourly_filehead']
    E3SM_profile_path = settings['E3SM_profile_path']

    #%% other settings
    
    # output height above ground. data will be interpolated into z_f
    z_f = np.hstack((np.arange(0, 500, 50), np.arange(500, 2000, 100), np.arange(2000, 5000, 300),
                     np.arange(5000, 10000, 500), np.arange(10000, 20001, 1000)))
    zlen = len(z_f)
    
    if not os.path.exists(E3SM_profile_path):
        os.makedirs(E3SM_profile_path)
    
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date, 'noleap')
    cday2 = yyyymmdd2cday(end_date, 'noleap')
    if start_date[0:4]!=end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
    
    #%% set variables for profiles
    variable_names = ['T', 'U', 'V', 'Q', 'RELHUM', 'RHW', 'RHI', 'CLOUD', 'LWC', 'IWC',
                'CLDLIQ', 'CLDICE', 'NUMLIQ', "AREI",  "AREL"]
    varlen = len(variable_names)
        
    if site=='SGP':
        E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    elif site=='ENA':
        E3SMdomain_range = '330e_to_335e_37n_to_42n'   
    else:
        raise ValueError('data for this site is not specified: ' + site)
    
    for mm in range(len(Model_List)):
        model = Model_List[mm]
        
        #%%  process data
        for cday in range(cday1, cday2 + 1):
            mmdd = cday2mmdd(cday)
            date = year0 + '-' + mmdd[0:2] + '-' + mmdd[2:4]
            
            # read in E3SM data
            variables = list()
            var_units = list()
            var_longnames = list()
             
            filename_input = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.' + date + '-00000.nc'
        
            (timem, lonm, timeunitm, lonmunit, lonmname) = read_E3SM(filename_input, 'lon_' + E3SMdomain_range)
            (timem, latm, timeunitm, latmunit, latmname) = read_E3SM(filename_input, 'lat_' + E3SMdomain_range)
            (timem, z3, timeunitm, zunit, zname) = read_E3SM(filename_input, 'Z3_' + E3SMdomain_range)
            
            x_idx = find_nearest(lonm, latm, lon0, lat0)
            zm = z3[:, :, x_idx]
            
            # read in all variables
            (timem, var2d, timeunitm, var2dunit, var2dlongname) = \
                         read_E3SM(filename_input, [a + '_' + E3SMdomain_range for a in variable_names])
              
            tlen = len(timem)
            for vv in range(varlen):
                var = var2d[vv][:, :, x_idx]
                var2 = np.full((tlen, zlen), np.nan)
                for tt in range(tlen):
                    # interpolate height above sea level to height above ground
                    var2[tt, :] = np.interp(z_f, np.flip(zm[tt, :]-zm[tt, -1]), np.flip(var[tt, :]))
                variables.append(var2) 
                var_units.append(var2dunit[vv])
                var_longnames.append(var2dlongname[vv])      
        
            cdaym = timeunit2cday(timeunitm, 'noleap')
            yearm = timeunitm.split(' ')[2][0:4]
            time = timem.data - 365*(int(year0)-int(yearm)) + cdaym
            
        
            # %% output extacted file
            outputname = 'Profile_vars_' + campaign + '_' + model + '.' + date + '.nc'
            print('output to this file: ' + E3SM_profile_path + outputname)
            
            # define filename
            f = Dataset(E3SM_profile_path + outputname, 'w', format = 'NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            z = f.createDimension('height', zlen)
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time", ))
            height_o = f.createVariable("height", "f8", ("height", ))
            lat_o = f.createVariable("lat", "f8", ())
            lon_o = f.createVariable("lon", "f8", ())
            var_o = list()
            for vv in range(varlen):
                var_o.append (f.createVariable(variable_names[vv], 'f8', ("time", "height")))
            
            # write data
            time_o[:] = time
            height_o[:] = z_f
            lat_o[:] = latm[x_idx]
            lon_o[:] = lonm[x_idx]
            for vv in range(varlen):
                var_o[vv][:] = np.array(variables[vv])
            
            # attributes
            time_o.units = "days since " + str(int(year0)-1) + "-12-31 00:00:00 UTC"
            lat_o.units = "latitude"
            lon_o.units = "longitude"
            height_o.units = "gpm above ground"
            for vv in range(varlen):
                var_o[vv].units = var_units[vv]
                var_o[vv].long_name = var_longnames[vv]
            
            # global attributes
            import time as ttt
            f.description = model + " extact vertical variables for " + campaign
            f.modeldata = filename_input
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
