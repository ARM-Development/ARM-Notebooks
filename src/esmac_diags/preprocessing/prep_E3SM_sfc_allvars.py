"""
# prepare E3SM surface aerosol properties at ARM sites
# input data is E3SM regional output
# output is surface variables at the nearest grid
"""

import os
import numpy as np
from ..subroutines.time_format_change import timeunit2cday, yyyymmdd2cday, cday2mmdd
from ..subroutines.read_netcdf import read_E3SM
from netCDF4 import Dataset

def find_nearest(xall, yall, x, y):
    distance = np.square(xall - x) + np.square(yall - y)
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
    E3SM_sfc_path = settings['E3SM_sfc_path']

    #%% other settings
        
    if not os.path.exists(E3SM_sfc_path):
        os.makedirs(E3SM_sfc_path)
    
    # change start date into calendar day
    cday1 = yyyymmdd2cday(start_date, 'noleap')
    cday2 = yyyymmdd2cday(end_date, 'noleap')
    if start_date[0:4] != end_date[0:4]:
        raise ValueError('currently not support multiple years. please set start_date and end_date in the same year')
    year0 = start_date[0:4]
        
    if site == 'SGP':
        E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    elif site == 'ENA':
        E3SMdomain_range = '330e_to_335e_37n_to_42n'   
    else:
        raise ValueError('data for this site is not specified: ' + site)
        
        
    #%% set variables
    for mm in range(len(Model_List)):
        model = Model_List[mm]
        variable1d_names = ['PS', 'PBLH', 'FLNT', 'FSNT', 'FLNS', 'FSNS', "LHFLX", "SHFLX",
                            'TREFHT','PRECT','PRECL', "TGCLDLWP", "TGCLDIWP"]
        variable2d_names = ['T', 'U', 'V', 'Q', 'RELHUM', 'RHW', 'RHI', 'CLOUD',
                            'CLDLIQ', 'CLDICE', 'NUMLIQ', 'NUMICE', 'CCN1', 'CCN3', 'CCN5', "AREI", "AREL", 
                            'bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', 'mom_a2', 'mom_a3', 'mom_a4',
                            'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
                            'soa_a1', 'soa_a2', 'soa_a3', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
                            'num_c1', 'num_c2', 'num_c3', 'num_c4', "dgnd_a01", "dgnd_a02", "dgnd_a03", "dgnd_a04",
                            "dgnw_a01", "dgnw_a02", "dgnw_a03", "dgnw_a04", 'EXTINCT', 'ABSORB']
        if model == 'NucSoaCond': # with so4 and soa in nucleation mode
            variable2d_names = variable2d_names + ['so4_a5', 'soa_a5', 'num_a5', 'num_c5', "dgnd_a05", "dgnw_a05"]
        elif model == 'Nuc':      # only with so4 in nucleation mode
            variable2d_names = variable2d_names + ['so4_a5','num_a5','num_c5', "dgnd_a05", "dgnw_a05"]
        var1dlen = len(variable1d_names)
        var2dlen = len(variable2d_names)
        varlen = var1dlen + var2dlen
        variable_names = variable1d_names + variable2d_names
        
        #%%  process data for each day
        for cday in range(cday1, cday2 + 1):
            mmdd = cday2mmdd(cday)
            date = year0 + '-' + mmdd[0:2] + '-' + mmdd[2:4]
            
            filename_input = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.' + date + '-00000.nc'
            
            # read in E3SM data
            variables = list()
            var_units = list()
            var_longnames = list()
            
            (timem, lonm, timeunitm, lonmunit, lonmname) = read_E3SM(filename_input, 'lon_' + E3SMdomain_range)
            (timem, latm, timeunitm, latmunit, latmname) = read_E3SM(filename_input, 'lat_' + E3SMdomain_range)
            x_idx = find_nearest(lonm, latm, lon0, lat0)
            
            
            (timem, var1d, timeunitm, var1dunit, var1dlongname) = \
                         read_E3SM(filename_input, [a + '_' + E3SMdomain_range for a in variable1d_names])
            (timem, var2d, timeunitm, var2dunit, var2dlongname) = \
                         read_E3SM(filename_input, [a + '_' + E3SMdomain_range for a in variable2d_names])
            for vv in range(var1dlen):
                variables.append(var1d[vv][:, x_idx])
            for vv in range(var2dlen):
                variables.append(var2d[vv][:, -1, x_idx])   # choose the lowest level
            var_units = var1dunit + var2dunit
            var_longnames = var1dlongname + var2dlongname
        
            cdaym = timeunit2cday(timeunitm, 'noleap')
            yearm = timeunitm.split(' ')[2][0:4]
            time = timem.data - 365*(int(year0) - int(yearm))  +  cdaym
            
            
            # %% output extacted file
            outputname = 'SFC_vars_' + campaign + '_' + model + '_' + date + '.nc'
            print('output to this file: ' + E3SM_sfc_path + outputname)
            
            # define filename
            f = Dataset(E3SM_sfc_path + outputname, 'w', format = 'NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time", ))
            lat_o = f.createVariable("lat", "f8", ())
            lon_o = f.createVariable("lon", "f8", ())
            var_o = list()
            for vv in range(varlen):
                var_o.append(f.createVariable(variable_names[vv], 'f8', ("time", )))
            
            # write data
            time_o[:] = time
            lat_o[:] = latm[x_idx]
            lon_o[:] = lonm[x_idx]
            for vv in range(varlen):
                var_o[vv][:] = np.array(variables[vv])
            
            # attributes
            time_o.units = "days since " + str(int(year0) - 1) + "-12-31 00:00:00 UTC"
            lat_o.units = "latitude"
            lon_o.units = "longitude"
            for vv in range(varlen):
                var_o[vv].units = var_units[vv]
                var_o[vv].long_name = var_longnames[vv]
            
            # global attributes
            import time as ttt
            f.description = model + " extact surface variables for " + campaign
            f.modeldata = filename_input
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
            
