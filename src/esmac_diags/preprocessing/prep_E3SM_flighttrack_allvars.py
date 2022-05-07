"""
# prepare E3SM aerosol variables for flight tracks
# input data is IWG measurements from aircraft and E3SM regional output
# output is aerosol variables for each flight
"""

import glob
import os
import numpy as np
from ..subroutines.time_format_change import hhmmss2sec, timeunit2cday
from ..subroutines.read_aircraft import read_iwg1, read_RF_NCAR
from ..subroutines.read_netcdf import read_E3SM
from netCDF4 import Dataset

def find_nearest(xall, yall, x, y):
    distance = np.square(xall-x) + np.square(yall-y)
    idx = distance.argmin()
    return(idx)

def run_prep(settings):
    #%% variables from settings
    campaign = settings['campaign']
    Model_List = settings['Model_List']
    E3SM_hourly_path = settings['E3SM_hourly_path']
    E3SM_hourly_filehead = settings['E3SM_hourly_filehead']
    E3SM_aircraft_path = settings['E3SM_aircraft_path']

    if campaign in ['HISCALE', 'ACEENA']:
        IOP = settings['IOP']
        iwgpath = settings['iwgpath']
    elif campaign in ['CSET', 'SOCRATES']:
        RFpath = settings['RFpath']
    else:
        raise ValueError('this aircraft campaign is not recognized: ' + campaign)

    #%% other settings
    
    if not os.path.exists(E3SM_aircraft_path):
        os.makedirs(E3SM_aircraft_path)
    
    
    #%% find all flight data
    if campaign == 'HISCALE':
        lst = glob.glob(iwgpath + '*a2.txt')
        lst.sort()
        if IOP == 'IOP1':
            lst = lst[0:17]
        elif IOP == 'IOP2':
            lst = lst[17:]
        elif IOP[0:4] == '2016':
            a = lst[0].split('_' + campaign + '_')
            lst = glob.glob(a[0] + '*' + IOP + '*')
            lst.sort()
    elif campaign == 'ACEENA':
        lst = glob.glob(iwgpath + '*a2.txt')
        lst.sort()
        if IOP == 'IOP1':
            lst = lst[0:20]
        elif IOP == 'IOP2':
            lst = lst[20:]
        elif IOP[0:4] == '2017' or IOP[0:4] == '2018':
            a = lst[0].split('_' + campaign + '_')
            lst = glob.glob(a[0] + '*' + IOP + '*')
            lst.sort()
    elif campaign in ['CSET', 'SOCRATES']:
        lst = glob.glob(RFpath + 'RF*.PNI.nc')
        lst.sort()
    else:
        raise ValueError('this aircraft campaign is not recognized: ' + campaign)
        
    print('total number of files:' + str(len(lst)))
    
    for filename in lst:
        
        fname = filename.split('.')
        #%% read in flight data
        if campaign in ['HISCALE', 'ACEENA']:
            date = fname[-3]
            print('input data for ' + date)
            # year = date[0:4]
            # month = date[4:6]
            
            (flight, flightvars) = read_iwg1(filename)
            timelen = len(flight)
            # get lat, lon, height, time
            lon = np.empty(timelen)
            lat = np.empty(timelen)
            height = np.empty(timelen)
            time = np.empty(timelen)
            if np.logical_and(campaign == 'ACEENA', date == '20180216a'):
                flight.insert(1403, list(flight[1403]))
                tstr = flight[1403][1]
                tstr = tstr[0:-1] + str(int(tstr[-1])-1)
                flight[1403][1] = tstr
                del flight[-1]
            for t in range(timelen):
                lat[t] = float(flight[t][2])
                lon[t] = float(flight[t][3]) + 360
                height[t] = float(flight[t][4])
                timestr = flight[t][1].split(' ')
                time[t] = hhmmss2sec(timestr[1])
        
        elif campaign in ['CSET', 'SOCRATES']:
            date = fname[-4]
            print('input data for ' + date)
            # year = date[0:4]
            # month = date[4:6]
            (time, height, timeunit, hunit, hlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'ALT')
            (time, lat, timeunit, latunit, latlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'LAT')
            (time, lon, timeunit, lonunit, lonlongname, cellsize, cellunit) = read_RF_NCAR(filename, 'LON')
            lon[lon<0] = lon[lon<0] + 360
            
        #%% set variables and model region
        for mm in range(len(Model_List)):
            model = Model_List[mm]
            variable_names = ['T', 'U', 'V', 'Q', 'RELHUM', 'RHW', 'RHI', 'CLOUD', 'CLDLIQ', 'NUMLIQ',
                              'CLDICE', 'NUMICE', 'AWNI', 'AWNC','AREI', 'AREL', 'FICE', 'IWC', 'LWC', 
                              'CCN1', 'CCN3', 'CCN5', 'bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3',
                              'mom_a1', 'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2', 'ncl_a3',
                              'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
                              'soa_a1', 'soa_a2', 'soa_a3', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
                              'num_c1', 'num_c2', 'num_c3', 'num_c4', "dgnd_a01", "dgnd_a02", "dgnd_a03", "dgnd_a04",
                              "dgnw_a01", "dgnw_a02", "dgnw_a03", "dgnw_a04", 'EXTINCT', 'ABSORB']
        
            if model == 'NucSoaCond': # with so4 and soa in nucleation mode
                variable_names = variable_names + ['so4_a5', 'soa_a5', 'num_a5', 'num_c5', "dgnd_a05", "dgnw_a05"]
            elif model == 'Nuc':      # only with so4 in nucleation mode
                variable_names = variable_names + ['so4_a5', 'num_a5', 'num_c5', "dgnd_a05", "dgnw_a05"]
            varlen = len(variable_names)
            
            if campaign == 'HISCALE':
                E3SMdomain_range = '260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
            elif campaign == 'ACEENA':
                E3SMdomain_range = '330e_to_335e_37n_to_42n'   
            elif campaign == 'CSET':
                E3SMdomain_range = '202e_to_240e_19n_to_40n'   
            elif campaign == 'SOCRATES':
                E3SMdomain_range = '133e_to_164e_42s_to_63s'   
            else:
                raise ValueError('this aircraft campaign is not recognized: ' + campaign)
            
            #%% read in E3SM data
            variables_out = list()
            pblh_out = list()
            for varname in variable_names:
                variables_out.append([])
            
            date2 = date[0:4] + '-' + date[4:6] + '-' + date[6:8]
            filename_input = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.' + date2 + '-00000.nc'
        
            (timem, lonm, timeunitm, lonmunit, lonmname) = read_E3SM(filename_input, 'lon_' + E3SMdomain_range)
            (timem, latm, timeunitm, latmunit, latmname) = read_E3SM(filename_input, 'lat_' + E3SMdomain_range)
            (timem, z3, timeunitm, zunit, zname) = read_E3SM(filename_input, 'Z3_' + E3SMdomain_range)
            (timem, pblh, timeunitm, pblhunit, pblhname) = read_E3SM(filename_input, 'PBLH_' + E3SMdomain_range)
            # read in all variables
            (timem, variables, timeunitm, var_units, var_longnames) = \
                         read_E3SM(filename_input, [a + '_' + E3SMdomain_range for a in variable_names])
                         
            # cdaym = timeunit2cday(timeunitm, 'noleap')
            # yearm = timeunitm.split(' ')[2][0:4]
            timem = 86400* (timem.data - int(timem[0]))
            
            for tt in range(len(time)):
                t_idx = np.abs(timem-time[tt]).argmin()
                x_idx = find_nearest(lonm, latm, lon[tt], lat[tt])
                z_idx = np.abs(z3[t_idx, :, x_idx]-height[tt]).argmin()
                for vv in range(varlen):
                    variables_out[vv].append(variables[vv][t_idx, z_idx, x_idx])
                pblh_out.append(pblh[t_idx, x_idx])
                
             # %% output extacted file
            outputname = 'Aircraft_vars_' + campaign + '_' + model + '_' + date + '.nc'
            print('output to this file: ' + E3SM_aircraft_path + outputname)
            
            # define filename
            f = Dataset(E3SM_aircraft_path + outputname, 'w', format = 'NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time",))
            height_o = f.createVariable("height", 'f8', ("time",))
            pblh_o = f.createVariable('PBLH', 'f8', ("time",))
            var_o = list()
            for vv in range(varlen):
                var_o.append (f.createVariable(variable_names[vv], 'f8', ("time", )))
            
            # write data
            time_o[:] = time
            height_o[:] = height
            pblh_o[:] = np.array(pblh_out)
            for vv in range(varlen):
                var_o[vv][:] = np.array(variables_out[vv])
            
            # attributes
            time_o.units = "Seconds since " + date2 + ' 00:00 UTC'
            height_o.units = 'm MSL'
            pblh_o.units = pblhunit
            pblh_o.long_name = pblhname
            for vv in range(varlen):
                var_o[vv].units = var_units[vv]
                var_o[vv].long_name = var_longnames[vv]
            
            # global attributes
            import time as ttt
            f.description = model + " extact for aircraft track for " + campaign
            f.aircraftfile = filename.split('\\')[-1]
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
            
