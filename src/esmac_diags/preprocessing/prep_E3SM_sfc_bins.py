"""
# prepare E3SM surface aerosol size distribution at ARM sites
# input data is E3SM regional output
# output is surface aerosol distribution at the nearest grid
"""

import os
import numpy as np
from ..subroutines.time_format_change import timeunit2cday, yyyymmdd2cday, cday2mmdd
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm
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
        
    #%% process data for each day
    for mm in range(len(Model_List)):
        model = Model_List[mm]
        
        for cday in range(cday1, cday2 + 1):
            mmdd = cday2mmdd(cday)
            date = year0 + '-' + mmdd[0:2] + '-' + mmdd[2:4]
            
            filename_input = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.' + date + '-00000.nc'
            
            (timem, lonm, timeunitm, lonmunit, lonmname) = read_E3SM(filename_input, 'lon_' + E3SMdomain_range)
            (timem, latm, timeunitm, latmunit, latmname) = read_E3SM(filename_input, 'lat_' + E3SMdomain_range)
            x_idx = find_nearest(lonm, latm, lon0, lat0)
            
            cdaym = timeunit2cday(timeunitm, 'noleap')
            yearm = timeunitm.split(' ')[2][0:4]
            time = timem.data - 365*(int(year0)-int(yearm)) + cdaym
                
            # do not use read_E3SM because hyam and hybm don't have units
            f = Dataset(filename_input, 'r')
            P0 = f.variables['P0'][:]
            hyam = f.variables['hyam'][:]
            hybm = f.variables['hybm'][:]
            T = f.variables['T_' + E3SMdomain_range][:]
            PS = f.variables['PS_' + E3SMdomain_range][:]
            num_a1 = f.variables['num_a1_' + E3SMdomain_range][:]
            num_a2 = f.variables['num_a2_' + E3SMdomain_range][:]
            num_a3 = f.variables['num_a3_' + E3SMdomain_range][:]
            num_a4 = f.variables['num_a4_' + E3SMdomain_range][:]
            dn1 = f.variables['dgnd_a01_' + E3SMdomain_range][:]
            dn2 = f.variables['dgnd_a02_' + E3SMdomain_range][:]
            dn3 = f.variables['dgnd_a03_' + E3SMdomain_range][:]
            dn4 = f.variables['dgnd_a04_' + E3SMdomain_range][:]
            if model[0:3] == 'Nuc':  # with nucleation mode
                num_a5 = f.variables['num_a5_' + E3SMdomain_range][:]
                dn5 = f.variables['dgnd_a05_' + E3SMdomain_range][:]
            f.close()
        
            Pres = np.nan*T
            zlen = T.shape[1]
            for kk in range(zlen):
                Pres[:, kk, :] = hyam[kk]*P0  +  hybm[kk]*PS
        
            numall = [num_a1[:, -1, x_idx], num_a2[:, -1, x_idx], num_a3[:, -1, x_idx], num_a4[:, -1, x_idx]]
            dnall = [dn1[:, -1, x_idx], dn2[:, -1, x_idx], dn3[:, -1, x_idx], dn4[:, -1, x_idx]]
            if model[0:3] == 'Nuc':  # with nucleation mode
                numall.append(num_a5[:, -1, x_idx])
                dnall.append(dn5[:, -1, x_idx])
               
            
            NCNall = calc_CNsize_cutoff_0_3000nm(dnall, numall, T[:, -1, x_idx], Pres[:, -1, x_idx])
    
            # calculate total CN concentration for CPC (>10nm) and CPCU (>3nm)
            NUCN = np.nansum(NCNall[3:, :], 0)    # >3nm
            NCN = np.nansum(NCNall[10:, :], 0)    # >10nm
            
            
        
            #%% output extacted file
            outputname = 'SFC_CNsize_' + campaign + '_' + model + '_' + date + '.nc'
            print('output to this file: ' + E3SM_sfc_path + outputname)
            
            # define filename
            f = Dataset(E3SM_sfc_path + outputname, 'w', format='NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            s = f.createDimension('size', 3000)  # unlimited
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time", ))
            size_o = f.createVariable("size", "f8", ("size", ))
            lat_o = f.createVariable("lat", "f8", ())
            lon_o = f.createVariable("lon", "f8", ())
            
            data_o = f.createVariable('NCNall', 'f8', ("size", "time"))
            ncn_o = f.createVariable("NCN", "f8", ("time", ))
            nucn_o = f.createVariable("NUCN", "f8", ("time", ))
            
            # write data
            time_o[:] = time
            lat_o[:] = latm[x_idx]
            lon_o[:] = lonm[x_idx]
            size_o[:] = np.arange(1, 3001)
            data_o[:, :] = NCNall
            ncn_o[:] = NCN
            nucn_o[:] = NUCN
            
            # attributes
            time_o.units = "days since " + str(int(year0)-1) + "-12-31 00:00:00 UTC"
            lat_o.units = "latitude"
            lon_o.units = "longitude"
            size_o.units = 'nm'
            size_o.long_name = "0 to 3000nm with 1nm increment"
            data_o.units = '#/m3'
            data_o.long_name = 'aerosol size distribution'
            ncn_o.units = '#/m3'
            ncn_o.long_name = 'aerosol number concentration for size >10nm'
            nucn_o.units = '#/m3'
            nucn_o.long_name = 'aerosol number concentration for size >3nm'
            
            # global attributes
            import time as ttt
            f.description = model + " extact surface aerosol size distribution for " + campaign
            f.modeldata = filename_input
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
