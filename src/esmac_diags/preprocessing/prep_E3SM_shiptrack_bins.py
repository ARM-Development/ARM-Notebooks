"""
# prepare E3SM surface aerosol size distribution at ARM ship-based field campaigns
# input data is E3SM regional output
# output is surface variables at the nearest grid of the ship track
"""

import os
import glob
import numpy as np
from ..subroutines.time_format_change import timeunit2cday, yyyymmdd2cday, cday2mmdd
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.read_ship import read_marmet
from ..subroutines.read_ARMdata import read_met
from ..subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm
from netCDF4 import Dataset

def find_nearest(xall, yall, x, y):
    distance = np.square(xall-x) + np.square(yall-y)
    idx = distance.argmin()
    return(idx)

def run_prep(settings):
    #%% variables from settings
    campaign = settings['campaign']
    shipmetpath = settings['shipmetpath']
    Model_List = settings['Model_List']
    E3SM_hourly_path = settings['E3SM_hourly_path']
    E3SM_hourly_filehead = settings['E3SM_hourly_filehead']
    E3SM_ship_path = settings['E3SM_ship_path']

    #%% other settings
        
    if not os.path.exists(E3SM_ship_path):
        os.makedirs(E3SM_ship_path)
        
    
    #%% get all ship data
    if campaign == 'MAGIC':
        lst = glob.glob(shipmetpath + 'marmet*.txt')
        E3SMdomain_range = '202e_to_243e_20n_to_35n'    # domain range in E3SM regional output
    elif campaign == 'MARCUS':
        lst = [1, 2, 3, 4]     # there are 4 ship trips (legs) for MARCUS
        E3SMdomain_range = '60e_to_160e_42s_to_70s'   
    else:
        raise ValueError('data for this field campaign is not specified: ' + campaign)
    lst.sort()
    print('total number of ship leg files:' + str(len(lst)))
    
       
    for filename in lst:
        
        
        #%% read in ship data
        
        if campaign == 'MAGIC':
            # for each ship leg
            legnum = filename[-6:-4]
            (shipdata, shipvarlist) = read_marmet(filename)
            year = [a[1] for a in shipdata]
            month = [a[2] for a in shipdata]
            day = [a[3] for a in shipdata]
            hh = [int(a[4]) for a in shipdata]
            mm = [int(a[5]) for a in shipdata]
            ss = [int(a[6]) for a in shipdata]
            lat = np.array([float(a[7]) for a in shipdata])
            lon = np.array([float(a[8]) for a in shipdata])
            
            # ymd = [year[i] + '-' + month[i] + '-' + day[i] for i in range(len(year))]   # yyyy-mm-dd
            yyyymmdd = [year[i] + month[i] + day[i] for i in range(len(year))]   # yyyymmdd
            ymd = list(set(yyyymmdd))  # unique date
            ymd.sort()
            
            
            time = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400. 
            for i in range(len(time)):
                cday0 = yyyymmdd2cday(yyyymmdd[i], 'noleap') 
                if year[i] == year[0]:
                    time[i] = time[i] + cday0
                else:
                    time[i] = time[i] + cday0 + 365  # next year
    
        elif campaign == 'MARCUS':
            legnum = str(filename)
            if legnum == '1':
                startdate = '2017-10-30'
                enddate = '2017-12-02'
            elif legnum == '2':
                startdate = '2017-12-13'
                enddate = '2018-01-11'
            elif legnum == '3':
                startdate = '2018-01-16'
                enddate = '2018-03-04'
            elif legnum == '4':
                startdate = '2018-03-09'
                enddate = '2018-03-22'
                
            cday1 = yyyymmdd2cday(startdate, 'noleap')
            cday2 = yyyymmdd2cday(enddate, 'noleap')
            if startdate[0:4]!=enddate[0:4]:
                cday2 = cday2 + 365  # cover two years
    
            time = np.empty(0)
            lon = np.empty(0)
            lat = np.empty(0)
            ymd = []
            for cc in range(cday1, cday2 + 1):
                if cc <= 365:
                    yyyymmdd = startdate[0:4] + cday2mmdd(cc)
                else:
                    yyyymmdd = enddate[0:4] + cday2mmdd(cc-365)
                    
                lst0 = glob.glob(shipmetpath + 'maraadmetX1.b1.' + yyyymmdd + '*')
                (time0, lon0, timeunit, lonunit, lon_long_name) = read_met(lst0[0], 'lon')
                (time0, lat0, timeunit, lonunit, lon_long_name) = read_met(lst0[0], 'lat')
                ymd0 = timeunit.split()[2]
                ymd.append(ymd0[0:4] + ymd0[5:7] + ymd0[8:10])
                
                time = np.hstack((time,  time0/86400. + cc))
                lat = np.hstack((lat, lat0))
                lon = np.hstack((lon, lon0))
    
        print('date for shipleg ' + legnum + ': ' + ymd[0] + '-' + ymd[-1])
        
        #%% read in E3SM data
        for mm in range(len(Model_List)):
            model = Model_List[mm]
            
            varname = ['T', 'Pres', 'num_a1', 'dgnd_a01', 'num_a2', 'dgnd_a02', 'num_a3', 'dgnd_a03', 'num_a4', 'dgnd_a04']
            if model == 'Nuc':  # with nucleation mode
                varname = varname + ['num_a5', 'dgnd_a05']
            variables = list()
            varlen = len(varname)
            
            # read all days in the ship leg
            for dd in range(len(ymd)):
                ymd2 = ymd[dd][0:4] + '-' + ymd[dd][4:6] + '-' + ymd[dd][6:8]
                print('read this date: ' + ymd2)
                filename_input = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.' + ymd2 + '-00000.nc'
                
                (timem, lonm, timeunitm, lonmunit, lonmname) = read_E3SM(filename_input, 'lon_' + E3SMdomain_range)
                (timem, latm, timeunitm, latmunit, latmname) = read_E3SM(filename_input, 'lat_' + E3SMdomain_range)
                # (timem, psm, timeunitm, psmunit, psmname) = read_E3SM(filename_input, 'PS_' + E3SMdomain_range)
                
                cdaym = timeunit2cday(timeunitm, 'noleap')
                yearm = timeunitm.split(' ')[2][0:4]
                timem2 = timem.data - 365*(int(ymd[0][0:4]) - int(yearm)) + cdaym
                
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
                if model == 'Nuc':  # with nucleation mode
                    num_a5 = f.variables['num_a5_' + E3SMdomain_range][:]
                    dn5 = f.variables['dgnd_a05_' + E3SMdomain_range][:]
                f.close()
            
                Pres = np.nan*T
                zlen = T.shape[1]
                for kk in range(zlen):
                    Pres[:, kk, :] = hyam[kk]*P0 + hybm[kk]*PS
            
                varall = [T, Pres, num_a1, dn1, num_a2, dn2, num_a3, dn3, num_a4, dn4]
                if model == 'Nuc':  # with nucleation mode
                    varall.append(num_a5)
                    varall.append(dn5)
                    
                # ship measurement times during the model day
                timeo = time[np.logical_and(time >= timem2[0], time < timem2[0] + 1)]
                lono = lon[np.logical_and(time >= timem2[0], time < timem2[0] + 1)]
                lato = lat[np.logical_and(time >= timem2[0], time < timem2[0] + 1)]
              
                    
                # allocation variables and attributes
                if dd == 0:
                    for vv in range(varlen):
                        variables.append([])
                        
                # extract the data at the time and location of ship
                for tt in range(len(timeo)):
                    t_idx = np.abs(timem2-timeo[tt]).argmin()
                    if lono[tt]<-900. or lato[tt]<-900:
                        for vv in range(varlen):
                            variables[vv].append(np.nan)
                    else:
                        x_idx = find_nearest(lonm, latm, lono[tt], lato[tt])
                        for vv in range(varlen):
                            variables[vv].append(varall[vv][t_idx, -1, x_idx])  # choose the lowest level
    
            numall = [np.array(a) for a in variables[2::2]]
            dnall = [np.array(a) for a in variables[3::2]]
            
            NCNall = calc_CNsize_cutoff_0_3000nm(dnall, numall, np.array(variables[0]), np.array(variables[1]))
    
            # calculate total CN concentration for CPC (>10nm) and CPCU (>3nm)
            NUCN = np.nansum(NCNall[3:, :], 0)    # >3nm
            NCN = np.nansum(NCNall[10:, :], 0)    # >10nm
            
            
            # %% output extacted file
            outputname = 'Ship_CNsize_' + campaign + '_' + model + '_shipleg' + legnum + '.nc'
            print('output to this file: ' + E3SM_ship_path + outputname)
            
            # define filename
            f = Dataset(E3SM_ship_path + outputname, 'w', format='NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            size = f.createDimension('size', 3000)
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time", ))
            size_o = f.createVariable("size", 'i8', ("size"))
            lat_o = f.createVariable("lat", "f8", ("time", ))
            lon_o = f.createVariable("lon", "f8", ("time", ))
            
            data_o = f.createVariable('NCNall', 'f8', ("size", "time"))
            ncn_o = f.createVariable("NCN", "f8", ("time", ))
            nucn_o = f.createVariable("NUCN", "f8", ("time", ))
            
            # write data
            time_o[:] = time
            size_o[:] = np.arange(1, 3001)
            lat[lat<-900] = -9999.
            lon[lon<-900] = -9999.
            lat_o[:] = lat
            lon_o[:] = lon
            NCNall[np.isnan(NCNall)] = -9999.
            NCN[np.isnan(NCN)] = -9999.
            NUCN[np.isnan(NUCN)] = -9999.
            data_o[:, :] = NCNall
            ncn_o[:] = NCN
            nucn_o[:] = NUCN
            
            # attributes
            time_o.units = "days since " + str(int(ymd[0][0:4])-1) + "-12-31 00:00:00 UTC"
            lat_o.units = "degree north"
            lon_o.units = "degree east"
            time_o.long_name = "Calendar Day"
            lat_o.long_name = "latitude"
            lon_o.long_name = "longitude"
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
            f.description = model + " calculated aerosol size distribution along ship tracks for " + campaign
            f.shiptrackdata = filename
            f.modeldata = E3SM_hourly_path[mm] + E3SM_hourly_filehead[mm] + '.cam.h3.*.nc'
            f.datanotes = 'variables are set as missing if GPS location is missing'
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()