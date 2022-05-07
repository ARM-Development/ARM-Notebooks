"""
# prepare E3SM surface variables at ARM ship-based field campaigns
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
    if campaign=='MAGIC':
        lst = glob.glob(shipmetpath+'marmet*.txt')
        E3SMdomain_range='202e_to_243e_20n_to_35n'    # domain range in E3SM regional output
    elif campaign=='MARCUS':
        lst = [1, 2, 3, 4]     # there are 4 ship trips (legs) for MARCUS
        E3SMdomain_range='60e_to_160e_42s_to_70s'   
    else:
        raise ValueError('data for this field campaign is not specified: ' + campaign)
    
    lst.sort()
    print('total number of ship leg files:'+str(len(lst)))
    
    
    for filename in lst:
    
        #%% read in ship data
        
        if campaign=='MAGIC':
            # for each ship leg
            legnum=filename[-6:-4]
            
            (shipdata, shipvarlist) = read_marmet(filename)
            year=[a[1] for a in shipdata]
            month=[a[2] for a in shipdata]
            day=[a[3] for a in shipdata]
            hh=[int(a[4]) for a in shipdata]
            mm=[int(a[5]) for a in shipdata]
            ss=[int(a[6]) for a in shipdata]
            lat=np.array([float(a[7]) for a in shipdata])
            lon=np.array([float(a[8]) for a in shipdata])
            
            # ymd = [year[i]+'-'+month[i]+'-'+day[i] for i in range(len(year))]   # yyyy-mm-dd
            yyyymmdd = [year[i]+month[i]+day[i] for i in range(len(year))]   # yyyymmdd
            ymd=list(set(yyyymmdd))  # unique date
            ymd.sort()
            
            
            time = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400. 
            for i in range(len(time)):
                cday0 = yyyymmdd2cday(yyyymmdd[i], 'noleap') 
                if year[i]==year[0]:
                    time[i]=time[i]+cday0
                else:
                    time[i]=time[i]+cday0+365  # next year
            
        elif campaign=='MARCUS':
            legnum=str(filename)
            if legnum=='1':
                startdate='2017-10-30'
                enddate='2017-12-02'
            elif legnum=='2':
                startdate='2017-12-13'
                enddate='2018-01-11'
            elif legnum=='3':
                startdate='2018-01-16'
                enddate='2018-03-04'
            elif legnum=='4':
                startdate='2018-03-09'
                enddate='2018-03-22'
                
            cday1=yyyymmdd2cday(startdate, 'noleap')
            cday2=yyyymmdd2cday(enddate, 'noleap')
            if startdate[0:4]!=enddate[0:4]:
                cday2=cday2+365  # cover two years
    
            time=np.empty(0)
            lon=np.empty(0)
            lat=np.empty(0)
            ymd=[]
            for cc in range(cday1, cday2+1):
                if cc<=365:
                    yyyymmdd=startdate[0:4]+cday2mmdd(cc)
                else:
                    yyyymmdd=enddate[0:4]+cday2mmdd(cc-365)
                    
                lst0 = glob.glob(shipmetpath+'maraadmetX1.b1.'+yyyymmdd+'*')
                (time0, lon0, timeunit, lonunit, lon_long_name)=read_met(lst0[0], 'lon')
                (time0, lat0, timeunit, lonunit, lon_long_name)=read_met(lst0[0], 'lat')
                ymd0 = timeunit.split()[2]
                ymd.append(ymd0[0:4]+ymd0[5:7]+ymd0[8:10])
                
                time = np.hstack((time, time0/86400. + cc))
                lat = np.hstack((lat, lat0))
                lon = np.hstack((lon, lon0))
        
        print('date for shipleg '+legnum+': '+ymd[0]+'-'+ymd[-1])
        
        #%% set variables to be read
        for mm in range(len(Model_List)):
            model=Model_List[mm]
            variable1d_names = ['PS', 'PBLH', 'FLNT', 'FSNT', 'FLNS', 'FSNS', "LHFLX", "SHFLX",
                                 'TREFHT', 'PRECT', 'PRECL', "TGCLDLWP", "TGCLDIWP"]
            variable2d_names = ['T', 'U', 'V', 'Q', 'RELHUM', 'RHW', 'RHI', 'CLOUD',
                             'CLDLIQ', 'CLDICE', 'NUMLIQ', 'NUMICE', 'CCN1', 'CCN3', 'CCN5', "AREI", "AREL",
                             'bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', 'mom_a2', 'mom_a3', 'mom_a4',
                             'ncl_a1', 'ncl_a2', 'ncl_a3', 'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
                             'soa_a1', 'soa_a2', 'soa_a3', 'num_a1', 'num_a2', 'num_a3', 'num_a4',
                             'num_c1', 'num_c2', 'num_c3', 'num_c4', "dgnd_a01", "dgnd_a02", "dgnd_a03", "dgnd_a04", 
                             "dgnw_a01", "dgnw_a02", "dgnw_a03", "dgnw_a04", 'EXTINCT', 'ABSORB']
            if model=='NucSoaCond': # with so4 and soa in nucleation mode
                variable2d_names=variable2d_names+['so4_a5', 'soa_a5', 'num_a5', 'num_c5', "dgnd_a05", "dgnw_a05"]
            elif model=='Nuc':      # only with so4 in nucleation mode
                variable2d_names=variable2d_names+['so4_a5', 'num_a5', 'num_c5', "dgnd_a05", "dgnw_a05"]
            var1dlen = len(variable1d_names)
            var2dlen = len(variable2d_names)
            variable_names = variable1d_names+variable2d_names
            varlen = var1dlen+var2dlen     
        
            #%% read in E3SM data
            variables = list()
            var_units = list()
            var_longnames = list()
            
            # read all days in the ship leg
            for dd in range(len(ymd)):
                ymd2 = ymd[dd][0:4]+'-'+ymd[dd][4:6]+'-'+ymd[dd][6:8]
                print('read this date: '+ymd2)
                filename_input = E3SM_hourly_path[mm]+E3SM_hourly_filehead[mm]+'.cam.h3.'+ymd2+'-00000.nc'
                
                (timem, lonm, timeunitm, lonmunit, lonmname)=read_E3SM(filename_input, 'lon_'+E3SMdomain_range)
                (timem, latm, timeunitm, latmunit, latmname)=read_E3SM(filename_input, 'lat_'+E3SMdomain_range)
                # (timem, psm, timeunitm, psmunit, psmname)=read_E3SM(filename_input, 'PS_'+E3SMdomain_range)
                
                cdaym = timeunit2cday(timeunitm, 'noleap')
                yearm = timeunitm.split(' ')[2][0:4]
                timem2 = timem.data-365*(int(ymd[0][0:4])-int(yearm)) + cdaym
                
                # ship measurement times during the model day
                timeo = time[np.logical_and(time>=timem2[0], time<timem2[0]+1)]
                lono = lon[np.logical_and(time>=timem2[0], time<timem2[0]+1)]
                lato = lat[np.logical_and(time>=timem2[0], time<timem2[0]+1)]
                
                (timem, var1d, timeunitm, var1dunit, var1dlongname) = \
                    read_E3SM(filename_input, [a+'_'+E3SMdomain_range for a in variable1d_names])
                (timem, var2d, timeunitm, var2dunit, var2dlongname) = \
                    read_E3SM(filename_input, [a+'_'+E3SMdomain_range for a in variable2d_names])
                
                # allocation variables and attributes
                if dd==0:
                    for vv in range(varlen):
                        variables.append([])
                    var_units=var1dunit+var2dunit
                    var_longnames=var1dlongname+var2dlongname
                    
                for tt in range(len(timeo)):
                    t_idx = np.abs(timem2-timeo[tt]).argmin()
                    if lono[tt]<-900. or lato[tt]<-900:
                        for vv in range(varlen):
                            variables[vv].append(-9999.)
                    else:
                        x_idx = find_nearest(lonm, latm, lono[tt], lato[tt])
                        for vv in range(var1dlen):
                            variables[vv].append(var1d[vv][t_idx, x_idx])
                        for vv in range(var2dlen):
                            variables[var1dlen+vv].append(var2d[vv][t_idx, -1, x_idx])  # choose the lowest level
    
            
            # %% output extacted file
            outputname = 'Ship_vars_'+campaign+'_'+model+'_shipleg'+legnum+'.nc'
            print('output to this file: '+E3SM_ship_path+outputname)
            
            # define filename
            f = Dataset(E3SM_ship_path+outputname, 'w', format='NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            
            # create variable list
            time_o = f.createVariable("time", "f8", ("time", ))
            lat_o = f.createVariable("lat", "f8", ("time", ))
            lon_o = f.createVariable("lon", "f8", ("time", ))
            var_o=list()
            for vv in range(varlen):
                var_o.append (f.createVariable(variable_names[vv], 'f8', ("time", )))
            
            # write data
            time_o[:] = time
            lat[lat<-900]=-9999.
            lon[lon<-900]=-9999.
            lat_o[:] = lat
            lon_o[:] = lon
            for vv in range(varlen):
                var_o[vv][:] = np.array(variables[vv])
            
            # attributes
            time_o.units = "days since "+str(int(ymd[0][0:4])-1)+"-12-31 00:00:00 UTC"
            lat_o.units = "degree north"
            lon_o.units = "degree east"
            time_o.long_name = "Calendar Day"
            lat_o.long_name = "latitude"
            lon_o.long_name = "longitude"
            for vv in range(varlen):
                var_o[vv].units = var_units[vv]
                var_o[vv].long_name = var_longnames[vv]
                var_o[vv].missing_value = -9999.
            
            # global attributes
            import time as ttt
            f.description = model+" extact variables along ship tracks for "+campaign
            f.shiptrackdata = filename
            f.modeldata = E3SM_hourly_path[mm]+E3SM_hourly_filehead[mm]+'.cam.h3.*.nc'
            f.datanotes = 'variables are set as missing if GPS location is missing'
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
            
