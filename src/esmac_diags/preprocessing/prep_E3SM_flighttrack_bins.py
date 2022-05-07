"""
# prepare E3SM aerosol size distribution for flight tracks
# input data is IWG measurements from aircraft and E3SM regional output
# output is  aerosol size distribution for each flight
"""

import glob
import os
import numpy as np
from ..subroutines.time_format_change import hhmmss2sec, timeunit2cday
from ..subroutines.read_aircraft import read_iwg1, read_RF_NCAR
from ..subroutines.read_netcdf import read_E3SM
from ..subroutines.CN_mode_to_size import calc_CNsize_cutoff_0_3000nm
from netCDF4 import Dataset


def find_nearest(xall,yall,x,y):
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
    
    
    if campaign=='HISCALE':
        E3SMdomain_range='260e_to_265e_34n_to_39n'    # domain range in E3SM regional output
    elif campaign=='ACEENA':
        E3SMdomain_range='330e_to_335e_37n_to_42n'   
    elif campaign=='CSET':
        E3SMdomain_range='202e_to_240e_19n_to_40n'   
    elif campaign=='SOCRATES':
        E3SMdomain_range='133e_to_164e_42s_to_63s'  
    else:
        raise ValueError('this aircraft campaign is not recognized: '+campaign)
    
    #%% find all flight data
    if campaign=='HISCALE':
        lst = glob.glob(iwgpath+'*a2.txt')
        lst.sort()
        if IOP=='IOP1':
            lst=lst[0:17]
        elif IOP=='IOP2':
            lst=lst[17:]
        elif IOP[0:4]=='2016':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    elif campaign=='ACEENA':
        lst = glob.glob(iwgpath+'*a2.txt')
        lst.sort()
        if IOP=='IOP1':
            lst=lst[0:20]
        elif IOP=='IOP2':
            lst=lst[20:]
        elif IOP[0:4]=='2017' or IOP[0:4]=='2018':
            a=lst[0].split('_'+campaign+'_')
            lst = glob.glob(a[0]+'*'+IOP+'*')
            lst.sort()
    elif campaign in ['CSET', 'SOCRATES']:
        lst = glob.glob(RFpath+'RF*.PNI.nc')
        lst.sort()
    else:
        raise ValueError('this aircraft campaign is not recognized: '+campaign)
        
    print('total number of files:'+str(len(lst)))
    
    for filename in lst:
        
        fname=filename.split('.')
        #%% read in flight data
        if campaign in ['HISCALE', 'ACEENA']:
            date=fname[-3]
            print('input data for '+date)
            # year=date[0:4]
            # month=date[4:6]
            
            (flight,flightvars)=read_iwg1(filename)
            timelen = len(flight)
            # get lat, lon, height, time
            lon=np.empty(timelen)
            lat=np.empty(timelen)
            height=np.empty(timelen)
            time=np.empty(timelen)
            if np.logical_and(campaign=='ACEENA', date=='20180216a'):
                flight.insert(1403,list(flight[1403]))
                tstr=flight[1403][1]
                tstr=tstr[0:-1]+str(int(tstr[-1])-1)
                flight[1403][1]=tstr
                del flight[-1]
            for t in range(timelen):
                lat[t]=float(flight[t][2])
                lon[t]=float(flight[t][3])+360
                height[t]=float(flight[t][4])
                timestr=flight[t][1].split(' ')
                time[t]=hhmmss2sec(timestr[1])
        
        elif campaign in ['CSET', 'SOCRATES']:
            date=fname[-4]
            print('input data for '+date)
            # year=date[0:4]
            # month=date[4:6]
            (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'ALT')
            (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LAT')
            (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'LON')
            lon[lon<0]=lon[lon<0]+360
            timelen = len(time)
        
            #%% read in E3SM data
        for mm in range(len(Model_List)):
            model=Model_List[mm]  
            date2 = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
            filename_input = E3SM_hourly_path[mm]+E3SM_hourly_filehead[mm]+'.cam.h3.'+date2+'-00000.nc'
        
            (timem,lonm,timeunitm,lonmunit,lonmname)=read_E3SM(filename_input,'lon_'+E3SMdomain_range)
            (timem,latm,timeunitm,latmunit,latmname)=read_E3SM(filename_input,'lat_'+E3SMdomain_range)
            (timem,z3,timeunitm,latmunit,latmname)=read_E3SM(filename_input,'Z3_'+E3SMdomain_range)
            # do not use read_E3SM because hyam and hybm don't have units
            f = Dataset(filename_input,'r')
            P0 = f.variables['P0'][:]
            hyam = f.variables['hyam'][:]
            hybm = f.variables['hybm'][:]
            T = f.variables['T_'+E3SMdomain_range][:]
            PS = f.variables['PS_'+E3SMdomain_range][:]
            num_a1 = f.variables['num_a1_'+E3SMdomain_range][:]
            num_a2 = f.variables['num_a2_'+E3SMdomain_range][:]
            num_a3 = f.variables['num_a3_'+E3SMdomain_range][:]
            num_a4 = f.variables['num_a4_'+E3SMdomain_range][:]
            dn1 = f.variables['dgnd_a01_'+E3SMdomain_range][:]
            dn2 = f.variables['dgnd_a02_'+E3SMdomain_range][:]
            dn3 = f.variables['dgnd_a03_'+E3SMdomain_range][:]
            dn4 = f.variables['dgnd_a04_'+E3SMdomain_range][:]
            if model[0:3]=='Nuc':  # with nucleation mode
                num_a5 = f.variables['num_a5_'+E3SMdomain_range][:]
                dn5 = f.variables['dgnd_a05_'+E3SMdomain_range][:]
            f.close()
            
            Pres = np.nan*T
            zlen=T.shape[1]
            for kk in range(zlen):
                Pres[:,kk,:] = hyam[kk]*P0 + hybm[kk]*PS
        
            #% find the nearest time and height of the aircraft measurements
            cdaym = timeunit2cday(timeunitm,'noleap')
            timem = 86400* (timem.data - int(timem[0]))
            NCNall=np.full((3000,timelen),np.nan)
            tzx0 = [0,0,0]
            t0 = 0
            for tt in range(timelen):
                t_idx = np.abs(timem-time[tt]).argmin()
                x_idx = find_nearest(lonm,latm,lon[tt],lat[tt])
                z_idx = np.abs(z3[t_idx,:,x_idx]-height[tt]).argmin()
                
                # copy the same grid to save time
                if [t_idx,x_idx,z_idx]==tzx0:
                    NCNall[:,tt] = NCNall[:,t0]
                else:
                    numall = [num_a1[t_idx,z_idx,x_idx],num_a2[t_idx,z_idx,x_idx],num_a3[t_idx,z_idx,x_idx],num_a4[t_idx,z_idx,x_idx]]
                    dnall = [dn1[t_idx,z_idx,x_idx],dn2[t_idx,z_idx,x_idx],dn3[t_idx,z_idx,x_idx],dn4[t_idx,z_idx,x_idx]]
                    if model[0:3]=='Nuc':  # with nucleation mode
                        numall.append(num_a5[t_idx,z_idx,x_idx])
                        dnall.append(dn5[t_idx,z_idx,x_idx])
                    NCNall[:,tt] = calc_CNsize_cutoff_0_3000nm(dnall,numall,T[t_idx,z_idx,x_idx],Pres[t_idx,z_idx,x_idx])
                    # update the time of this unique grid
                    tzx0=[t_idx,x_idx,z_idx]
                    t0=tt
    
            # calculate total CN concentration for CPC (>10nm) and CPCU (>3nm)
            NUCN = np.nansum(NCNall[3:,:],0)    # >3nm
            NCN = np.nansum(NCNall[10:,:],0)    # >10nm
            
                
            #%% output extacted file
            outputname = 'Aircraft_CNsize_'+campaign+'_'+model+'_'+date+'.nc'
            print('output to this file: '+E3SM_aircraft_path+outputname)
            
            # define filename
            f = Dataset(E3SM_aircraft_path+outputname, 'w', format='NETCDF4')
            
            # define dimensions
            t = f.createDimension('time', None)  # unlimited
            size=f.createDimension('size',3000)
            
            # create variable list
            time_o = f.createVariable("time","f8",("time",))
            height_o = f.createVariable("height",'f8',("time",))
            size_o = f.createVariable("size",'i8',("size"))
            
            data_o = f.createVariable('NCNall','f8',("size","time"))
            ncn_o = f.createVariable("NCN","f8",("time",))
            nucn_o = f.createVariable("NUCN","f8",("time",))
            
            # write data
            time_o[:] = time
            height_o[:] = height
            size_o[:] = np.arange(1,3001)
            data_o[:,:]=NCNall
            ncn_o[:]=NCN
            nucn_o[:]=NUCN
            
            # attributes
            time_o.units = "Seconds since "+date2+' 00:00 UTC'
            height_o.units = 'm MSL'
            size_o.units = 'nm'
            size_o.long_name="0 to 3000nm with 1nm increment"
            data_o.units = '#/m3'
            data_o.long_name = 'aerosol size distribution'
            ncn_o.units = '#/m3'
            ncn_o.long_name = 'aerosol number concentration for size >10nm'
            nucn_o.units = '#/m3'
            nucn_o.long_name = 'aerosol number concentration for size >3nm'
            
            # global attributes
            import time as ttt
            f.description = model+" extact for aircraft track for "+campaign
            f.aircraftfile = filename.split('\\')[-1]
            f.create_time = ttt.ctime(ttt.time())
            
            f.close()
    
