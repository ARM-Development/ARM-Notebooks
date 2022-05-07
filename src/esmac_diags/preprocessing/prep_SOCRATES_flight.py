"""
prepare aircraft data from ACEENA
options of average data into coarser resolution
"""

import glob
import os
import numpy as np
import time as ttt
from netCDF4 import Dataset
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, avg_time_2d, \
                    median_time_1d, median_time_2d
from esmac_diags.subroutines.read_aircraft import read_RF_NCAR, read_ccn_socrates
from esmac_diags.subroutines.quality_control import qc_cn_max, qc_remove_neg, \
                     qc_mask_cloudflag, qc_uhsas_RF_NCAR
from esmac_diags.subroutines.specific_data_treatment import lwc2cflag


# RFpath = '../../../data/SOCRATES/obs/aircraft/aircraft_lowrate/'
# ccnpath = '../../../data/SOCRATES/obs/aircraft/CCN/'
# prep_data_path = 'C:/Users/tang357/Downloads/SOCRATES/'

# dt=60

            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN(ccnpath, RFpath, prep_data_path, dt=60):
    """
    prepare CCN number concentration from NCAR research flight (RF) low-rate data
    
    Parameters
    ----------
    CCNpath : str
        input path for CCN dta from NCAR RF measurements
    RFpath : str
        input path for NCAR RF data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    #%% find all data
    lst = glob.glob(RFpath+'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLON')
        (time,p,timeunit,punit,plongname,cellsize,cellunit)=read_RF_NCAR(filename,'PSXC')
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        p = p.filled(np.nan)
    
        # read ccn data
        filename_ccn1=glob.glob(ccnpath+'CCNscanning_SOCRATES_GV_RF*'+date[0:8]+'_R0.ict')
        filename_ccn2=glob.glob(ccnpath+'CCNspectra_SOCRATES_GV_RF*'+date[0:8]+'_R0.ict')
        if len(filename_ccn2)==0:
            print('does not find any CCN data, skip this date...')
            continue
            
        (data1,ccn1list)=read_ccn_socrates(filename_ccn1[0])
        time_scan = data1[0,:]
        ccn_scan = data1[1,:]
        SS_scan = data1[3,:]
        (data2,ccn2list)=read_ccn_socrates(filename_ccn2[0])
        time1 = data2[0,:]
        time2 = data2[1,:]
        ccn_1 = data2[2,:]
        ccn_2 = data2[6,:]
        ccn_3 = data2[10,:]
        ccn_5 = data2[18,:]
        time_spec = 0.5*(time1+time2)
            
        # quality controls
        p = qc_remove_neg(p)
        ccn_scan = qc_remove_neg(ccn_scan)
        SS_scan = qc_remove_neg(SS_scan)
        ccn_1 = qc_remove_neg(ccn_1)
        ccn_2 = qc_remove_neg(ccn_2)
        ccn_3 = qc_remove_neg(ccn_3)
        ccn_5 = qc_remove_neg(ccn_5)
        
        ccn_1s = np.array(ccn_scan)
        ccn_2s = np.array(ccn_scan)
        ccn_3s = np.array(ccn_scan)
        ccn_5s = np.array(ccn_scan)
        ccn_1s[np.logical_or(SS_scan<0.05, SS_scan>0.15)] = np.nan
        ccn_2s[np.logical_or(SS_scan<0.15, SS_scan>0.25)] = np.nan
        ccn_3s[np.logical_or(SS_scan<0.25, SS_scan>0.35)] = np.nan
        ccn_5s[np.logical_or(SS_scan<0.45, SS_scan>0.55)] = np.nan
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        ccn1s_1 = median_time_1d(time_scan, ccn_1s, time_new)
        ccn2s_1 = median_time_1d(time_scan, ccn_2s, time_new)
        ccn3s_1 = median_time_1d(time_scan, ccn_3s, time_new)
        ccn5s_1 = median_time_1d(time_scan, ccn_5s, time_new)
        # ccn1_1 = np.interp(time_new, time_spec, ccn_1)
        # ccn2_1 = np.interp(time_new, time_spec, ccn_2)
        # ccn3_1 = np.interp(time_new, time_spec, ccn_3)
        # ccn5_1 = np.interp(time_new, time_spec, ccn_5)
        ccn1_1 = median_time_1d(time_spec, ccn_1, time_new)
        ccn2_1 = median_time_1d(time_spec, ccn_2, time_new)
        ccn3_1 = median_time_1d(time_spec, ccn_3, time_new)
        ccn5_1 = median_time_1d(time_spec, ccn_5, time_new)
     
        # for SORCRATES, CCN data is measured in CCN chamber where pressure is set as 400hPa
        # convert to ambient pressure. 
        # no information about temperature so assume it is ambient temperature
        p_new = median_time_1d(time, p, time_new)
        ccn1s_p = ccn1s_1 * p_new/400
        ccn2s_p = ccn2s_1 * p_new/400
        ccn3s_p = ccn3s_1 * p_new/400
        ccn5s_p = ccn5s_1 * p_new/400
        ccn1_p = ccn1_1 * p_new/400
        ccn2_p = ccn2_1 * p_new/400
        ccn3_p = ccn3_1 * p_new/400
        ccn5_p = ccn5_1 * p_new/400
        
        #%% output data
        outfile = prep_data_path + 'CCN_SOCRATES_' + date + '.nc'
                
        if not os.path.exists(prep_data_path):
            os.makedirs(prep_data_path)
            
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        c1_o = f.createVariable('CCN1_spec', 'f8', ("time",))
        c2_o = f.createVariable('CCN2_spec', 'f8', ("time",))
        c3_o = f.createVariable('CCN3_spec', 'f8', ("time",))
        c5_o = f.createVariable('CCN5_spec', 'f8', ("time",))
        c1s_o = f.createVariable('CCN1_scan', 'f8', ("time",))
        c2s_o = f.createVariable('CCN2_scan', 'f8', ("time",))
        c3s_o = f.createVariable('CCN3_scan', 'f8', ("time",))
        c5s_o = f.createVariable('CCN5_scan', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        c1_o[:] = ccn1_p
        c2_o[:] = ccn2_p
        c3_o[:] = ccn3_p
        c5_o[:] = ccn5_p
        c1s_o[:] = ccn1s_p
        c2s_o[:] = ccn2s_p
        c3s_o[:] = ccn3s_p
        c5s_o[:] = ccn5s_p
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        c1_o.units = '#/cm3'
        c1_o.long_name = 'CCN number concentration'
        c1_o.description = 'extracted from CCN spectra data for SS=0.1%'
        c2_o.units = '#/cm3'
        c2_o.long_name = 'CCN number concentration'
        c2_o.description = 'extracted from CCN spectra data for SS=0.2%'
        c3_o.units = '#/cm3'
        c3_o.long_name = 'CCN number concentration'
        c3_o.description = 'extracted from CCN spectra data for SS=0.3%'
        c5_o.units = '#/cm3'
        c5_o.long_name = 'CCN number concentration'
        c5_o.description = 'extracted from CCN spectra data for SS=0.5%'
        c1s_o.units = '#/cm3'
        c1s_o.long_name = 'CCN number concentration'
        c1s_o.description = 'extracted from scanning CCN data for SS between 0.05% and 0.15%'
        c2s_o.units = '#/cm3'
        c2s_o.long_name = 'CCN number concentration'
        c2s_o.description = 'extracted from scanning CCN data for SS between 0.15% and 0.25%'
        c3s_o.units = '#/cm3'
        c3s_o.long_name = 'CCN number concentration'
        c3s_o.description = 'extracted from scanning CCN data for SS between 0.25% and 0.35%'
        c5s_o.units = '#/cm3'
        c5s_o.long_name = 'CCN number concentration'
        c5s_o.description = 'extracted from scanning CCN data for SS between 0.45% and 0.55%'
    
        # global attributes
        f.title = "CCN number concentration"
        f.description = 'CCN is measured in CCN chamber where pressure is set as 400hPa. '+\
            'The measurements are converted into ambient pressure. '+\
            'Temperature is assumed in ambient condition.'
        f.description_scan = 'CCN*_scan is median value of each time window as the raw data has a resolution of 10s. '+\
            'only supersaturation falls within +/- 0.05% of the targetted SS is used'
        f.description_spectra = 'CCN*_spec save data in the nearest time point as the raw data has a resolution of 5min'
        f.input_file = [filename_ccn1[0].split('/')[-1], filename_ccn2[0].split('/')[-1]]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize(RFpath, prep_data_path, dt=60):
    """
    prepare aerosol size distribution from NCAR research flight (RF) low-rate data
    
    Parameters
    ----------
    RFpath : str
        input path for NCAR RF data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    #%% find all data
    lst = glob.glob(RFpath+'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLON')
        (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
        (time,uhsas,timeunit,dataunit,long_name,sizeh,cellunit)=read_RF_NCAR(filename,'CUHSAS_LWII')
        # calculate cloud flag based on LWC
        cflag=lwc2cflag(lwc,lwcunit)
        sizeh=sizeh*1000.
        sizel = np.hstack((2*sizeh[0]-sizeh[1],  sizeh[0:-1]))
        size = (sizeh+sizel)/2.
    
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        uhsas = qc_remove_neg(uhsas[:,0,:].data)
        uhsas = qc_mask_cloudflag(uhsas, cflag)
        uhsas = qc_uhsas_RF_NCAR(uhsas)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        # uhsas1 = avg_time_2d(time, uhsas, time_new)
        uhsas1 = median_time_2d(time, uhsas, time_new)
        # uhsas1 = median_time_forflight_2d(time, uhsas, time_new, height)
        
#     #%% 
#     import matplotlib.pyplot as plt
#     # fig = plt.figure(figsize=(8,2))
#     # ax1 = fig.add_subplot(1, 1, 1)
#     # ax1.plot(time/3600, lwcobs)
#     # ax1.plot(time_new/3600, lwc, color='r', marker='.',linewidth=2)
#     # ax1.set_ylim(-0.01, 0.8)
#     # ax2 = fig.add_subplot(2, 1, 2)
#     # ax2.plot(time/3600, uhsas100)
#     # ax2.plot(time_new/3600, uhsas, color='r', marker='.',linewidth=2)
#     # ax2.set_ylim(-50, 5000)
    
#     fig = plt.figure(figsize=(8,5))
#     ax1 = fig.add_subplot(2, 1, 1)
#     h1=ax1.contourf(time/3600, size, uhsas.T, np.arange(0,160, 10))
#     fig.colorbar(h1)
#     ax2 = fig.add_subplot(2, 1, 2)
#     h2=ax2.contourf(time_new/3600, size, uhsas1.T, np.arange(0,160, 10))
#     ax1.set_yscale('log')
#     ax2.set_yscale('log')
#     fig.colorbar(h2)
    
        #%% output data
        outfile = prep_data_path + 'UHSASsize_SOCRATES_' + date + '.nc'
        
        if not os.path.exists(prep_data_path):
            os.makedirs(prep_data_path)
            
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new)) 
        f.createDimension('size', len(size))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        merge_o = f.createVariable('size_distribution_uhsas', 'f8', ("time", "size"))
        
        # write data
        time_o[:] = time_new
        size_o[:] = size
        sizeh_o[:] = sizeh
        sizel_o[:] = sizel
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        merge_o[:] = uhsas1
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        size_o.units = 'nm'
        size_o.long_name = 'center of size bin'
        sizeh_o.units = 'nm'
        sizeh_o.long_name = 'upper bound of size bin'
        sizel_o.units = 'nm'
        sizel_o.long_name = 'lower bound of size bin'
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        merge_o.units = '#/cm3'
        merge_o.long_name = 'aerosol size distribution'
    
        # global attributes
        f.title = "Aerosol size distribution from UHSAS"
        f.description = 'median value of each time window'
        f.input_file = filename.split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CN(RFpath, prep_data_path, dt=60):
    """
    prepare aerosol number concentration from NCAR research flight (RF) low-rate data
    
    Parameters
    ----------
    RFpath : str
        input path for NCAR RF data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    #%% find all data
    lst = glob.glob(RFpath+'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLON')
        (time,cpc10,timeunit,cpc10unit,cpc10longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCN')
        (time,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_LWII')
        (time,uhsas100cvi,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_CVIU')
        (time,lwc,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
        # calculate cloud flag based on LWC
        cflag=lwc2cflag(lwc,lwcunit)
        
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        cpc10 = qc_remove_neg(cpc10.data)
        cpc10 = qc_cn_max(cpc10, 10)
        cpc10 = qc_mask_cloudflag(cpc10, cflag)
        uhsas100 = qc_remove_neg(uhsas100.data)
        uhsas100 = qc_cn_max(uhsas100, 100)
        uhsas100 = qc_mask_cloudflag(uhsas100, cflag)
        
        uhsas100cvi = qc_remove_neg(uhsas100cvi.data)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        cpc = median_time_1d(time, cpc10, time_new)
        uhsas = median_time_1d(time, uhsas100, time_new)
        # cpc = median_time_forflight_1d(time, cpc10, time_new, height)
        # uhsas = median_time_forflight_1d(time, uhsas100, time_new, height)
    
        #%% output data
        outfile = prep_data_path + 'CN_SOCRATES_' + date + '.nc'
        
        if not os.path.exists(prep_data_path):
            os.makedirs(prep_data_path)
            
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        c_o = f.createVariable('CPC10', 'f8', ("time",))
        u_o = f.createVariable('UHSAS100', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        c_o[:] = cpc
        u_o[:] = uhsas
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        c_o.units = cpc10unit
        c_o.long_name = 'CPC measured aerosol number (size>10nm)'
        u_o.units = uhsas100unit
        u_o.long_name = 'UHSAS measured aerosol number (size>100nm)'
    
        # global attributes
        f.title = "Aerosol number concentration"
        f.description = 'median value of each time window, in ambient condition'
        f.input_file = filename.split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LWC(RFpath, prep_data_path, dt=60):
    """
    prepare liqiud water content from NCAR research flight (RF) low-rate data
    
    Parameters
    ----------
    RFpath : str
        input path for NCAR RF data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    #%% find all data
    lst = glob.glob(RFpath+'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLON')
        (time,lwcobs,timeunit,lwcunit,lwclongname,cellsize,cellunit)=read_RF_NCAR(filename,'PLWCC')
    
        
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        lwcobs = qc_remove_neg(lwcobs.data)
        
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        lwc = avg_time_1d(time, lwcobs, time_new)
        # lwc = median_time_1d(time, lwcobs, time_new)
        # lwc = median_time_forflight_1d(time, lwcobs, time_new, height)
     
        #%% output data
        outfile = prep_data_path + 'LWC_SOCRATES_' + date + '.nc'
        
        if not os.path.exists(prep_data_path):
            os.makedirs(prep_data_path)
            
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        lwc_o = f.createVariable('LWC', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        lwc_o[:] = lwc
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        lwc_o.units = 'g/m3'
        lwc_o.long_name = 'Liquid water content'
    
        # global attributes
        f.title = "cloud water content from PMS-King probe"
        f.description = 'mean value of each time window'
        f.input_file = filename.split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd(RFpath, prep_data_path, dt=60):
    """
    prepare cloud droplet number and size data from NCAR research flight (RF) low-rate data
    
    Parameters
    ----------
    RFpath : str
        input path for NCAR RF data
    prep_data_path : str
        output path
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.
    
    """    
    #%% find all data
    lst = glob.glob(RFpath+'RF*.PNI.nc')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = filename.split('.')
        date = fname[-4]
        print(date)
        
        (time,height,timeunit,hunit,hlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGALT')
        (time,lat,timeunit,latunit,latlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLAT')
        (time,lon,timeunit,lonunit,lonlongname,cellsize,cellunit)=read_RF_NCAR(filename,'GGLON')
        (time,c1dc,timeunit,c1dcunit,c1dclongname,c1dcsize,cellunit)=read_RF_NCAR(filename,'C1DC_RWOI')
        (time,cdp,timeunit,cdpunit,cdplongname,cdpsize,cellunit)=read_RF_NCAR(filename,'CCDP_RWIO')
        # 2-DC include all-particle (c2dca) and round-particle (C2DCR). use all-particle here
        (time,c2dc,timeunit,c2dcunit,c2dclongname,c2dcsize,cellunit)=read_RF_NCAR(filename,'C2DCA_RWOI')
        (time,cpip,timeunit,cpipunit,cpiplongname,cpipsize,cellunit)=read_RF_NCAR(filename,'CPIP_RWII')
        # 2-DS include C2DSA_2H and C2DSA_2V with all-particle and round-particle. 
        # C2DSA_2V looks off the other instruments in size distribution. use C2DSA_2H
        (time,c2ds2h,timeunit,c2ds2hunit,c2ds2hlongname,c2ds2hsize,cellunit)=read_RF_NCAR(filename,'C2DSA_2H')
        (time,c2ds2v,timeunit,c2ds2vunit,c2ds2vlongname,c2ds2vsize,cellunit)=read_RF_NCAR(filename,'C2DSA_2V')
    
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        c1dc = qc_remove_neg(c1dc[:,0,:].data)
        cdp = qc_remove_neg(cdp[:,0,:].data)
        c2dc = qc_remove_neg(c2dc[:,0,:].data)
        cpip = qc_remove_neg(cpip[:,0,:].data)
        c2ds2h = qc_remove_neg(c2ds2h[:,0,:].data)
        c2ds2v = qc_remove_neg(c2ds2v[:,0,:].data)
        
        cpipsize = (cpipsize[0:-1]+cpipsize[1:])/2
        
        # change all units to #/L
        cdp = cdp * 1000        #cdp unit is #/cm3
    
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        nd_1dc = avg_time_2d(time, c1dc, time_new)
        nd_cdp = avg_time_2d(time, cdp, time_new)
        nd_2dc = avg_time_2d(time, c2dc, time_new)
        nd_pip = avg_time_2d(time, cpip, time_new)
        nd_2ds_h = avg_time_2d(time, c2ds2h, time_new)
        nd_2ds_v = avg_time_2d(time, c2ds2v, time_new)
    
            
        #%% 
        # import matplotlib.pyplot as plt
        # fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize=(8,9))
        # h1=ax1.contourf(time_new, cdpsize, np.log(nd_cdp.T), np.arange(1,10))
        # h2=ax2.contourf(time_new, c2ds2hsize, np.log(nd_2ds_h.T), np.arange(1,10))
        # h3=ax3.contourf(time_new, c2ds2vsize, np.log(nd_2ds_v.T), np.arange(1,10))
        # h4=ax4.contourf(time_new, c1dcsize, np.log(nd_1dc.T), np.arange(-6,3))
        # h5=ax5.contourf(time_new, c2dcsize, np.log(nd_2dc.T), np.arange(-6,3))
        # h6=ax6.contourf(time_new, cpipsize, np.log(nd_pip.T), np.arange(-6,3))
        # ax1.set_ylabel('cdp')
        # ax2.set_ylabel('2DS-2H')
        # ax3.set_ylabel('2DS-2V')
        # ax4.set_ylabel('1DC')
        # ax5.set_ylabel('2DC')
        # ax6.set_ylabel('PIP')
        # ax1.set_yscale('log')
        # ax2.set_yscale('log')
        # ax3.set_yscale('log')
        # ax4.set_yscale('log')
        # ax5.set_yscale('log')
        # ax6.set_yscale('log')
        # ax1.set_title(date)
        # cax = plt.axes([0.92, 0.55, 0.02, 0.3])
        # fig.colorbar(h1,cax=cax)
        # cax = plt.axes([0.92, 0.15, 0.02, 0.3])
        # fig.colorbar(h6,cax=cax)
        # ax1.set_ylim(1,5000)
        # ax2.set_ylim(1,5000)
        # ax3.set_ylim(1,5000)
        # ax4.set_ylim(1,5000)
        # ax5.set_ylim(1,5000)
        # ax6.set_ylim(1,5000)
        
        # # fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))
        # # h1=ax1.contourf(time, c1dcsize, np.log(c1dc.T), np.arange(-6,3))
        # # h2=ax2.contourf(time_new, c1dcsize, np.log(nd_1dc.T), np.arange(-6,3))
        # # h1=ax1.contourf(time, cdpsize, np.log(cdp.T), np.arange(-6,3))
        # # h2=ax2.contourf(time_new, cdpsize, np.log(nd_cdp.T), np.arange(-6,3))
        # # h1=ax1.contourf(time, c2dcsize, np.log(c2dc.T), np.arange(-6,3))
        # # h2=ax2.contourf(time_new, c2dcsize, np.log(nd_2dc.T), np.arange(-6,3))
        # # h1=ax1.contourf(time_new, c2ds2hsize, np.log(nd_2ds_h.T))
        # # h2=ax2.contourf(time_new, c2ds2vsize, np.log(nd_2ds_v.T))
        # # ax1.set_yscale('log')
        # # ax2.set_yscale('log')
        # # ax1.set_title(date)
        # # cax = plt.axes([0.92, 0.2, 0.02, 0.6])
        # # fig.colorbar(h1,cax=cax)
        # # ax1.set_ylim(1,5000)
        # # ax2.set_ylim(1,5000)
        
        # dp_1dc = c1dcsize[1:]-c1dcsize[0:-1]
        # dp_2dc = c2dcsize[1:]-c2dcsize[0:-1]
        # dp_pip = cpipsize[1:]-cpipsize[0:-1]
        # dp_2ds = c2ds2hsize[1:]-c2ds2hsize[0:-1]
        # dp_cdp = cdpsize[1:]-cdpsize[0:-1]
        # dp_cdp = np.hstack((1, dp_cdp))
        
        # fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
        # ax1.plot(c1dcsize, np.nanmean(c1dc,0)/dp_1dc[0],'r')
        # ax1.plot(cdpsize, np.nanmean(cdp,0)/dp_cdp,'m')
        # ax1.plot(c2dcsize, np.nanmean(c2dc,0)/dp_2dc[0],'b')
        # ax1.plot(cpipsize, np.nanmean(cpip,0)/dp_pip[0],'g')
        # ax1.plot(c2ds2hsize, np.nanmean(c2ds2h,0)/dp_2ds[0],'k')
        # ax1.plot(c2ds2vsize, np.nanmean(c2ds2v,0)/dp_2ds[0],'gray')
        # ax1.set_xscale('log')
        # ax1.set_yscale('log')
        # ax2.plot(c1dcsize, np.nanmean(nd_1dc,0)/dp_1dc[0],'r',label='C1DC_RWOI')
        # ax2.plot(cdpsize, np.nanmean(nd_cdp,0)/dp_cdp,'m',label='CCDP_RWIO')
        # ax2.plot(c2dcsize, np.nanmean(nd_2dc,0)/dp_2dc[0],'b',label='C2DCA_RWOI')
        # ax2.plot(cpipsize, np.nanmean(nd_pip,0)/dp_pip[0],'g',label='CPIP_RWII')
        # ax2.plot(c2ds2hsize, np.nanmean(nd_2ds_h,0)/dp_2ds[0],'k',label='C2DSA_2H')
        # ax2.plot(c2ds2vsize, np.nanmean(nd_2ds_v,0)/dp_2ds[0],'gray',label='C2DSA_2V')
        # ax2.set_xscale('log')
        # ax2.set_yscale('log')
        # ax2.legend()
        # ax2.set_xlabel('size(um)')
        # ax2.set_ylabel('dNd/Dp (#/L)')
        # ax2.set_title(date)
        # e
    
        #%% output data
        outfile = prep_data_path + 'Nd_size_SOCRATES_' + date + '.nc'
        
        if not os.path.exists(prep_data_path):
            os.makedirs(prep_data_path)
            
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size_1dc', len(c1dcsize)) 
        f.createDimension('size_2dc', len(c2dcsize)) 
        f.createDimension('size_2ds', len(c2ds2hsize)) 
        f.createDimension('size_cdp', len(cdpsize)) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size1_o = f.createVariable("size_1dc", "f8", ("size_1dc", ))
        size2_o = f.createVariable("size_2dc", "f8", ("size_2dc", ))
        size2h_o = f.createVariable("size_2ds", "f8", ("size_2ds", ))
        sizec_o = f.createVariable("size_cdp", "f8", ("size_cdp", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        nd1_o = f.createVariable('Nd_1dc', 'f8', ("time", "size_1dc"))
        nd2_o = f.createVariable('Nd_2dc', 'f8', ("time", "size_2dc"))
        nd2h_o = f.createVariable('Nd_2ds', 'f8', ("time", "size_2ds"))
        ndc_o = f.createVariable('Nd_cdp', 'f8', ("time", "size_cdp"))
        
        # write data
        time_o[:] = time_new
        size1_o[:] = c1dcsize
        size2_o[:] = c2dcsize
        size2h_o[:] = c2ds2hsize
        sizec_o[:] = cdpsize
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        nd1_o[:, :] = nd_1dc
        nd2_o[:, :] = nd_2dc
        nd2h_o[:, :] = nd_2ds_h
        ndc_o[:, :] = nd_cdp
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        size1_o.units = 'um'
        size1_o.long_name = 'upper bound of size bin for 1DC'
        size2_o.units = 'um'
        size2_o.long_name = 'upper bound of size bin for 2DC'
        size2h_o.units = 'um'
        size2h_o.long_name = 'upper bound of size bin for 2DS'
        sizec_o.units = 'um'
        sizec_o.long_name = 'upper bound of size bin for CDP'
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        nd1_o.units = '#/L'
        nd1_o.long_name = 'cloud droplet number concentration in each bin from 1DC'
        nd2_o.units = '#/L'
        nd2_o.long_name = 'cloud droplet number concentration in each bin from 2DC'
        nd2h_o.units = '#/L'
        nd2h_o.long_name = 'cloud droplet number concentration in each bin from 2DS'
        ndc_o.units = '#/L'
        ndc_o.long_name = 'cloud droplet number concentration in each bin from CDP'
    
        
        # global attributes
        f.title = "cloud droplet size distribution from CDP, 2-DS, 2-DC or 1-DC, in ambient condition"
        f.description = 'average of each time window'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()