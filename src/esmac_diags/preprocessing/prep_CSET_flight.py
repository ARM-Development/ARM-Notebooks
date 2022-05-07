"""
prepare aircraft data from ACEENA
options of average data into coarser resolution
"""

import glob
import os
import re
import numpy as np
import time as ttt
from netCDF4 import Dataset
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, avg_time_2d, \
                    median_time_1d, median_time_2d
from esmac_diags.subroutines.read_aircraft import read_RF_NCAR
from esmac_diags.subroutines.quality_control import qc_cn_max, qc_remove_neg, \
                     qc_mask_cloudflag, qc_uhsas_RF_NCAR
from esmac_diags.subroutines.specific_data_treatment import lwc2cflag


# RFpath = '../../../data/CSET/obs/aircraft/aircraft_lowrate/'
# prep_data_path = 'C:/Users/tang357/Downloads/CSET/'

# dt=60

# if not os.path.exists(prep_data_path):
#     os.makedirs(prep_data_path)
    
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
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
        
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
        (time,uhsas,timeunit,dataunit,long_name,sizeh,cellunit)=read_RF_NCAR(filename,'CUHSAS_RWOOU')
        (time,pcasp,timeunit,dataunit,long_name,sizep,cellunit)=read_RF_NCAR(filename,'CS200_RWOOP')
        # calculate cloud flag based on LWC
        cflag=lwc2cflag(lwc,lwcunit)
        sizeh=sizeh*1000.
        sizel = np.hstack((2*sizeh[0]-sizeh[1],  sizeh[0:-1]))
        size = (sizeh+sizel)/2.
        sizep=sizep*1000.
        sizep2 = np.hstack((2*sizep[0]-sizep[1],  sizep[0:-1]))
        size_pcasp = (sizep+sizep2)/2.
    
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        uhsas = qc_remove_neg(uhsas[:,0,:].data)
        uhsas = qc_mask_cloudflag(uhsas, cflag)
        uhsas = qc_uhsas_RF_NCAR(uhsas)
        # PCASP for CSET are all missing
        pcasp = qc_remove_neg(pcasp[:,0,:].data)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        # uhsas1 = avg_time_2d(time, uhsas, time_new)
        uhsas1 = median_time_2d(time, uhsas, time_new)
        # pcasp1 = median_time_2d(time, pcasp, time_new)
        
        #%% output data
        outfile = prep_data_path + 'UHSASsize_CSET_' + date + '.nc'
        
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
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
        
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
        (time,uhsas100,timeunit,uhsas100unit,uhsas100longname,cellsize,cellunit)=read_RF_NCAR(filename,'CONCU100_RWOOU')
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
        outfile = prep_data_path + 'CN_CSET_' + date + '.nc'
        
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
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
        
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
        outfile = prep_data_path + 'LWC_CSET_' + date + '.nc'
        
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
    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
        
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
        (time,c1dc,timeunit,c1dcunit,c1dclongname,c1dcsize,cellunit)=read_RF_NCAR(filename,'C1DC_LWOO')
        (time,cdp,timeunit,cdpunit,cdplongname,cdpsize,cellunit)=read_RF_NCAR(filename,'CCDP_LWOI')
        # 2-DC include all-particle (c2dca) and round-particle (C2DCR). use all-particle here
        (time,c2dca,timeunit,c2dcaunit,c2dcalongname,c2dcasize,cellunit)=read_RF_NCAR(filename,'C2DCA_LWOO')
        (time,c2dcr,timeunit,c2dcrunit,c2dcrlongname,c2dcrsize,cellunit)=read_RF_NCAR(filename,'C2DCR_LWOO')
    
        
        time = time.filled(np.nan)
        height = height.filled(np.nan)
        lon = lon.filled(np.nan)
        lat = lat.filled(np.nan)
        
        # quality controls
        c1dc = qc_remove_neg(c1dc[:,0,:].data)
        cdp = qc_remove_neg(cdp[:,0,:].data)
        c2dca = qc_remove_neg(c2dca[:,0,:].data)
        c2dcr = qc_remove_neg(c2dcr[:,0,:].data)
        
        # change all units to #/L
        cdp = cdp * 1000        #cdp unit is #/cm3
    
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        nd_1dc = avg_time_2d(time, c1dc, time_new)
        nd_cdp = avg_time_2d(time, cdp, time_new)
        nd_2dca = avg_time_2d(time, c2dca, time_new)
        nd_2dcr = avg_time_2d(time, c2dcr, time_new)
    
        #%% output data
        outfile = prep_data_path + 'Nd_size_CSET_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size_1dc', len(c1dcsize)) 
        f.createDimension('size_2dc', len(c2dcasize)) 
        f.createDimension('size_cdp', len(cdpsize)) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size1_o = f.createVariable("size_1dc", "f8", ("size_1dc", ))
        size2_o = f.createVariable("size_2dc", "f8", ("size_2dc", ))
        sizec_o = f.createVariable("size_cdp", "f8", ("size_cdp", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        nd1_o = f.createVariable('Nd_1dc', 'f8', ("time", "size_1dc"))
        nd2_o = f.createVariable('Nd_2dc', 'f8', ("time", "size_2dc"))
        ndc_o = f.createVariable('Nd_cdp', 'f8', ("time", "size_cdp"))
        
        # write data
        time_o[:] = time_new
        size1_o[:] = c1dcsize
        size2_o[:] = c2dcasize
        sizec_o[:] = cdpsize
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        nd1_o[:, :] = nd_1dc
        nd2_o[:, :] = nd_2dca
        ndc_o[:, :] = nd_cdp
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        size1_o.units = 'um'
        size1_o.long_name = 'upper bound of size bin for 1DC'
        size2_o.units = 'um'
        size2_o.long_name = 'upper bound of size bin for 2DC'
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
        ndc_o.units = '#/L'
        ndc_o.long_name = 'cloud droplet number concentration in each bin from CDP'
    
        
        # global attributes
        f.title = "cloud droplet size distribution from CDP, 2-DC or 1-DC, in ambient condition"
        f.description = 'average of each time window'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
     