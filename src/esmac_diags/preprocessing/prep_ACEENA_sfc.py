"""
prepare surface data from ACEENA
options of output data into coarser resolution
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d,\
                avg_time_1d, avg_time_2d
from esmac_diags.subroutines.read_ARMdata import read_uhsas
from esmac_diags.subroutines.quality_control import qc_remove_neg, qc_mask_qcflag, \
                qc_mask_qcflag_cpc
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_ARM


# #%% 
# import matplotlib.pyplot as plt
# fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))
# ax1.plot(time, reff)
# ax1.plot(time_new, reff_new, 'r.')
# ax2.plot(time, reff)
# ax2.plot(time_new, reff_new, 'r.')
# ax1.set_xlim(17333, 17372)
# ax2.set_xlim(17550, 17595)
# # ax1.set_ylim(0, 2e4)
# # e
    


# acsmpath = '../../../data/ACEENA/obs/surface/arm_acsm/'
# armbepath = '../../../data/ACEENA/obs/profile/armbe/'
# arsclpath = '../../../data/ACEENA/obs/profile/arscl/'
# ccnpath = '../../../data/ACEENA/obs/surface/arm_aosccn1/'
# cpcpath = '../../../data/ACEENA/obs/surface/arm_cpcf/'
# uhsaspath = '../../../data/ACEENA/obs/surface/arm_uhsas/'
# mfrsrpath = '../../../data/ACEENA/obs/surface/arm_mfrsr/'
# Nd_WUpath = '../../../data/ACEENA/obs/surface/Wu_etal_retrieval/'
# ndroppath = '../../../data/ACEENA/obs/surface/enandrop/'
# predatapath = 'C:/Users/tang357/Downloads/ACEENA/'
# dt=3600

# height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
#                 1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
#                 3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
#                 10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ACSM(acsmpath, predatapath, dt=3600):
    """
    prepare acsm data

    Parameters
    ----------
    acsmpath : str
        input datapath
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                           
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(acsmpath, 'enaaosacsmC1.b2*.nc'))
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    org = obsdata['total_organics'].load()
    qc_org = obsdata['qc_total_organics'].load()
    so4 = obsdata['sulfate'].load()
    qc_so4 = obsdata['qc_sulfate'].load()
    nh4 = obsdata['ammonium'].load()
    qc_nh4 = obsdata['qc_ammonium'].load()
    no3 = obsdata['nitrate'].load()
    qc_no3 = obsdata['qc_nitrate'].load()
    chl = obsdata['chloride'].load()
    qc_chl = obsdata['qc_chloride'].load()
    obsdata.close()
    
    # quality controls
    org = qc_mask_qcflag(org,qc_org)
    no3 = qc_mask_qcflag(no3,qc_no3)
    so4 = qc_mask_qcflag(so4,qc_so4)
    nh4 = qc_mask_qcflag(nh4,qc_nh4)
    chl = qc_mask_qcflag(chl,qc_chl)
    org = qc_remove_neg(org)
    no3 = qc_remove_neg(no3)
    so4 = qc_remove_neg(so4)
    nh4 = qc_remove_neg(nh4)
    chl = qc_remove_neg(chl)
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    org_new = median_time_1d(time, org, time_new)
    no3_new = median_time_1d(time, no3, time_new)
    so4_new = median_time_1d(time, so4, time_new)
    nh4_new = median_time_1d(time, nh4, time_new)
    chl_new = median_time_1d(time, chl, time_new)
    
    
    #%% output file
    outfile = predatapath + 'sfc_ACSM_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'org': (['time'], np.float32(org_new)),
                    'so4': (['time'], np.float32(so4_new)),
                    'nh4': (['time'], np.float32(nh4_new)),
                    'no3': (['time'], np.float32(no3_new)),
                    'chl': (['time'], np.float32(chl_new)),
                    },
                      coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['org'].attrs["long_name"] = org.long_name
    ds['org'].attrs["units"] = org.units
    ds['so4'].attrs["long_name"] = so4.long_name
    ds['so4'].attrs["units"] = so4.units
    ds['nh4'].attrs["long_name"] = nh4.long_name
    ds['nh4'].attrs["units"] = nh4.units
    ds['no3'].attrs["long_name"] = no3.long_name
    ds['no3'].attrs["units"] = no3.units
    ds['chl'].attrs["long_name"] = chl.long_name
    ds['chl'].attrs["units"] = chl.units
    
    ds.attrs["title"] = 'Aerosol composition from surface ACSM'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_ccn(ccnpath, predatapath, dt=3600):
    """
    prepare surface CCN data. 
    two IOPs are different .dat files, save them separately

    Parameters
    ----------
    ccnpath : char
        input datapath of CCN data
    predatapath : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    lst = glob.glob(os.path.join(ccnpath, '*.nc'))
    ccndata = xr.open_mfdataset(lst, combine='by_coords')
    time = ccndata['time'].load()
    ccn = ccndata['N_CCN'].load()
    qc_ccn = ccndata['qc_N_CCN'].load()
    ss = ccndata['supersaturation_calculated'].load()
    ccndata.close()
    
    ccn = qc_mask_qcflag(ccn, qc_ccn)
    
    # 0.1%
    idx = np.logical_and(ss>0.05, ss<0.15)
    ccn1 = ccn[idx]
    time1 = time[idx]
    ss1 = ss[idx]
    # 0.2%
    idx = np.logical_and(ss>0.15, ss<0.25)
    ccn2 = ccn[idx]
    time2 = time[idx]
    ss2 = ss[idx]
    # 0.5%
    idx = np.logical_and(ss>0.45, ss<0.55)
    ccn5 = ccn[idx]
    time5 = time[idx]
    ss5 = ss[idx]
    # 0.6%
    idx = np.logical_and(ss>0.55, ss<0.65)
    ccn6 = ccn[idx]
    time6 = time[idx]
    ss6 = ss[idx]
          
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    ccn1_new = median_time_1d(time1, ccn1, time_new)
    ss1_i = median_time_1d(time1, ss1, time_new)
    ccn2_new = median_time_1d(time2, ccn2, time_new)
    ss2_i = median_time_1d(time2, ss2, time_new)
    ccn5_new = median_time_1d(time5, ccn5, time_new)
    ss5_i = median_time_1d(time5, ss5, time_new)
    ccn6_new = median_time_1d(time6, ccn6, time_new)
    ss6_i = median_time_1d(time6, ss6, time_new)
        
    
    #%% output file
    outfile = predatapath + 'sfc_CCN_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'CCN1': ('time', np.float32(ccn1_new)),
                     'CCN2': ('time', np.float32(ccn2_new)),
                     'CCN5': ('time', np.float32(ccn5_new)),
                     'CCN6': ('time', np.float32(ccn6_new)),
                     'ss1': ('time', np.float32(ss1_i)),
                     'ss2': ('time', np.float32(ss2_i)),
                     'ss5': ('time', np.float32(ss5_i)),
                     'ss6': ('time', np.float32(ss6_i)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['CCN1'].attrs["long_name"] = "0.1% Cloud Condensation Nuclei - measured"
    ds['CCN1'].attrs["units"] = "cm-3"
    ds['CCN1'].attrs["description"] = "ARM-measured CCN targetted to 0.1% SS. see SS1 for actual measured SS"
    ds['ss1'].attrs["long_name"] = "Actual Supersaturation targetted to 0.1%"
    ds['ss1'].attrs["units"] = "%"
    ds['ss1'].attrs["description"] = "measured SS that is closest to 0.1%. ccn1_m is measured at this SS"
    ds['CCN2'].attrs["long_name"] = "0.2% Cloud Condensation Nuclei"
    ds['CCN2'].attrs["units"] = "cm-3"
    ds['CCN2'].attrs["description"] = "ARM-measured CCN targetted to 0.2% SS. see SS2 for actual measured SS"
    ds['ss2'].attrs["long_name"] = "Actual Supersaturation targetted to 0.2%"
    ds['ss2'].attrs["units"] = "%"
    ds['ss2'].attrs["description"] = "measured SS that is closest to 0.2%. ccn2_m is measured at this SS"
    ds['CCN5'].attrs["long_name"] = "0.5% Cloud Condensation Nuclei"
    ds['CCN5'].attrs["units"] = "cm-3"
    ds['CCN5'].attrs["description"] = "ARM-measured CCN targetted to 0.5% SS. see SS5 for actual measured SS"
    ds['ss5'].attrs["long_name"] = "Actual Supersaturation targetted to 0.5%"
    ds['ss5'].attrs["units"] = "%"
    ds['ss5'].attrs["description"] = "measured SS that is closest to 0.5%. ccn5_m is measured at this SS"
    ds['CCN6'].attrs["long_name"] = "0.6% Cloud Condensation Nuclei"
    ds['CCN6'].attrs["units"] = "cm-3"
    ds['CCN6'].attrs["description"] = "ARM-measured CCN targetted to 0.6% SS. see SS6 for actual measured SS"
    ds['ss6'].attrs["long_name"] = "Actual Supersaturation targetted to 0.6%"
    ds['ss6'].attrs["units"] = "%"
    ds['ss6'].attrs["description"] = "measured SS that is closest to 0.6%. ccn6_m is measured at this SS"
    
    ds.attrs["title"] = 'Surface CCN number concentration'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CPC(cpcpath,  predatapath, dt=3600):
    """
    prepare CPC data

    Parameters
    ----------
    cpcpath : str
        input datapath for CPC (10nm)
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst1 = glob.glob(os.path.join(cpcpath, '*.nc'))
    lst1.sort()
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    time10 = obsdata['time']
    cpc10 = obsdata['concentration'].load()
    qc_cpc10 = obsdata['qc_concentration'].load()
    obsdata.close()
    
    # quality controls
    cpc10 = qc_mask_qcflag_cpc(cpc10, qc_cpc10.astype(int).data)
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time10[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time10[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    cpc10_new = median_time_1d(time10, cpc10, time_new)
    
    #%% output file
    outfile = predatapath + 'sfc_CPC_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cpc10': (['time'], np.float32(cpc10_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cpc10'].attrs["long_name"] = 'CPC aerosol number concentration (>10nm)'
    ds['cpc10'].attrs["units"] = '1/cm3'
    
    ds.attrs["title"] = 'Aerosol number concentration from surface CPC'
    ds.attrs["inputfile_sample"] = lst1[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_UHSAS(uhsaspath, predatapath, dt=3600):
    """
    prepare UHSAS data

    Parameters
    ----------
    uhsaspath : str
        input datapath for UHSAS
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(uhsaspath, '*.nc'))
    lst.sort()
    
    # UHSAS@ACEENA changes NETCDF filehead in 2017-06-20. 
    (time1, dmin, dmax, uhsas, timeunit, dataunit, long_name) = read_uhsas(lst[0])
    time = np.array([np.datetime64(timeunit[14:24]) + np.timedelta64(int(x),'s') for x in time1])
    for filename in lst[1:]:
        (time1, dmin2, dmax2, uhsas2, timeunit, dataunit, long_name) = read_uhsas(filename)
        time2 = np.array([np.datetime64(timeunit[14:24]) + np.timedelta64(int(x),'s') for x in time1])
        time = np.hstack((time, time2))
        uhsas = np.vstack((uhsas, uhsas2))
    # quality controls
    uhsas = qc_remove_neg(uhsas)
    
    size = (dmin+dmax)/2
    idx100 = dmin>=100
    uhsas100 = np.nansum(uhsas[:,idx100], 1)
    uhsas100[uhsas100==0] = np.nan
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    uhsas_new = median_time_2d(time, uhsas, time_new)
    uhsas100_new = median_time_1d(time, uhsas100, time_new)
        
    
    #%% output file
    outfile = predatapath + 'sfc_UHSAS_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'size_low': (['size'], dmin),
                    'size_high': (['size'], dmax),
                    'uhsas_all': (['time', 'size'], uhsas_new),
                    'uhsas100': (['time'], uhsas100_new),
                    },
                      coords={'time': ('time', time_new), 'size': ('size', size)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['size_low'].attrs["long_name"] = "lower bound of size bin"
    ds['size_low'].attrs["units"] = "nm"
    ds['size_high'].attrs["long_name"] = "upper bound of size bin"
    ds['size_high'].attrs["units"] = "nm"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['uhsas_all'].attrs["long_name"] = 'aerosol number size distribution'
    ds['uhsas_all'].attrs["units"] = '1/cm3'
    ds['uhsas100'].attrs["long_name"] = 'aerosol number concentration for size >100nm'
    ds['uhsas100'].attrs["units"] = '1/cm3'
    
    ds.attrs["title"] = 'Aerosol number concentration and size distribution from UHSAS'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_cloud_2d(armbepath, predatapath, height_out, dt=3600):
    """
    prepare cloud fraction data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    height_out : numpy array
        vertical dimension of output data. will average the original ARSCL resolution to it.
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldrad*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    height = obsdata['height'].load()
    cloud = obsdata['cld_frac'].load()
    qc_cloud = obsdata['qc_cld_frac'].load()
    obsdata.close()    
    
    #%% re-shape the data into coarser resolution
    # startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    cloud_i = np.full((len(time_new),len(height)), np.nan)
    for kk in range(len(height)):
        # quality controls. For ARMBE cloud fraction, remove data with <30% valid points within 1-hr window 
        cl = cloud[:,kk]
        cl[qc_cloud[:,kk]>=2] = np.nan
        # interpolate into standard time
        cloud_i[:,kk] = np.interp(time_new, time, cl)
        
    cloud_o = avg_time_2d(height,cloud_i.T,height_out).T
        
    
    #%% output file
    outfile = predatapath + 'cloud_2d_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cloud': (['time','height'], np.float32(cloud_o))
                    },
                     coords={'time': ('time', time_new), 'height':('height',np.float32(height_out))})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cloud'].attrs["long_name"] = "Cloud Fraction"
    ds['cloud'].attrs["units"] = "%"
    ds['cloud'].attrs["description"] = "Cloud Fraction based on radar and MPL"
    
    ds.attrs["title"] = 'ARMBE hourly 2-d cloud fraction data derived from ARSCL data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'interpolated into each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_cloudheight_ARSCL(arsclpath, predatapath, dt=3600):
    """
    prepare cloud base and top height data at ARM sites from ARSCL
    include multi-layer clouds
    
    Parameters
    ----------
    arsclpath : char
        input datapath.  
    predatapath : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
            
    #%% read in data
    lst1 = glob.glob(os.path.join(arsclpath, 'enaarsclkazrbnd1kolliasC1.c0*.nc'))
    lst1.sort()
    arscldata = xr.open_mfdataset(lst1, combine='by_coords')
    arscltime = arscldata['time'].load()
    cbh = arscldata['cloud_base_best_estimate'].load()
    cbhs = arscldata['cloud_layer_base_height'].load()
    cths = arscldata['cloud_layer_top_height'].load()
    arscldata.close()
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # data treatments
    cbh = qc_remove_neg(cbh, remove_zero='True')
    for ll in range(cbhs.shape[1]):
        cbhs[:,ll] = qc_remove_neg(cbhs[:,ll], remove_zero='True')
        cths[:,ll] = qc_remove_neg(cths[:,ll], remove_zero='True')
    cth = np.nanmax(cths,axis=1)  # cloud top height for all clouds
        
    #%% re-shape the data into coarser resolution
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    cbh_new = avg_time_1d(arscltime, cbh, time_new)
    cth_new = avg_time_1d(arscltime, cth, time_new)
    cbhs_new = avg_time_2d(arscltime, cbhs, time_new)
    cths_new = avg_time_2d(arscltime, cths, time_new)
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # output file
    outfile = predatapath + 'cloudheight_ARSCL_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cbh': ('time', np.float32(cbh_new)),
                     'cth': ('time', np.float32(cth_new)),
                     'cbhs': (['time', 'layer'], np.float32(cbhs_new)),
                     'cths': (['time', 'layer'], np.float32(cths_new)),
                    },
                     coords={'time': ('time', time_new), 'layer': ('layer', np.arange(cbhs.shape[1]))})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cbh'].attrs["long_name"] = 'cloud base height, best estimate from ceilometer and MPL'
    ds['cbh'].attrs["units"] = "m"
    ds['cth'].attrs["long_name"] = 'cloud top height, the highest level of cths'
    ds['cth'].attrs["units"] = "m"
    ds['cbhs'].attrs["long_name"] = 'Base height of hydrometeor layers for up to 10 layers'
    ds['cbhs'].attrs["units"] = "m"
    ds['cths'].attrs["long_name"] = 'Top height of hydrometeor layers for up to 10 layers'
    ds['cths'].attrs["units"] = "m"
    
    
    ds.attrs["title"] = "cloud heights for up to 10 cloud layers from ARSCL data"
    ds.attrs["description"] = 'average value of each time window'
    ds.attrs["inputfile_sample"] = lst1[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LWP(armbepath, mfrsrpath, predatapath, dt=3600):
    """
    prepare liquid water path
    Although LWP is measured by microwave radiometer (MWR), it is processed in 
    different ways and stored in different datastreams.
    Here we get LWP from two datasets: ARMBE and MFRSR

    Parameters
    ----------
    armbepath : str
        input datapath for ARMBE (MWR data)
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% read in data
    lst1 = glob.glob(os.path.join(armbepath, '*armbecldradC1*.nc'))
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    
    time1 = obsdata['time']
    lwp = obsdata['lwp']
    qc_lwp = obsdata['qc_lwp']
    obsdata.close()
    
    lwp.load()
    qc_lwp.load()
    
    # quality controls. For ARMBE lwp, remove data with <30% valid points within 1-hr window 
    lwp[qc_lwp>=2] = np.nan
    
    #%% read in MFRSR LWP for comparison
    lst2 = glob.glob(os.path.join(mfrsrpath, '*.cdf'))
    lst2.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst2[0])
    time2 = mfrsrdata['time']
    lwp2 = mfrsrdata['lwp']
    qc_lwp2 = mfrsrdata['qc_lwp']
    mfrsrdata.close()
    for file in lst2[1:]:
        mfrsrdata = xr.open_dataset(file)
        time2 = xr.concat([time2, mfrsrdata['time']], dim="time")
        lwp2 = xr.concat([lwp2, mfrsrdata['lwp']], dim="time")
        qc_lwp2 = xr.concat([qc_lwp2, mfrsrdata['qc_lwp']], dim="time")
        mfrsrdata.close()
        
    lwp2.load()
    qc_lwp2.load()
    lwp2 = qc_mask_qcflag(lwp2, qc_lwp2)
    lwp2 = lwp2*1000. # change unit from kg/m2 (mm) to g/m2
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lwp_new = avg_time_1d(time1, lwp, time_new)
    lwp2_new = avg_time_1d(time2, lwp2, time_new)
    
    #%% output file
    outfile = predatapath + 'LWP_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwp_armbe': ('time', np.float32(lwp_new)),
                    'lwp_mfrsr': ('time', np.float32(lwp2_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwp_armbe'].attrs["long_name"] = "liquid water path"
    ds['lwp_armbe'].attrs["units"] = "g/m2"
    ds['lwp_armbe'].attrs["description"] = "liquid water path from ARMBE data based on MWR measurements"
    ds['lwp_mfrsr'].attrs["long_name"] = "liquid water path"
    ds['lwp_mfrsr'].attrs["units"] = "g/m2"
    ds['lwp_mfrsr'].attrs["description"] = "liquid water path from MFRSR data based on MWR measurements"
    
    ds.attrs["title"] = 'surface-retrieved cloud liquid water path'
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["description"] = 'mean value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_LTS(armbepath, predatapath, dt=3600):
    """
    prepare lower tropospheric stability (potential temperature difference between 700hPa and surface) from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                   
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbeatmC1*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    pres = obsdata['pressure'].load()
    T = obsdata['temperature_p'].load()
    Ts = obsdata['temperature_sfc'].load()
    ps = obsdata['pressure_sfc'].load()
    obsdata.close()    
    
    # in ARMBEATM, pres[12] is 700hPa. check if not
    p700_idx = 12
    if np.abs(pres[p700_idx]-700)>0.1:
        raise ValueError('this level is not 700hPa, check: '+str(pres[p700_idx]))
    p850_idx = 6
    if np.abs(pres[p850_idx]-850)>0.1:
        raise ValueError('this level is not 850hPa, check: '+str(pres[p850_idx]))
    
    # remove one timestep that temperature is incorrect
    idx = T[:,30]>300
    T[idx.data,:]=np.nan
    
    # only extract valid data when sounding is available
    LTS700_valid  = np.empty((0))
    LTS850_valid  = np.empty((0))
    time700_valid = np.empty(0,dtype='datetime64')
    time850_valid = np.empty(0,dtype='datetime64')

    for tt in range(len(time)):
        if ~np.isnan(T[tt,p700_idx]):
            time700_valid = np.hstack((time700_valid,time[tt]))
            theta_s = Ts[tt].data * (1000/ps[tt].data)**0.286
            theta_700 = T[tt,p700_idx].data * (1000/700)**0.286       
            LTS700_valid = np.append(LTS700_valid, theta_700-theta_s )
        if ~np.isnan(T[tt,p850_idx]):
            time850_valid = np.hstack((time850_valid,time[tt]))
            theta_s = Ts[tt].data * (1000/ps[tt].data)**0.286
            theta_850 = T[tt,p850_idx].data * (1000/850)**0.286       
            LTS850_valid = np.append(LTS850_valid, theta_850-theta_s )
        
    #%% re-shape the data into coarser resolution
    # startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    LTS700_new = avg_time_1d(time700_valid, LTS700_valid, time_new)
    LTS850_new = avg_time_1d(time850_valid, LTS850_valid, time_new)
        
    #%% output file
    outfile = predatapath + 'LTS_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'LTS700': (['time'], np.float32(LTS700_new)),
                    'LTS850': (['time'], np.float32(LTS850_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['LTS700'].attrs["long_name"] = "Lower Tropospheric Stability"
    ds['LTS700'].attrs["units"] = "K"
    ds['LTS700'].attrs["description"] = "potential temperature difference between 700hPa and surface, " + \
        "calculated from radiosonde  and surface data"
    ds['LTS850'].attrs["long_name"] = "Lower Tropospheric Stability"
    ds['LTS850'].attrs["units"] = "K"
    ds['LTS850'].attrs["description"] = "potential temperature difference between 850hPa and surface, " + \
        "calculated from radiosonde  and surface data"
    
    ds.attrs["title"] = 'Lower Tropospheric Stability from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_cod(mfrsrpath,  predatapath, dt=3600):
    """
    prepare cloud optical depth data from MFRSR

    Parameters
    ----------
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(mfrsrpath, '*.cdf'))
    lst.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst[0])
    time = mfrsrdata['time']
    cod = mfrsrdata['optical_depth_instantaneous']
    qc_cod = mfrsrdata['qc_optical_depth_instantaneous']
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
        qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
        mfrsrdata.close()
    
    # quality controls
    cod.load()
    qc_cod.load()
    cod = qc_mask_qcflag(cod, qc_cod)
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    cod_new = avg_time_1d(time, cod, time_new)
    
    #%% output file
    outfile = predatapath + 'cod_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cod': (['time'], np.float32(cod_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cod'].attrs["long_name"] = 'cloud optical depth'
    ds['cod'].attrs["units"] = cod.units
    
    ds.attrs["title"] = 'cloud optical depth from surface MFRSR'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mfrsr_Reff(mfrsrpath,  predatapath, dt=3600):
    """
    prepare cloud effective radius data from MFRSR

    Parameters
    ----------
    mfrsrpath : str
        input datapath for MFRSR
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(mfrsrpath, '*.cdf'))
    lst.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst[0])
    time = mfrsrdata['time']
    reff = mfrsrdata['effective_radius_instantaneous']
    qc_reff = mfrsrdata['qc_effective_radius_instantaneous']
    mfrsrdata.close()
    for file in lst[1:]:
        mfrsrdata = xr.open_dataset(file)
        time = xr.concat([time, mfrsrdata['time']], dim="time")
        reff = xr.concat([reff, mfrsrdata['effective_radius_instantaneous']], dim="time")
        qc_reff = xr.concat([qc_reff, mfrsrdata['qc_effective_radius_instantaneous']], dim="time")
        mfrsrdata.close()
    
    # quality controls
    reff.load()
    qc_reff.load()
    reff = qc_mask_qcflag(reff, qc_reff)
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    reff_new = median_time_1d(time, reff, time_new)
    
    #%% output file
    outfile = predatapath + 'reff_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'reff': (['time'], np.float32(reff_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['reff'].attrs["long_name"] = 'cloud effective radius'
    ds['reff'].attrs["units"] = reff.units
    
    ds.attrs["title"] = 'cloud effective radius from surface MFRSR'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_precip(armbepath, predatapath, dt=3600):
    """
    prepare surface precipitation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                   
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbeatmC1*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    precip = obsdata['precip_rate_sfc'].load()
    obsdata.close()    
                
    #%% re-shape the data into coarser resolution
    # startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    precip_new = avg_time_1d(time, precip, time_new)
        
    #%% output file
    outfile = predatapath + 'precip_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'precip': (['time'], np.float32(precip_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['precip'].attrs["long_name"] = "Surface Precipitation Rate"
    ds['precip'].attrs["units"] = "mm/hr"
    
    ds.attrs["title"] = 'surface precipitation data from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'precipitation at ENA is measured by ARM Present Weather Detector (PWD). mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_radiation(armbepath, predatapath, dt=3600):
    """
    prepare surface radiation data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                   
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldradC1*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    lwdn = obsdata['lwdn'].load()
    lwup = obsdata['lwup'].load()
    swdn = obsdata['swdn'].load()
    swup = obsdata['swup'].load()
    obsdata.close()    
                
    #%% re-shape the data into coarser resolution
    # startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    lwdn_new = avg_time_1d(time, lwdn, time_new)
    swdn_new = avg_time_1d(time, swdn, time_new)
    lwup_new = avg_time_1d(time, lwup, time_new)
    swup_new = avg_time_1d(time, swup, time_new)
        
    #%% output file
    outfile = predatapath + 'sfc_radiation_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwdn': (['time'], np.float32(lwdn_new)),
                    'swdn': (['time'], np.float32(swdn_new)),
                    'lwup': (['time'], np.float32(lwup_new)),
                    'swup': (['time'], np.float32(swup_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwdn'].attrs["long_name"] = "Surface downward longwave flux"
    ds['lwdn'].attrs["units"] = "W/m2"
    ds['swdn'].attrs["long_name"] = "Surface downward shortwave flux"
    ds['swdn'].attrs["units"] = "W/m2"
    ds['lwup'].attrs["long_name"] = "Surface upward longwave flux"
    ds['lwup'].attrs["units"] = "W/m2"
    ds['swup'].attrs["long_name"] = "Surface upward shortwave flux"
    ds['swup'].attrs["units"] = "W/m2"
    
    ds.attrs["title"] = 'surface radiative flux data from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["description"] = 'mean of each time window'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_totcld(armbepath, predatapath, dt=3600):
    """
    prepare total cloud fraction data from ARMBE

    Parameters
    ----------
    armbepath : str
        input datapath. use hourly-averaged ARMBE data
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                       
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(armbepath, '*armbecldrad*.nc'))
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time']
    cf_arscl = obsdata['tot_cld']
    qc_cf_arscl = obsdata['qc_tot_cld']
    cf_tsi = obsdata['tot_cld_tsi']
    qc_cf_tsi = obsdata['qc_tot_cld_tsi']
    cf_visst = obsdata['cld_tot']
    obsdata.close()
    
    cf_arscl.load()
    cf_tsi.load()
    cf_visst.load()
    qc_cf_arscl.load()
    qc_cf_tsi.load()
    
    # quality controls. For ARMBE cloud fraction, <30% valid points within 1-hr window are flagged and removed
    cf_arscl[qc_cf_arscl>=2] = np.nan
    cf_tsi[qc_cf_tsi>=2] = np.nan
    
    # change unit from 1 to %
    cf_arscl = cf_arscl*100
    cf_tsi = cf_tsi*100
                
    #%% re-shape the data into coarser resolution
    # startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    # enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    # time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    time_new = pd.date_range(start='2017-06-01', end='2018-03-01', freq=str(int(dt))+"s")  # ACEENA time period
    
    cf_arscl_new = avg_time_1d(time, cf_arscl, time_new)
    cf_tsi_new = avg_time_1d(time, cf_tsi, time_new)
    cf_visst_new = avg_time_1d(time, cf_visst, time_new)
    
    #%% output file
    outfile = predatapath + 'totcld_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'tot_cld_arscl': ('time', np.float32(cf_arscl_new)),
                    'tot_cld_tsi': ('time', np.float32(cf_tsi_new)),
                    'tot_cld_visst': ('time', np.float32(cf_visst_new))
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['tot_cld_arscl'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_arscl'].attrs["units"] = "%"
    ds['tot_cld_arscl'].attrs["description"] = "Total cloud fraction based on ARSCL radar/MPL measurements"
    ds['tot_cld_tsi'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_tsi'].attrs["units"] = "%"
    ds['tot_cld_tsi'].attrs["description"] = "Total cloud fraction based on total sky imager, 100 degree FOV"
    ds['tot_cld_visst'].attrs["long_name"] = "cloud_area_fraction"
    ds['tot_cld_visst'].attrs["units"] = "%"
    ds['tot_cld_visst'].attrs["description"] = "Total cloud fraction based on VISST satellite product"
    
    
    ds.attrs["title"] = 'total cloud fraction from ARMBE hourly data'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Ndrop(ndroppath, predatapath, dt=3600):
    """
    prepare cloud deoplet number concentration from Ndrop data
    
    Parameters
    ----------
    ndroppath : char
        input datapath.  
    predatapath : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst = glob.glob(os.path.join(ndroppath, '*ndropmfrsrC1.c1*.nc'))
    lst.sort()
    obsdata = xr.open_mfdataset(lst, combine='by_coords')
    time = obsdata['time'].load()
    nd = obsdata['drop_number_conc'].load()
    qc_nd = obsdata['qc_drop_number_conc'].load()
    obsdata.close()
    
    # quality control
    nd = qc_mask_qcflag(nd,qc_nd)
    nd = nd*1e-6   # m-3 to cm-3
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1].data))[:10]
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, nd, time_new)
    # nd_new = avg_time_1d(time, nd, time_new)
    
    #%% output file
    outfile = predatapath + 'Ndrop_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cdnc': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cdnc'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['cdnc'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_ARMretrieval(mfrsrpath, arsclpath, predatapath, dt=3600):
    """
    prepare cloud deoplet number concentration (Nd) data at ARM sites
    input data is cloud optical depth from MFRSR, LWP from MWR (in the MFRSR data), 
    and cloud top/base height from ARSCL
    using  ARM surface-based retrievals algorithm
    https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf 
    assuming cloud is homogeneous between cloud top and cloud base in the column
    
    Parameters
    ----------
    mfrsrpath : char
        input datapath. 
    arsclpath : char
        input datapath.  
    predatapath : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    lst1 = glob.glob(os.path.join(arsclpath, 'enaarsclkazrbnd1kolliasC1.c0*.nc'))
    lst1.sort()
    arscldata = xr.open_mfdataset(lst1, combine='by_coords')
    arscltime = arscldata['time']
    cbh = arscldata['cloud_base_best_estimate']
    cbhs = arscldata['cloud_layer_base_height']
    cths = arscldata['cloud_layer_top_height']
    arscldata.close()
    
    lst2 = glob.glob(os.path.join(mfrsrpath, '*.cdf'))
    lst2.sort()
    # first data
    mfrsrdata = xr.open_dataset(lst2[0])
    mfrsrtime = mfrsrdata['time']
    lwp = mfrsrdata['lwp']
    qc_lwp = mfrsrdata['qc_lwp']
    cod = mfrsrdata['optical_depth_instantaneous']
    qc_cod = mfrsrdata['qc_optical_depth_instantaneous']
    mfrsrdata.close()
    for file in lst2[1:]:
        mfrsrdata = xr.open_dataset(file)
        mfrsrtime = xr.concat([mfrsrtime, mfrsrdata['time']], dim="time")
        lwp = xr.concat([lwp, mfrsrdata['lwp']], dim="time")
        qc_lwp = xr.concat([qc_lwp, mfrsrdata['qc_lwp']], dim="time")
        cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
        qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
        mfrsrdata.close()
    
    lwp.load()
    qc_lwp.load()
    cod.load()
    qc_cod.load()
    cbh.load()
    cbhs.load()
    cths.load()
    
    # lwp = qc_mask_qcflag(lwp,qc_lwp)  # do not mask LWP since clearsky (LWP=0) is flagged
    cod = qc_mask_qcflag(cod,qc_cod)
    cbh = qc_remove_neg(cbh, remove_zero='True')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # calculate CDNC
    # calculate cloud depth for the lowest cloud layer in ARSCL
    cth = cths[:,0]
    H = cths[:,0] - cbh
    H.load()
    H = qc_remove_neg(H, remove_zero='true')
    H[cbhs[:,1] > 0] = np.nan
    
    # filter some data samples
    H[cth>5000.] = np.nan   # remove deep clouds with cloud top >5km
    lwp[np.logical_or(lwp<0.02, lwp>0.3)] = np.nan
    cod[np.logical_or(cod<2, cod>60)] = np.nan
    
    # calculate CDNC first then average into 1hr
    time = mfrsrtime.data
    H_tmp = np.interp(np.int64(time), np.int64(arscltime), H)
    nd = calc_cdnc_ARM(lwp, cod, H_tmp)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, nd, time_new)
        
    #%% output file
    outfile = predatapath + 'Nd_ARMretrieval_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cdnc': ('time', np.float32(nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cdnc'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['cdnc'].attrs["units"] = "cm-3"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from ARM retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf '
    ds.attrs["inputfile_sample"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_WU(Nd_WUpath, predatapath, dt=3600):
    """
    prepare cloud deoplet number concentration (Nd) and effective radius data for ACEENA
    retrieval algorithm developed by Peng Wu, Xiquan Dong and Baike Xi
    reference: https://doi.org/10.1175/jcli-d-20-0272.1
    
    Parameters
    ----------
    Nd_WUpath : char
        input datapath. 
    predatapath : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
        
    #%% read in data
    time = np.array([],dtype='datetime64[s]')
    ndall = np.array([])
    reall = np.array([])
    lst = glob.glob(os.path.join(Nd_WUpath, 'ENA_micro_sfc_retrieval*.nc'))
    lst.sort()
    for fname in lst:
        f = Dataset(fname,'r')
        hr = f.variables['Time'][:]
        nc = f.variables['NC'][:]
        re0 = f.variables['RC'][:]
        re = np.nanmedian(re0, axis=0)
        day = fname[-11:-7]+'-'+fname[-7:-5]+'-'+fname[-5:-3]
        f.close()
        
        time1 = np.array([np.datetime64(day) + np.timedelta64(int(x*3600),'s') for x in hr])
        time = np.hstack((time, time1))
        ndall = np.hstack((ndall, nc))
        reall = np.hstack((reall, re))
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, ndall, time_new)
    re_new = median_time_1d(time, reall, time_new)
    
    #%% output file
    outfile = predatapath + 'Nd_Reff_Wu_etal_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                     'cdnc': ('time', np.float32(nd_new)),
                     'reff': ('time', np.float32(re_new)),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cdnc'].attrs["long_name"] = "cloud droplet number concentration (layer-mean)"
    ds['cdnc'].attrs["units"] = "cm-3"
    ds['reff'].attrs["long_name"] = "cloud droplet effective radius"
    ds['reff'].attrs["units"] = "um"
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from Peng Wu et al's retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://doi.org/10.1175/jcli-d-20-0272.1'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
