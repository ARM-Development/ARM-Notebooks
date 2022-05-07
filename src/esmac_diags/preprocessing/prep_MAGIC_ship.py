"""
prepare ship data from MAGIC
options of output data into coarser resolution
"""

import glob
import os
import xarray as xr
import pandas as pd
import numpy as np
import time as ttt
from netCDF4 import Dataset
import esmac_diags
# from esmac_diags.subroutines.read_ship import read_marmet
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, median_time_1d, median_time_2d
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag, qc_cn_max, qc_ccn_max

# shipmetpath = '../../../data/MAGIC/obs/ship/magmarinemetM1.b1/'
# mwrpath = '../../../data/MAGIC/obs/ship/magmwrret1liljclouM1.s2/'
# cpcpath = '../../../data/MAGIC/obs/ship/magaoscpcfM1.a1/'
# ccnpath = '../../../data/MAGIC/obs/ship/magaosccn100M1.a1/'
# uhsaspath = '../../../data/MAGIC/obs/ship/magaosuhsasM1.a1/'
# Ndpath = '../../../data/MAGIC/obs/ship/Cloud_Micro_Retrieval/'
# prep_data_path = 'C:/Users/tang357/Downloads/MAGIC/'

# dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN(shipmetpath, ccnpath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    ccnpath : str
        input path of CCN data
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
    
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load().data
    lon = shipdata['lon_mean_gps'].load().data
    qc_lon = shipdata['qc_lon_mean_gps'].load().data
    lat = shipdata['lat_mean_gps'].load().data
    qc_lat = shipdata['qc_lat_mean_gps'].load().data
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    
    #%% read in data
    lst2 = glob.glob(ccnpath+'magaosccn100M1.a1.*')
    
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    ccn = obsdata['N_CCN'].load().data
    SS = obsdata['CCN_ss_set'].load().data
    obsdata.close()
    
    ccn=qc_ccn_max(ccn,SS)
    
    ccn_1s = np.array(ccn)
    ccn_2s = np.array(ccn)
    ccn_3s = np.array(ccn)
    ccn_5s = np.array(ccn)
    ccn_6s = np.array(ccn)
    ccn_1s[np.logical_or(SS<0.05, SS>0.15)] = np.nan
    ccn_2s[np.logical_or(SS<0.15, SS>0.25)] = np.nan
    ccn_3s[np.logical_or(SS<0.25, SS>0.35)] = np.nan
    ccn_5s[np.logical_or(SS<0.45, SS>0.55)] = np.nan
    ccn_6s[np.logical_or(SS<0.55, SS>0.65)] = np.nan
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    ccn1 = median_time_1d(time2, ccn_1s, time_new)
    ccn2 = median_time_1d(time2, ccn_2s, time_new)
    ccn3 = median_time_1d(time2, ccn_3s, time_new)
    ccn5 = median_time_1d(time2, ccn_5s, time_new)
    ccn6 = median_time_1d(time2, ccn_6s, time_new)
    
    #%% 
    # import matplotlib.pyplot as plt
    
    # fig = plt.figure(figsize=(8,6))
    # ax1 = fig.add_subplot(3, 1, 1)
    # ax1.plot(time2, ccn_1s)
    # ax1.plot(time_new, ccn1, color='r', marker='.',linewidth=2)
    # # ax1.set_ylim(.5, .7)
    # ax2 = fig.add_subplot(3, 1, 2)
    # ax2.plot(time2, ccn_6s)
    # ax2.plot(time_new, ccn6, color='r', marker='.',linewidth=2)
    # # ax2.set_ylim(-50, 2000)
    # ax3= fig.add_subplot(3, 1, 3)
    # ax3.plot(time2, ccn_2s)
    # ax3.plot(time_new, ccn2, color='r', marker='.',linewidth=2)
    # # ax2.set_ylim(-50, 2000)
    
    # # fig = plt.figure(figsize=(8,3))
    # # ax1 = fig.add_subplot(1, 1, 1)
    # # h1=ax1.contourf(time_new, size.data, uhsas1.T, np.arange(0,200,20))
    # # fig.colorbar(h1)
    # # ax1.set_yscale('log')
    # # import matplotlib.dates as mdates
    # # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    # r
    
    #%% output file
    outfile = prep_data_path + 'CCN_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'CCN1': (['time'], ccn1),
                    'CCN2': (['time'], ccn2),
                    'CCN3': (['time'], ccn3),
                    'CCN5': (['time'], ccn5),
                    'CCN6': (['time'], ccn6),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CCN1'].attrs["long_name"] = 'CCN concentration for SS between 0.05% and 0.15%'
    ds['CCN1'].attrs["units"] = "1/cm3"
    ds['CCN2'].attrs["long_name"] = 'CCN concentration for SS between 0.15% and 0.25%'
    ds['CCN2'].attrs["units"] = "1/cm3"
    ds['CCN3'].attrs["long_name"] = 'CCN concentration for SS between 0.25% and 0.35%'
    ds['CCN3'].attrs["units"] = "1/cm3"
    ds['CCN5'].attrs["long_name"] = 'CCN concentration for SS between 0.45% and 0.55%'
    ds['CCN5'].attrs["units"] = "1/cm3"
    ds['CCN6'].attrs["long_name"] = 'CCN concentration for SS between 0.55% and 0.65%'
    ds['CCN6'].attrs["units"] = "1/cm3"
    
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    cpcpath : str
        input path of CPC data
    uhsaspath : str
        input path of UHSAS data
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
    
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon_mean_gps'].load()
    qc_lon = shipdata['qc_lon_mean_gps'].load()
    lat = shipdata['lat_mean_gps'].load()
    qc_lat = shipdata['qc_lat_mean_gps'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    
    #%% read in data
    lst1 = glob.glob(cpcpath+'magaoscpcfM1.a1.*')
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    time1 = obsdata['time'].load()
    cpc = obsdata['concentration'].load()
    obsdata.close()
    
    lst2 = glob.glob(uhsaspath+'magaosuhsasM1.a1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    raw_count = obsdata['size_distribution'].load()
    flow_rate = obsdata['sampling_volume'].load()/60.    # cc/min to cc/s
    obsdata.close()
    
    sample_time = 10    # sample interval is 10s
    uhsas=np.full(raw_count.shape, np.nan)
    for bb in range(uhsas.shape[1]):
        uhsas[:, bb] = raw_count[:, bb].data /flow_rate.data /sample_time
    dataunit='1/cm3'
    uhsas = qc_remove_neg(uhsas)
    uhsas100 = np.nansum(uhsas[:, dmin>100], axis=1)
    
    cpc = qc_cn_max(cpc,10)
    uhsas100 = qc_cn_max(uhsas100,100)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    cpc1 = median_time_1d(time1, cpc, time_new)
    uhsas1 = median_time_1d(time2, uhsas100, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CN_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'CPC10': (['time'], cpc1),
                    'UHSAS100': (['time'], uhsas1),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CPC10'].attrs["long_name"] = 'CPC measured aerosol number (size>10nm)'
    ds['CPC10'].attrs["units"] = "1/cm3"
    ds['UHSAS100'].attrs["long_name"] = 'UHSAS measured aerosol number (size>100nm)'
    ds['UHSAS100'].attrs["units"] = "1/cm3"
    
    ds.attrs["input data_example"] = [lst1[0].split('/')[-1], lst2[0].split('/')[-1]]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize(shipmetpath, uhsaspath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    cpcpath : str
        input path of CPC data
    uhsaspath : str
        input path of UHSAS data
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
    
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon_mean_gps'].load()
    qc_lon = shipdata['qc_lon_mean_gps'].load()
    lat = shipdata['lat_mean_gps'].load()
    qc_lat = shipdata['qc_lat_mean_gps'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    
    #%% read in data
    
    lst2 = glob.glob(uhsaspath+'magaosuhsasM1.a1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    raw_count = obsdata['size_distribution'].load()
    flow_rate = obsdata['sampling_volume'].load()/60.    # cc/min to cc/s
    obsdata.close()
    
    sample_time = 10    # sample interval is 10s
    uhsas=np.full(raw_count.shape, np.nan)
    for bb in range(uhsas.shape[1]):
        uhsas[:, bb] = raw_count[:, bb].data /flow_rate.data /sample_time
    # dataunit='1/cm3'
    
    uhsas = qc_remove_neg(uhsas)
    size = (dmin+dmax)/2
    
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    uhsas1 = median_time_2d(time2, uhsas, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CNsize_UHSAS_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'size_low': (['size'], dmin.data),
                    'size_high': (['size'], dmax.data),
                    'size_distribution_uhsas': (['time', 'size'], uhsas1),
                    },
                     coords={'time': ('time', time_new), 'size': ('size', size.data)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['size_low'].attrs["long_name"] = "lower bound of size bin"
    ds['size_low'].attrs["units"] = "nm"
    ds['size_high'].attrs["long_name"] = "upper bound of size bin"
    ds['size_high'].attrs["units"] = "nm"
    ds['size'].attrs["long_name"] = "aerosol size"
    ds['size'].attrs["units"] = "nm"
    ds['size'].attrs["description"] = "middle of bin: 0.5*(dmin+dmax)"
    ds['size_distribution_uhsas'].attrs["long_name"] = 'aerosol number size distribution'
    ds['size_distribution_uhsas'].attrs["units"] = '1/cm3'
    
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_MWR(shipmetpath, mwrpath, prep_data_path, dt=3600):
    """
    prepare LWP retrieved from MWR
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    mwrpath : str
        input path for MWR data
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
    
    lst = glob.glob(shipmetpath+'magmarinemetM1*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon_mean_gps'].load()
    qc_lon = shipdata['qc_lon_mean_gps'].load()
    lat = shipdata['lat_mean_gps'].load()
    qc_lat = shipdata['qc_lat_mean_gps'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in MWR
    lst2 = glob.glob(mwrpath+'magmwrret1liljclouM1.s2.*')
    
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    lwp = obsdata['be_lwp'].load()
    qc_lwp = obsdata['qc_be_lwp'].load()
    obsdata.close()
    
    lwp = qc_mask_qcflag(lwp, qc_lwp)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = avg_time_1d(time, lon, time_new)
    lat1 = avg_time_1d(time, lat, time_new)
    lwp1 = avg_time_1d(time2, lwp, time_new)
    lwp1 = qc_remove_neg(lwp1)
    
    #%% output file
    outfile = prep_data_path + 'LWP_MAGIC.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'lwp': (['time'], lwp1),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['lwp'].attrs["long_name"] = lwp.long_name
    ds['lwp'].attrs["units"] = lwp.units
    
    ds.attrs["input data_example"] = lst2[0].split('/')[-1]
    ds.attrs["description"] = 'average into '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_Nd_Wu_etal(Ndpath, prep_data_path, dt=3600):
    """
    prepare cloud deoplet number concentration (Nd) and effective radius data for MAGIC
    retrieval algorithm developed by Peng Wu and Xiquan Dong 
    reference: https://doi.org/10.1175/jcli-d-20-0272.1 (method)
            https://doi.org/10.1029/2020EA001588 (data analysis)
    
    
    Parameters
    ----------
    Ndpath : char
        input datapath. 
    prep_data_path : char
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """

    if not os.path.exists(prep_data_path):
        os.makedirs(prep_data_path)
        
    #%% read in data
    time = np.array([],dtype='datetime64[s]')
    ndall = np.array([])
    reall = np.array([])
    lst = glob.glob(os.path.join(Ndpath, 'MAG_*.nc'))
    lst.sort()
    for fname in lst:
        f = Dataset(fname,'r')
        hr = f.variables['time'][:]
        nc = f.variables['nc'][:]
        re0 = f.variables['rc'][:]
        re = np.nanmedian(re0, axis=0)
        day = fname[-36:-32]+'-'+fname[-32:-30]+'-'+fname[-30:-28]
        f.close()
        
        time1 = np.array([np.datetime64(day) + np.timedelta64(int(x*3600),'s') for x in hr])
        time = np.hstack((time, time1))
        ndall = np.hstack((ndall, nc))
        reall = np.hstack((reall, re))
    
    ndall = qc_remove_neg(ndall)
    reall = qc_remove_neg(reall)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time[0]))[:10]
    enddate = np.datetime_as_string(np.datetime64(time[-1]))[:10]
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    nd_new = median_time_1d(time, ndall, time_new)
    re_new = median_time_1d(time, reall, time_new)
    
    #%% output file
    outfile = prep_data_path + 'Nd_Reff_WuDong_MAGIC.nc'
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
    
    
    ds.attrs["title"] = "cloud droplet number concentration timeseries from Peng Wu and Xiquan Dong's retrieval"
    ds.attrs["description"] = 'median value of each time window'
    ds.attrs["reference"] = 'https://doi.org/10.1175/jcli-d-20-0272.1'
    ds.attrs["inputfile_sample"] = lst[0].split('/')[-1]
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

        