"""
prepare ship data from MARCUS
options of output data into coarser resolution
"""

import glob
import os
import xarray as xr
import pandas as pd
import numpy as np
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, median_time_1d, median_time_2d
from esmac_diags.subroutines.quality_control import  qc_remove_neg, qc_mask_qcflag

# shipmetpath = '../../../data/MARCUS/obs/ship/maraadmetX1.b1/'
# mwrpath = '../../../data/MARCUS/obs/ship/marmwrret1liljclouM1.s2/'
# cpcpath = '../../../data/MARCUS/obs/ship/maraoscpcf1mM1.b1/'
# ccnpath = '../../../data/MARCUS/obs/ship/maraosccn1colavgM1.b1/'
# uhsaspath = '../../../data/MARCUS/obs/ship/maraosuhsasM1.a1/'
# exhaustfreepath = '../../../data/MARCUS/obs/ship/ship_exhaustfree/'
# prep_data_path = 'C:/Users/tang357/Downloads/MARCUS/'

# dt=3600

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
exhaustfreepath : str
    input path of exhaust-free data
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    filename_exhaustfree = exhaustfreepath+'CCN_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    ccn1s = obsdata['CCN1'].load()
    ccn2s = obsdata['CCN2'].load()
    ccn5s = obsdata['CCN5'].load()
    obsdata.close()
    
    ccn1s = qc_remove_neg(ccn1s)
    ccn2s = qc_remove_neg(ccn2s)
    ccn5s = qc_remove_neg(ccn5s)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    ccn1 = median_time_1d(time2, ccn1s, time_new)
    ccn2 = median_time_1d(time2, ccn2s, time_new)
    ccn5 = median_time_1d(time2, ccn5s, time_new)
    
    #%% 
    # import matplotlib.pyplot as plt
    
    # fig = plt.figure(figsize=(8,6))
    # ax1 = fig.add_subplot(3, 1, 1)
    # ax1.plot(time2, ccn1s)
    # ax1.plot(time_new, ccn1, color='r', marker='.',linewidth=2)
    # # ax1.set_ylim(.5, .7)
    # ax2 = fig.add_subplot(3, 1, 2)
    # ax2.plot(time2, ccn5s)
    # ax2.plot(time_new, ccn5, color='r', marker='.',linewidth=2)
    # # ax2.set_ylim(-50, 2000)
    # ax3= fig.add_subplot(3, 1, 3)
    # ax3.plot(time2, ccn2s)
    # ax3.plot(time_new, ccn2, color='r', marker='.',linewidth=2)
    # # ax2.set_ylim(-50, 2000)
    
    # # fig = plt.figure(figsize=(8,3))
    # # ax1 = fig.add_subplot(1, 1, 1)
    # # h1=ax1.contourf(time_new, size.data, uhsas1.T, np.arange(0,200,20))
    # # fig.colorbar(h1)
    # # ax1.set_yscale('log')
    # # import matplotlib.dates as mdates
    # # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    #%% output file
    outfile = prep_data_path + 'CCN_MARCUS_exhaustfree.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lat': (['time'], lat1),
                    'lon': (['time'], lon1),
                    'CCN1': (['time'], ccn1.data),
                    'CCN2': (['time'], ccn2.data),
                    'CCN5': (['time'], ccn5.data),
                    },
                     coords={'time': ('time', time_new)})
    
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lat'].attrs["long_name"] = "latitude"
    ds['lat'].attrs["units"] = "degree_north"
    ds['lon'].attrs["long_name"] = "longitude"
    ds['lon'].attrs["units"] = "degree_east"
    ds['CCN1'].attrs["long_name"] = ccn1s.long_name
    ds['CCN1'].attrs["description"] = ccn1s.description
    ds['CCN1'].attrs["units"] = "1/cm3"
    ds['CCN2'].attrs["long_name"] = ccn2s.long_name
    ds['CCN2'].attrs["description"] = ccn2s.description
    ds['CCN2'].attrs["units"] = "1/cm3"
    ds['CCN5'].attrs["long_name"] = ccn5s.long_name
    ds['CCN5'].attrs["description"] = ccn5s.description
    ds['CCN5'].attrs["units"] = "1/cm3"
    
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    lst2 = glob.glob(ccnpath+'maraosccn1colavgM1.b1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    ccn = obsdata['N_CCN'].load()
    qc_ccn = obsdata['qc_N_CCN'].load()
    SS = obsdata['supersaturation_calculated'].load()
    obsdata.close()
    ccn = qc_mask_qcflag(ccn, qc_ccn)
    
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
    
    #%% output file
    outfile = prep_data_path + 'CCN_MARCUS.nc'
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    lst1 = glob.glob(cpcpath+'maraoscpcf1mM1.b1.*')
    obsdata = xr.open_mfdataset(lst1, combine='by_coords')
    time1 = obsdata['time'].load()
    cpc = obsdata['concentration'].load()
    qc_cpc = obsdata['qc_concentration'].load()
    obsdata.close()
    cpc = qc_remove_neg(cpc.data)
    cpc = qc_mask_qcflag(cpc, qc_cpc)
    
    lst2 = glob.glob(uhsaspath+'maraosuhsasM1.a1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    uhsas = obsdata['concentration'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    obsdata.close()
    uhsas = qc_remove_neg(uhsas.data)
    uhsas100 = np.nansum(uhsas[:, dmin>100], axis=1)
    uhsas100 = qc_remove_neg(uhsas100, remove_zero='True')
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    cpc1 = median_time_1d(time1, cpc, time_new)
    uhsas1 = median_time_1d(time2, uhsas100, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CN_MARCUS.nc'
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
def prep_CN_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    exhaustfreepath : str
        input path of exhaust-free data
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in data
    filename_exhaustfree = exhaustfreepath + 'CPC_UHSAS_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    cpc = obsdata['CPC'].load()
    uhsas = obsdata['UHSAS100'].load()
    obsdata.close()
    
    cpc = qc_remove_neg(cpc.data)
    uhsas = qc_remove_neg(uhsas.data)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    cpc1 = median_time_1d(time2, cpc, time_new)
    uhsas1 = median_time_1d(time2, uhsas, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CN_MARCUS_exhaustfree.nc'
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
    
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
    ds.attrs["description"] = 'median value of '+str(int(dt))+'sec resolution'
    ds.attrs["creation_date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CNsize_exhaustfree(shipmetpath, exhaustfreepath, prep_data_path, dt=3600):
    """
    prepare surface aerosol size distribution
    
    Parameters
    ----------
    shipmetpath : str
        input path for ship location data
    exhaustfreepath : str
        input path of exhaust-free UHSAS size distribution data
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in UHSAS data
    filename_exhaustfree = exhaustfreepath + 'CPC_UHSAS_exhaustfree_1hr.nc'
    obsdata = xr.open_dataset(filename_exhaustfree)
    time2 = obsdata['time'].load()
    size = obsdata['size'].load()
    dmin = obsdata['size_low'].load()
    dmax = obsdata['size_high'].load()
    uhsas = obsdata['UHSAS'].load()
    obsdata.close()
    
    uhsas = qc_remove_neg(uhsas.data)
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    uhsas1 = median_time_2d(time2, uhsas, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CNsize_UHSAS_MARCUS_exhaustfree.nc'
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
    
    ds.attrs["input data_example"] = filename_exhaustfree.split('/')[-1]
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
    uhsaspath : str
        input path for UHSAS size distribution data
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in UHSAS data
    lst2 = glob.glob(uhsaspath+'maraosuhsasM1.a1.*')
    obsdata = xr.open_mfdataset(lst2, combine='by_coords')
    time2 = obsdata['time'].load()
    uhsas = obsdata['concentration'].load()
    dmin = obsdata['lower_size_limit'][0,:].load()
    dmax = obsdata['upper_size_limit'][0,:].load()
    obsdata.close()
    
    uhsas = qc_remove_neg(uhsas.data)
    size = (dmin+dmax)/2
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(time2[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(time2[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    lon1 = median_time_1d(time, lon, time_new)
    lat1 = median_time_1d(time, lat, time_new)
    uhsas1 = median_time_2d(time2, uhsas, time_new)
    
    #%% output file
    outfile = prep_data_path + 'CNsize_UHSAS_MARCUS.nc'
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
    
    lst = glob.glob(shipmetpath+'maraadmetX1.b1.*')
    if len(lst)==0:
        raise ValueError('cannot find any data')
    
    shipdata = xr.open_mfdataset(lst, combine='by_coords')
    time = shipdata['time'].load()
    lon = shipdata['lon'].load()
    qc_lon = shipdata['qc_lon'].load()
    lat = shipdata['lat'].load()
    qc_lat = shipdata['qc_lat'].load()
    shipdata.close()
    
    lat = qc_mask_qcflag(lat, qc_lat)
    lon = qc_mask_qcflag(lon, qc_lon)
    
    #%% read in MWR
    lst2 = glob.glob(mwrpath+'marmwrret1liljclouM1.s2.*')
    lst2.sort()
    
    # first data
    obsdata = xr.open_dataset(lst2[0])
    time2 = obsdata['time'].load()
    lwp = obsdata['be_lwp'].load()
    qc_lwp = obsdata['qc_be_lwp'].load()
    obsdata.close()
    for file in lst2[1:]:
        obsdata = xr.open_dataset(file)
        if file[-18:] == '20171223.000018.nc':
            t1 = obsdata['time'].load()
            lwp1 = obsdata['be_lwp'].load()
            qc_lwp1 = obsdata['qc_be_lwp'].load()
            t2 = xr.concat([t1[:358], t1[4082:], t1[366:4082]], dim="time")
            lwp2 = xr.concat([lwp1[:358], lwp1[4082:], lwp1[366:4082]], dim="time")
            qc_lwp2 = xr.concat([qc_lwp1[:358], qc_lwp1[4082:], qc_lwp1[366:4082]], dim="time")
            time2 = xr.concat([time2, t2], dim="time")
            lwp = xr.concat([lwp, lwp2], dim="time")
            qc_lwp = xr.concat([qc_lwp, qc_lwp2], dim="time")
        elif file[-18:] == '20171227.000000.nc':
            t1 = obsdata['time'].load()
            lwp1 = obsdata['be_lwp'].load()
            qc_lwp1 = obsdata['qc_be_lwp'].load()
            time2 = xr.concat([time2, t1[0:3417]], dim="time")
            lwp = xr.concat([lwp, lwp1[0:3417]], dim="time")
            qc_lwp = xr.concat([qc_lwp, qc_lwp1[0:3417]], dim="time")
        else:
            time2 = xr.concat([time2, obsdata['time']], dim="time")
            lwp = xr.concat([lwp, obsdata['be_lwp']], dim="time")
            qc_lwp = xr.concat([qc_lwp, obsdata['qc_be_lwp']], dim="time")
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
    outfile = prep_data_path + 'LWP_MARCUS.nc'
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


        