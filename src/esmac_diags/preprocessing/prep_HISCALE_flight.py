"""
prepare aircraft data from HISCALE
options of output data into coarser resolution
"""

import glob
import os
import re
import numpy as np
import time as ttt
from netCDF4 import Dataset
import esmac_diags
from esmac_diags.subroutines.time_format_change import hhmmss2sec
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d,avg_time_1d,avg_time_2d
from esmac_diags.subroutines.read_aircraft import read_ams, read_fims, read_fims_bin, read_iwg1, \
                    read_pcasp, read_cvi_hiscale, read_ccn_hiscale, read_beasd, read_cpc, read_wcm, read_mergedSD
from esmac_diags.subroutines.quality_control import qc_fims_bin, qc_remove_neg, \
                    qc_mask_qcflag, qc_mask_cloudflag


#%% test settings
# iwgpath = '../../../data/HISCALE/obs/aircraft/mei-iwg1/'
# amspath = '../../../data/HISCALE/obs/aircraft/shilling-ams/'
# beasdpath = '../../../data/HISCALE/obs/aircraft/pekour-aafbe/'
# ccnpath = '../../../data/HISCALE/obs/aircraft/mei-ccn/'
# wcmpath = '../../../data/HISCALE/obs/aircraft/matthews-wcm/'
# fimspath = '../../../data/HISCALE/obs/aircraft/wang-fims/'
# pcasppath = '../../../data/HISCALE/obs/aircraft/tomlinson-pcasp/'
# cvipath = '../../../data/HISCALE/obs/aircraft/pekour-cvi/'
# cpcpath = '../../../data/HISCALE/obs/aircraft/mei-cpc/'
# mergeSDpath = '../../../data/HISCALE/obs/aircraft/mergedSD/'
# merged_size_path = 'C:/Users/tang357/Downloads/HISCALE/'
# prep_data_path = 'C:/Users/tang357/Downloads/HISCALE/'

# dt=60

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_AMS(amspath, iwgpath, prep_data_path, dt=60):
    """
    prepare aerosol composition from AMS
    
    Parameters
    ----------
    amspath : str
        input path for aerosol composition data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()

    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in AMS data
        lst2 = glob.glob(amspath+'*'+date[0:8]+'*')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (ams,amslist)=read_ams(lst2[flightidx-1])
            time2=ams[0,:]
            flag=ams[-1,:]
            orgaaf=ams[1,:]
            no3aaf=ams[3,:]
            so4aaf=ams[5,:]
            nh4aaf=ams[7,:]
            chlaaf=ams[9,:]
        elif len(lst2)==0:
            print('does not find any AMS data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')

        # quality controls
        orgaaf=qc_mask_qcflag(orgaaf,flag)
        no3aaf=qc_mask_qcflag(no3aaf,flag)
        so4aaf=qc_mask_qcflag(so4aaf,flag)
        nh4aaf=qc_mask_qcflag(nh4aaf,flag)
        chlaaf=qc_mask_qcflag(chlaaf,flag)
        orgaaf=qc_remove_neg(orgaaf)
        no3aaf=qc_remove_neg(no3aaf)
        so4aaf=qc_remove_neg(so4aaf)
        nh4aaf=qc_remove_neg(nh4aaf)
        chlaaf=qc_remove_neg(chlaaf)
        
        # change values from standardize condition to ambient condition
        T_ams = np.interp(time2,time,T_amb)
        P_ams = np.interp(time2,time,p_amb)
        orgaaf = orgaaf * (296.15/(T_ams+273.15)) * (P_ams/1013.25)
        no3aaf = no3aaf * (296.15/(T_ams+273.15)) * (P_ams/1013.25)
        so4aaf = so4aaf * (296.15/(T_ams+273.15)) * (P_ams/1013.25)
        nh4aaf = nh4aaf * (296.15/(T_ams+273.15)) * (P_ams/1013.25)
        chlaaf = chlaaf * (296.15/(T_ams+273.15)) * (P_ams/1013.25)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        org = median_time_1d(time2, orgaaf, time_new)
        no3 = median_time_1d(time2, no3aaf, time_new)
        so4 = median_time_1d(time2, so4aaf, time_new)
        nh4 = median_time_1d(time2, nh4aaf, time_new)
        chl = median_time_1d(time2, chlaaf, time_new)
        
        #%%
        # import matplotlib.pyplot as plt
        # fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))
        # ax1.plot(time/3600, lwcobs)
        # ax1.plot(time_new/3600, lwc, color='r', marker='.', linewidth=2)
        # # ax1.set_ylim(0, 2e4)
        # ax2.plot(time/3600, twcobs)
        # ax2.plot(time_new/3600, twc, color='r', marker='.', linewidth=2)
        # ax1.set_ylim(-0.01, 0.1)
        # ax1.set_title(date)
        # e
        
        # fig,ax = plt.subplots(1,1,figsize=(6,4))
        # ax.plot(d_merge, np.nanmean(conc_merge,axis=0))
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_xlim(10,2000)
        # plt.grid(True,linestyle=':')

        #%% output data
        outfile = prep_data_path + 'AMS_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        org_o = f.createVariable('ORG', 'f8', ("time",))
        no3_o = f.createVariable('NO3', 'f8', ("time",))
        so4_o = f.createVariable('SO4', 'f8', ("time",))
        nh4_o = f.createVariable('NH4', 'f8', ("time",))
        chl_o = f.createVariable('CHL', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        org_o[:] = org
        no3_o[:] = no3
        so4_o[:] = so4
        nh4_o[:] = nh4
        chl_o[:] = chl
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        org_o.units = 'ug/m3'
        org_o.long_name = 'Organic Matter'
        no3_o.units = 'ug/m3'
        no3_o.long_name = 'Nitrate'
        so4_o.units = 'ug/m3'
        so4_o.long_name = 'Sulfate'
        nh4_o.units = 'ug/m3'
        nh4_o.long_name = 'Ammonium'
        chl_o.units = 'ug/m3'
        chl_o.long_name = 'Chloride'

        # global attributes
        f.title = "Aerosol Composition from AMS"
        f.description = 'median value of each time window, in ambient condition'
        f.input_file = lst2[flightidx-1].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_beasd(beasdpath,iwgpath, prep_data_path, dt=60):
    """
    prepare best-estimate aerosol size distribution data
    
    Parameters
    ----------
    beasdpath : str
        input path for BEASD data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in BEASD data
        lst2 = glob.glob(beasdpath + 'BEASD_G1_' + date[0:8] +'*R2_HISCALE_001s.ict')
        lst2.sort()
        
        if len(lst2)==0:
            print('does not find any BEASD data, skip this date...')
            continue
        
        filename2 = lst2[flightidx-1]
            
        (data2,varlist, dia_merge_h, dia_merge_l, d_merge)=read_beasd(filename2)
        time2 = data2[0,:]
        Cloud_flag = data2[-3,:]
        CVI_flag = data2[-4,:]
        FCDP_flag = data2[-5,:]
        CAS_flag = data2[-6,:]
        PCASP_flag = data2[-7,:]
        FIMS_flag = data2[-8,:]
        conc_merge = data2[1:-8, :]
        
        # quality controls
        conc_merge = qc_remove_neg(conc_merge)
        all_flag = Cloud_flag+CVI_flag+PCASP_flag+FIMS_flag
        conc_merge = qc_mask_qcflag(conc_merge.T, all_flag)
    
        # remove the lowest few bins that FIMS has many false zeros
        conc_merge[:, 0:3] = np.nan
    
        if not all(time2 == time):
            raise ValueError('BEASD time is inconsistent with IWG')
            
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        merge1 = median_time_2d(time, conc_merge, time_new)
        # merge1 = median_time_forflight_2d(time, conc_merge, time_new, height)
    
        #%% output data
        outfile = prep_data_path + 'beasd_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size', len(d_merge)) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        merge_o = f.createVariable('size_distribution_merged', 'f8', ("time", "size"))
        
        # write data
        time_o[:] = time_new
        size_o[:] = d_merge
        sizeh_o[:] = dia_merge_h
        sizel_o[:] = dia_merge_l
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        merge_o[:, :] = merge1
        
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
        merge_o.long_name = 'merged size distribution'
    
        # global attributes
        f.title = "Best Estimate Aerosol Size Distribution"
        f.description = 'median value of each time window, in ambient condition'
        f.description2 = 'The first 3 bins (size<13.7nm) are removed since FIMS report many false zeros'
        f.input_file = filename2.split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CCN(ccnpath, iwgpath, prep_data_path, dt=60):
    """
    prepare CCN number concentration
    
    Parameters
    ----------
    ccnpath : str
        input path for CCN data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in PCASP data
        lst2 = glob.glob(ccnpath+'CCN_G1_'+date[0:8]+'*R2_HiScale001s.*')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (data0,ccnlist)=read_ccn_hiscale(lst2[flightidx-1])
            flag = data0[7,:]
            time2 = data0[0,:]
            ccna = data0[10,:]
            ccnb = data0[11,:]
            SSa = data0[2,:]
            SSb = data0[5,:]
        elif len(lst2)==0:
            print('does not find any CCN data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        if time2.shape != time.shape:
            raise ValueError('CCN time is inconsistent with IWG')
            
        # quality controls
        ccna = qc_mask_qcflag(ccna,flag)
        ccnb = qc_mask_qcflag(ccnb,flag)
        ccna = qc_mask_cloudflag(ccna, cldflag)
        ccnb = qc_mask_cloudflag(ccnb, cldflag)
        SSa=qc_remove_neg(SSa)
        SSb=qc_remove_neg(SSb)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        ccn2 = median_time_1d(time, ccna, time_new)
        ccn5 = median_time_1d(time, ccnb, time_new)
        # ccn2 = median_time_forflight_1d(time, ccna, time_new, height)
        # ccn5 = median_time_forflight_1d(time, ccnb, time_new, height)
        
        #%% output data
        outfile = prep_data_path + 'CCN_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        ccn2_o = f.createVariable('CCN2', 'f8', ("time",))
        ccn5_o = f.createVariable('CCN5', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        ccn2_o[:] = ccn2
        ccn5_o[:] = ccn5
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        ccn2_o.units = '#/cm3'
        ccn2_o.long_name = 'CCN number for SS=0.2%'
        ccn5_o.units = '#/cm3'
        ccn5_o.long_name = 'CCN number for SS=0.5%'
    
        # global attributes
        f.title = "CCN number concentration"
        f.description = 'median value of each time window, in ambient condition'
        f.description2 = 'The actual SS for CCN2 is 0.24%; the actual SS for CCN5 is 0.46%'
        f.input_file = lst2[flightidx-1].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_CPC(cpcpath, iwgpath, prep_data_path, dt=60):
    """
    prepare CPC data
    
    Parameters
    ----------
    cpcpath : str
        input path for CPC data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in CPC data
        lst2 = glob.glob(cpcpath + 'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (cpc,cpclist)=read_cpc(lst2[flightidx-1])
            if date=='20160425a':
                cpc=np.insert(cpc,0,cpc[:,0],axis=1)
                cpc[0,0]=cpc[0,0]-1
            time2 = cpc[0,:]
            cpc10_o = cpc[1,:]
            cpc3_o = cpc[2,:]
        elif len(lst2)==0:
            print('does not find any CPC data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        if time2.shape != time.shape:
            raise ValueError('CPC time is inconsistent with IWG')
            
        # quality checks
        # (cpc3,cpc10) = qc_cpc_air(cpc3,cpc10)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        cpc3 = median_time_1d(time, cpc3_o, time_new)
        cpc10 = median_time_1d(time, cpc10_o, time_new)
        # cpc3 = median_time_forflight_1d(time, cpc3_o, time_new, height)
        # cpc10 = median_time_forflight_1d(time, cpc10_o, time_new, height)
    
        #%% output data
        outfile = prep_data_path + 'CPC_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        c3_o = f.createVariable('cpc3', 'f8', ("time",))
        c10_o = f.createVariable('cpc10', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        c3_o[:] = cpc3
        c10_o[:] = cpc10
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        c3_o.units = '#/cm3'
        c3_o.long_name = 'aerosol number for size > 3nm'
        c10_o.units = '#/cm3'
        c10_o.long_name = 'aerosol number for size > 10nm'
    
        # global attributes
        f.title = "Aerosol number concentration from CPC"
        f.description = 'median value of each time window, in ambient condition'
        f.input_file = lst2[flightidx-1].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mergeSD(mergeSDpath, iwgpath, prep_data_path, dt=60):
    """

    Parameters
    ----------
    mergeSDpath : str
        input path for merged cloud droplet size data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
        
        #%% read in merged SD data
        lst2 = glob.glob(mergeSDpath + 'aaf.g1.hiscale.mergedSD.'+date[0:8]+'*.txt')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (time2, Nd, Ndsize, dmean, dmin, dmax) = read_mergedSD(lst2[flightidx-1])
        elif len(lst2)==0:
            print('does not find any data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        # quality check
        Nd = qc_remove_neg(Nd, remove_zero='True')
        Ndsize = qc_remove_neg(Ndsize)
    
        # unit change from dNd/Dp to Nd
        Dp = dmax-dmin
        for tt in range(Ndsize.shape[1]):   
            Ndsize[:,tt] = Ndsize[:,tt]*Dp
    
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        Nd1 = avg_time_1d(time2, Nd, time_new)
        merge1 = avg_time_2d(time2, Ndsize.T, time_new)
        
            
        #%% 
        # import matplotlib.pyplot as plt
        # fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))
        # h1=ax1.contourf(time2, dmean, np.log(Ndsize), np.arange(-10,13))
        # h2=ax2.contourf(time_new, dmean, np.log(merge1.T), np.arange(-10,13))
        # ax1.set_yscale('log')
        # ax2.set_yscale('log')
        # ax1.set_title(date)
        # cax = plt.axes([0.92, 0.2, 0.02, 0.6])
        # fig.colorbar(h1,cax=cax)
        # ax1.set_xlim(ax2.get_xlim())
        
        # fig,ax=plt.subplots(figsize=(4,3))
        # ax.plot(dmean,np.nanmean(Ndsize,1))
        # ax.plot(dmean,np.nanmean(merge1,0),'r.')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_title(date)
    
        #%% output data
        outfile = prep_data_path + 'mergedSD_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size', len(dmean)) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        nd_o = f.createVariable("Nd", 'f8', ("time", ))
        merge_o = f.createVariable('Nd_bin', 'f8', ("time", "size"))
        
        # write data
        time_o[:] = time_new
        size_o[:] = dmean
        sizeh_o[:] = dmax
        sizel_o[:] = dmin
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        nd_o[:] = Nd1
        merge_o[:, :] = merge1
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        size_o.units = 'um'
        size_o.long_name = 'center of size bin'
        sizeh_o.units = 'um'
        sizeh_o.long_name = 'upper bound of size bin'
        sizel_o.units = 'um'
        sizel_o.long_name = 'lower bound of size bin'
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        nd_o.units = '#/L'
        nd_o.long_name = 'Total cloud droplet number concentration'
        merge_o.units = '#/L'
        merge_o.long_name = 'cloud droplet number concentration by bins'
    
        
        # global attributes
        f.title = "ARM MergedSD product of droplet size distribution from FCDP, 2-DS and HVPS, in ambient condition"
        f.description = 'average of each time window'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
            
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mergesize_HISCALE(fimspath, pcasppath, iwgpath, cvipath, merged_size_path, dt=60):
    """

    Parameters
    ----------
    fimspath : str
        input path for FIMS
    pcasppath : str
        input path for PCASP
    iwgpath : str
        input path for IWG
    cvipath : str
        input path for CVI
    merged_size_path : str
        output path for final merged size
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    if not os.path.exists(merged_size_path):
        os.makedirs(merged_size_path)
            
    # read in fims bin
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'HISCALE_FIMS_bins_R1.dat')
    dlnDp_f = np.empty(len(d_fims))
    for bb in range(len(d_fims)):
        dlnDp_f[bb] = np.log(dmax_f[bb]/dmin_f[bb])
    dlnDp_f = np.mean(dlnDp_f)
    
    # %% find all data
    # lst = glob.glob(iwgpath + 'aaf.iwg1001s.g1.hiscale.20160830*.a2.txt')
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in data
        # IWG
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
    
        # FIMS
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*' + str(flightidx) + '_HISCALE_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            T2 = data0[-1,:]
            p2 = data0[-2,:]*1000.  # kPa to Pa
            # remove some unrealistic data and change data from #/dlnDp to number
            data2 = qc_fims_bin(data0[1:-2, :]) * dlnDp_f
            data2 = qc_remove_neg(data2)
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time > time_fims[-1], time < time_fims[0])
            fims[:, idx] = np.nan
            for dd in range(20):   # remove zeros in small size range
                fims[dd, fims[dd,:]==0] = np.nan
            
            # convert to environment T and p. FIMS measurements are in instrument temperature and pressure
            for tt in range(len(time)):
                fims[:,tt] = fims[:,tt]*((p_amb[tt]/p2[tt])*((T2[tt]+273.15)/(T_amb[tt]+273.15)))
        elif len(filename_f) == 0:
            time_fims = time
            fims = np.nan*np.empty([len(d_fims), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_f)
    
        # PCASP    
        filename_p = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + str(flightidx) + '_hiscale001s.ict.txt')
        # PCASP bins are different in IOP1 and IOP2
        if date[4:6] == '04' or date[4:6] == '05':
            # binlen = 27
            dmax_p = [130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        elif date[4:6] == '08' or date[4:6] == '09':
            # binlen = 30
            dmax_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        # read in data
        if len(filename_p) == 1:
            (data0, pcasplist) = read_pcasp(filename_p[0])
            time_pcasp = data0[0, :]
            d_pcasp = [float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
            # treat flag=1 (caution) as good data
            flag[flag==1]=0
            # remove some questionable data
            pcasp2 = qc_remove_neg(pcasp.T)
            pcasp2 = qc_mask_qcflag(pcasp2, flag)
            pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
            pcasp = pcasp2.T
            if not all(time_pcasp == time):
                raise ValueError('PCASP time is inconsistent with IWG')
        elif len(filename_p) == 0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.nan*np.empty([len(d_pcasp), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_p)
        # !! convert PCASP data to environmental T and p (PCASP standard T/p: [(1013.25/Pamb)*(Tamb/293.15)])
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
            
            
        # CVI
        filename_c = glob.glob(cvipath + 'CVI_G1_' + date[0:8] + '*R4_HISCALE_001s.ict.txt')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1 or len(filename_c) == 2:
            (cvi, cvilist) = read_cvi_hiscale(filename_c[flightidx-1])
            time_cvi = cvi[0, :]
            cvi_inlet = cvi[-1, :]
            # enhance_factor = qc_remove_neg(cvi[2, :])
            # dilution_factor = qc_remove_neg(cvi[3, :])
            cvi_mode = cvi[4, :]
            cvi_qc = cvi[5, :]
            if not all(time_cvi == time):
                raise ValueError('CVI time is inconsistent with IWG')
        elif len(filename_c) == 0:
            time_cvi = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvi_mode = np.nan*np.empty([len(time)])
            cvi_qc = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_c)
        cvi_mode = qc_mask_qcflag(cvi_mode, cvi_qc)
            
        #%% now merge fims and pcasp
        timelen = len(time)
        nbin_merge = 44
        nbin_fims = len(d_fims)
        nbin_pcasp = len(d_pcasp)
        # low and high range of each bin
        dia_merge_l = np.empty(nbin_merge)
        dia_merge_h = np.empty(nbin_merge)
        for n in range(nbin_fims):
            dia_merge_l[n] = dmin_f[n]
            dia_merge_h[n] = dmax_f[n]
        idx = dmax_p.index(500)
        # use upper range (425) of FIMS as low bound and 500 nm of PCASP as high bound
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        for n in range(idx + 1, nbin_pcasp):
            dia_merge_l[nbin_fims + n-idx] = dmin_p[n]
            dia_merge_h[nbin_fims + n-idx] = dmax_p[n]
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration
        conc_merge = np.empty([timelen, nbin_merge])
        for k in range(timelen):
            # use fims data up to d_fims[23] (~0.19 um)
            for n in range(23 + 1):
                if cvi_inlet[k] == 0: 
                    fims[n, k] = np.nan
                conc_merge[k, n] = fims[n, k]
            # overlapping bins
            idx = dmin_p.index(200)   # start merging size. corresponding to 10 in IOP2
            if fims[24, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.95
                    pfac = 0.05
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 24] = (fims[24, k]*ffac + (pcasp[idx, k]*1.0 + pcasp[idx + 1, k]*0.25)*pfac)
            if fims[25, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.8
                    pfac = 0.2
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 25] = (fims[25, k]*ffac + (pcasp[idx + 1, k]*0.75 + pcasp[idx + 2, k]*0.8)*pfac)
            if fims[26, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.65
                    pfac = 0.35
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 26] = (fims[26, k]*ffac + (pcasp[idx + 2, k]*0.2 + pcasp[idx + 3, k]*1.0 + pcasp[idx + 4, k]*0.5)*pfac)
            if fims[27, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.35
                    pfac = 0.65
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 27] = (fims[27, k]*ffac + (pcasp[idx + 4, k]*0.5 + pcasp[idx + 5, k]*0.25)*pfac)
            if fims[28, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.2
                    pfac = 0.8
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 28] = (fims[28, k]*ffac + (pcasp[idx + 5, k]*0.5)*pfac)
            if fims[29, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.05
                    pfac = 0.95
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 29] = (fims[29, k]*ffac + (pcasp[idx + 5, k]*0.25 + pcasp[idx + 6, k]*0.25)*pfac)
            conc_merge[k, 30] = pcasp[idx + 6, k]*0.75
            # using PCASP for upper bins
            nn = 31
            for n in range(idx + 7, nbin_pcasp):
                conc_merge[k, nn] = pcasp[n, k]
                nn = nn + 1
            
        #%% remove the first few bins since FIMS report many false zeros
        conc_merge[:,0:3] = np.nan
        
        #%% re-shape the data into coarser resolution
        conc_merge = qc_remove_neg(conc_merge, remove_zero='False')
        # remove all CVI inlet
        conc_merge[cvi_inlet!=1, :] = np.nan
        # remove all cloud detection
        conc_merge[cldflag!=0, :] = np.nan
        conc_merge[np.isnan(fims[0,:]*pcasp[0,:]), :] = np.nan
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        merge1 = median_time_2d(time, conc_merge, time_new)
        # merge1 = median_time_forflight_2d(time, conc_merge, time_new, height)
        
        
        #%% output data
        outfile = merged_size_path + 'merged_bin_fims_pcasp_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size', nbin_merge) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        merge_o = f.createVariable('size_distribution_merged', 'f8', ("time", "size"))
        
        # write data
        time_o[:] = time_new
        size_o[:] = d_merge
        sizeh_o[:] = dia_merge_h
        sizel_o[:] = dia_merge_l
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        merge_o[:, :] = merge1
        
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
        merge_o.long_name = 'merged size distribution'
    
        
        # global attributes
        f.title = "Merged size distribution from FIMS and PCASP, in ambient condition"
        f.description = 'remove all cloud flags and CVI inlet measurements. median value of each time window'
        f.description2 = 'The first 3 bins (size<13.7nm) are removed since FIMS report many false zeros'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mergesize_HISCALE_withCPC(cpcpath, fimspath, pcasppath, iwgpath, cvipath, merged_size_path, dt=60):
    """

    Parameters
    ----------
    cpcpath : str
        input path for CPC
    fimspath : str
        input path for FIMS
    pcasppath : str
        input path for PCASP
    iwgpath : str
        input path for IWG
    cvipath : str
        input path for CVI
    merged_size_path : str
        output path for final merged size
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
    if not os.path.exists(merged_size_path):
        os.makedirs(merged_size_path)
            
    # read in fims bin
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'HISCALE_FIMS_bins_R1.dat')
    dlnDp_f = np.empty(len(d_fims))
    for bb in range(len(d_fims)):
        dlnDp_f[bb] = np.log(dmax_f[bb]/dmin_f[bb])
    dlnDp_f = np.mean(dlnDp_f)
    
    # %% find all data
    # lst = glob.glob(iwgpath + 'aaf.iwg1001s.g1.hiscale.20160830*.a2.txt')
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in data
        # IWG
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
    
        # FIMS
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*' + str(flightidx) + '_HISCALE_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            T2 = data0[-1,:]
            p2 = data0[-2,:]*1000.  # kPa to Pa
            # remove some unrealistic data and change data from #/dlnDp to number
            data2 = qc_fims_bin(data0[1:-2, :]) * dlnDp_f
            data2 = qc_remove_neg(data2)
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time > time_fims[-1], time < time_fims[0])
            fims[:, idx] = np.nan
            fims[:, fims[20,:]==0] = np.nan
            # convert to environment T and p. FIMS measurements are in instrument temperature and pressure
            for tt in range(len(time)):
                fims[:,tt] = fims[:,tt]*((p_amb[tt]/p2[tt])*((T2[tt]+273.15)/(T_amb[tt]+273.15)))
        elif len(filename_f) == 0:
            time_fims = time
            fims = np.nan*np.empty([len(d_fims), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_f)
    
        # PCASP    
        filename_p = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + str(flightidx) + '_hiscale001s.ict.txt')
        # PCASP bins are different in IOP1 and IOP2
        if date[4:6] == '04' or date[4:6] == '05':
            # binlen = 27
            dmax_p = [130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        elif date[4:6] == '08' or date[4:6] == '09':
            # binlen = 30
            dmax_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        # read in data
        if len(filename_p) == 1:
            (data0, pcasplist) = read_pcasp(filename_p[0])
            time_pcasp = data0[0, :]
            d_pcasp = [float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
            # treat flag=1 (caution) as good data
            flag[flag==1]=0
            # remove some questionable data
            pcasp2 = qc_remove_neg(pcasp.T)
            pcasp2 = qc_mask_qcflag(pcasp2, flag)
            pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
            pcasp = pcasp2.T
            if not all(time_pcasp == time):
                raise ValueError('PCASP time is inconsistent with IWG')
        elif len(filename_p) == 0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.nan*np.empty([len(d_pcasp), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_p)
        # !! convert PCASP data to environmental T and p (PCASP standard T/p: [(1013.25/Pamb)*(Tamb/293.15)])
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
            
            
        # CVI
        filename_c = glob.glob(cvipath + 'CVI_G1_' + date[0:8] + '*R4_HISCALE_001s.ict.txt')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1 or len(filename_c) == 2:
            (cvi, cvilist) = read_cvi_hiscale(filename_c[flightidx-1])
            time_cvi = cvi[0, :]
            cvi_inlet = cvi[-1, :]
            # enhance_factor = qc_remove_neg(cvi[2, :])
            # dilution_factor = qc_remove_neg(cvi[3, :])
            cvi_mode = cvi[4, :]
            cvi_qc = cvi[5, :]
            if not all(time_cvi == time):
                raise ValueError('CVI time is inconsistent with IWG')
        elif len(filename_c) == 0:
            time_cvi = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvi_mode = np.nan*np.empty([len(time)])
            cvi_qc = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_c)
        cvi_mode = qc_mask_qcflag(cvi_mode, cvi_qc)
            
        #%% read in CPC data
        lst2 = glob.glob(cpcpath + 'CPC_G1_'+date[0:8]+'*R2_HiScale001s.ict.txt')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (cpc,cpclist)=read_cpc(lst2[flightidx-1])
            if date=='20160425a':
                cpc=np.insert(cpc,0,cpc[:,0],axis=1)
                cpc[0,0]=cpc[0,0]-1
            time2 = cpc[0,:]
            cpc10 = cpc[1,:]
            cpc3 = cpc[2,:]
        elif len(lst2)==0:
            time2 = time
            cpc10 = np.nan*np.empty([len(time)])
            cpc3 = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find too many files')
    
        if time2.shape != time.shape:
            raise ValueError('CPC time is inconsistent with IWG')
        # quality controls
        cpc10 = qc_remove_neg(cpc10)
        cpc3 = qc_mask_cloudflag(cpc3, cldflag)
            
        #%% now merge fims and pcasp
        timelen = len(time)
        nbin_merge = 44
        nbin_fims = len(d_fims)
        nbin_pcasp = len(d_pcasp)
        # low and high range of each bin
        dia_merge_l = np.empty(nbin_merge)
        dia_merge_h = np.empty(nbin_merge)
        for n in range(nbin_fims):
            dia_merge_l[n] = dmin_f[n]
            dia_merge_h[n] = dmax_f[n]
        idx = dmax_p.index(500)
        # use upper range (425) of FIMS as low bound and 500 nm of PCASP as high bound
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        for n in range(idx + 1, nbin_pcasp):
            dia_merge_l[nbin_fims + n-idx] = dmin_p[n]
            dia_merge_h[nbin_fims + n-idx] = dmax_p[n]
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration
        conc_merge = np.empty([timelen, nbin_merge])
        for k in range(timelen):
            # use fims data up to d_fims[23] (~0.19 um)
            for n in range(23 + 1):
                if cvi_inlet[k] == 0: 
                    fims[n, k] = np.nan
                conc_merge[k, n] = fims[n, k]
            # overlapping bins
            idx = dmin_p.index(200)   # start merging size. corresponding to 10 in IOP2
            if fims[24, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.95
                    pfac = 0.05
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 24] = (fims[24, k]*ffac + (pcasp[idx, k]*1.0 + pcasp[idx + 1, k]*0.25)*pfac)
            if fims[25, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.8
                    pfac = 0.2
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 25] = (fims[25, k]*ffac + (pcasp[idx + 1, k]*0.75 + pcasp[idx + 2, k]*0.8)*pfac)
            if fims[26, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.65
                    pfac = 0.35
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 26] = (fims[26, k]*ffac + (pcasp[idx + 2, k]*0.2 + pcasp[idx + 3, k]*1.0 + pcasp[idx + 4, k]*0.5)*pfac)
            if fims[27, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.35
                    pfac = 0.65
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 27] = (fims[27, k]*ffac + (pcasp[idx + 4, k]*0.5 + pcasp[idx + 5, k]*0.25)*pfac)
            if fims[28, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.2
                    pfac = 0.8
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 28] = (fims[28, k]*ffac + (pcasp[idx + 5, k]*0.5)*pfac)
            if fims[29, k] > 0:
                if cvi_inlet[k] == 1:
                    ffac = 0.05
                    pfac = 0.95
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 29] = (fims[29, k]*ffac + (pcasp[idx + 5, k]*0.25 + pcasp[idx + 6, k]*0.25)*pfac)
            conc_merge[k, 30] = pcasp[idx + 6, k]*0.75
            # using PCASP for upper bins
            nn = 31
            for n in range(idx + 7, nbin_pcasp):
                conc_merge[k, nn] = pcasp[n, k]
                nn = nn + 1
            
        #%% remove the first few bins and use CPC to construct 3-10nm and 10-*nm bins
        dia_merge_l = np.insert(dia_merge_l[4:], 0, [3, 10])
        dia_merge_h = np.insert(dia_merge_h[3:], 0, 10)
        d_merge = (dia_merge_l+dia_merge_h)/2
        nbin_merge = len(d_merge)
        
        conc_merge = qc_remove_neg(conc_merge, remove_zero='False')
        for dd in range(20):
            conc_merge[fims[dd,:]<=0, dd] = np.nan
        conc_merge2 = np.array(conc_merge[:, 4:])
        conc_merge3 = np.insert(conc_merge2, 0, cpc10-np.nansum(conc_merge2,axis=1), axis=1)
        conc_merge4 = np.insert(conc_merge3, 0, cpc3-cpc10, axis=1)
        
        #%% re-shape the data into coarser resolution
        conc_merge = qc_remove_neg(conc_merge4, remove_zero='False')
        # remove all CVI inlet
        conc_merge[cvi_inlet!=1, :] = np.nan
        # remove all cloud detection
        conc_merge[cldflag!=0, :] = np.nan
        conc_merge[np.isnan(fims[0,:]*pcasp[0,:]), :] = np.nan
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        merge1 = median_time_2d(time, conc_merge, time_new)
        # merge1 = median_time_forflight_2d(time, conc_merge, time_new, height)
        
        
        #%% output data
        outfile = merged_size_path + 'merged_bin_cpc_fims_pcasp_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        f.createDimension('size', nbin_merge) 
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        merge_o = f.createVariable('size_distribution_merged', 'f8', ("time", "size"))
        
        # write data
        time_o[:] = time_new
        size_o[:] = d_merge
        sizeh_o[:] = dia_merge_h
        sizel_o[:] = dia_merge_l
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        merge_o[:, :] = merge1
        
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
        merge_o.long_name = 'merged size distribution'
    
        
        # global attributes
        f.title = "Merged size distribution from CPC, FIMS and PCASP, in ambient condition"
        f.description = 'remove all cloud flags and CVI inlet measurements. median value of each time window'
        f.description2 = 'size bin 3-10nm is calculated by CPC3nm - CPC10nm; size bin 10-15.6nm is calculated by '+\
            'CPC10nm - accumulation of CN for size>15.6nm'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_PCASP100(pcasppath, iwgpath, prep_data_path, dt=60):
    """
    prepare aerosol number concentration for size>100nm from PCASP
    
    Parameters
    ----------
    pcasppath : str
        input path for PCASP data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in PCASP data
        lst2 = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*_hiscale001s.ict.txt')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (data0, pcasplist) =read_pcasp(lst2[flightidx-1])
            time2 = data0[0,:]
            d_pcasp = [1000*float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
        elif len(lst2)==0:
            print('does not find any PCASP data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        if time2.shape != time.shape:
            raise ValueError('PCASP time is inconsistent with IWG')
        
        # quality control
        pcasp2 = qc_remove_neg(pcasp.T)
        flag[flag==1] = 0
        pcasp2 = qc_mask_qcflag(pcasp2, flag)
        pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
        pcasp = pcasp2.T
        
        # calculate PCASP for aerosol size>100
        if date[4:6] == '04' or date[4:6] == '05':   # for HISCALE IOP1, PCASP size starts from 120nm
            pcasp100_o = np.sum(pcasp, axis=0)
        elif date[4:6] == '08' or date[4:6] == '09':   # for HISCALE IOP2, PCASP size starts from 90nm
            pcasp100_o = np.sum(pcasp[1:, :], axis=0)
        
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        for tt in range(len(time)):
            pcasp100_o[tt] = pcasp100_o[tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
        
        # pcasp100_o[pcasp100_o==0]=np.nan
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        pcasp100 = median_time_1d(time, pcasp100_o, time_new)
        # pcasp100 = median_time_forflight_1d(time, pcasp100_o, time_new, height)
        
        #%% output data
        outfile = prep_data_path + 'PCASP100_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        p_o = f.createVariable('pcasp100', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        p_o[:] = pcasp100
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        p_o.units = '#/cm3'
        p_o.long_name = 'aerosol number for size > 100nm'
    
        # global attributes
        f.title = "Aerosol number concentration measured by PCASP"
        f.description = 'median value of each time window, in ambient condition'
        f.input_file = lst2[flightidx-1].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_WCM(wcmpath, iwgpath, prep_data_path, dt=60):
    """
    prepare cloud liquid water and total water content data
    
    Parameters
    ----------
    wcmpath : str
        input path for water content measurement data
    iwgpath : str
        input path for IWG
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
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('hiscale.|.a2', filename)
        date = fname[-2]
        print(date)
        if date[-1] == 'a':
            flightidx = 1
        else:
            flightidx = 2
        
        #%% read in IWG data
        (iwg, iwgvars) = read_iwg1(filename)
        timelen = len(iwg)
        # get lat, lon, height, time
        lon = np.empty(timelen)
        lat = np.empty(timelen)
        height = np.empty(timelen)
        time = np.empty(timelen)
        cldflag = np.empty(timelen)
        legnum = np.empty(timelen)
        T_amb = np.empty(timelen)
        p_amb = np.empty(timelen)
        for t in range(timelen):
            lat[t] = float(iwg[t][2])
            lon[t] = float(iwg[t][3])
            height[t] = float(iwg[t][4])
            T_amb[t] = float(iwg[t][20])
            p_amb[t] = float(iwg[t][23])
            cldflag[t] = int(iwg[t][35])
            legnum[t] = int(iwg[t][-1])
            timestr = iwg[t][1].split(' ')
            time[t] = hhmmss2sec(timestr[1])
            
        #%% read in WCM data
        lst2 = glob.glob(wcmpath+'WCM_G1_'+date[0:8]+'*')
        lst2.sort()
        
        if len(lst2)==1 or len(lst2)==2: # some days have two flights
            (wcm,wcmlist)=read_wcm(lst2[flightidx-1])
            time2=wcm[0,:]
            flag=wcm[-1,:]
            twcobs=wcm[1,:]
            lwcobs=wcm[2,:]
        elif len(lst2)==0:
            print('does not find any data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        if not all(time2 == time):
            raise ValueError('WCM time is inconsistent with IWG')
                
        # quality controls
        twcobs=qc_mask_qcflag(twcobs,flag)
        lwcobs=qc_mask_qcflag(lwcobs,flag)
        twcobs=qc_remove_neg(twcobs)
        lwcobs=qc_remove_neg(lwcobs)
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        twc = avg_time_1d(time, twcobs, time_new)
        lwc = avg_time_1d(time, lwcobs, time_new)
        # twc = median_time_forflight_1d(time2, twcobs, time_new, height)
        # lwc = median_time_forflight_1d(time2, lwcobs, time_new, height)
    
        #%% output data
        outfile = prep_data_path + 'WCM_HISCALE_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        twc_o = f.createVariable('TWC', 'f8', ("time",))
        lwc_o = f.createVariable('LWC', 'f8', ("time",))
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        twc_o[:] = twc
        lwc_o[:] = lwc
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        twc_o.units = 'g/m3'
        twc_o.long_name = 'Total water content'
        lwc_o.units = 'g/m3'
        lwc_o.long_name = 'Liquid water content'
    
        # global attributes
        f.title = "cloud water content from WCM"
        f.description = 'median value of each time window'
        f.input_file = lst2[flightidx-1].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()

#%% main for all        
# prep_AMS(amspath, iwgpath, prep_data_path)        
# prep_beasd(beasdpath,iwgpath, prep_data_path)        
# prep_CCN(ccnpath, iwgpath, prep_data_path)
# prep_CPC(cpcpath, iwgpath, prep_data_path)
# prep_mergeSD(mergeSDpath, iwgpath, prep_data_path)
# prep_mergesize_HISCALE(fimspath, pcasppath, iwgpath, cvipath, merged_size_path)
# prep_mergesize_HISCALE_withCPC(cpcpath, fimspath, pcasppath, iwgpath, cvipath, merged_size_path)
# prep_PCASP100(pcasppath, iwgpath, prep_data_path)
# prep_WCM(wcmpath, iwgpath, prep_data_path)