"""
prepare aircraft data from ACEENA
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
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, avg_time_2d, \
                    median_time_1d, median_time_2d
from esmac_diags.subroutines.read_aircraft import read_ams, read_fims, read_fims_bin, read_iwg1, \
                    read_pcasp, read_cpc, read_wcm, read_opc, read_mergedSD
from esmac_diags.subroutines.read_ARMdata import read_cvi_aceena, read_ccn
from esmac_diags.subroutines.quality_control import qc_fims_bin, qc_remove_neg, \
                    qc_mask_qcflag, qc_mask_cloudflag

# iwgpath = '../../../data/ACEENA/obs/aircraft/IWG/'
# amspath = '../../../data/ACEENA/obs/aircraft/shilling-hrfams/'
# ccnpath = '../../../data/ACEENA/obs/aircraft/ccn_aaf/'
# wcmpath = '../../../data/ACEENA/obs/aircraft/wcm_ACEENA/'
# fimspath = '../../../data/ACEENA/obs/aircraft/FIMS/'
# pcasppath = '../../../data/ACEENA/obs/aircraft/pcasp_g1/'
# cvipath = '../../../data/ACEENA/obs/aircraft/inletcvi/'
# cpcpath = '../../../data/ACEENA/obs/aircraft/cpc_aaf/'
# opcpath = '../../../data/ACEENA/obs/aircraft/opciso/'
# mergeSDpath = '../../../data/ACEENA/obs/aircraft/mergedSD/'
# merged_size_path = 'C:/Users/tang357/Downloads/ACEENA/'
# prep_data_path = 'C:/Users/tang357/Downloads/ACEENA/'

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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        lst2 = glob.glob(amspath+'AceEnaAMS_G1_'+date[0:8]+'*')
        lst2.sort()
        
        if len(lst2)==1:
            (ams,amslist)=read_ams(lst2[0])
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
        # height2 = np.interp(time2,time,height)
        # org = median_time_forflight_1d(time2, orgaaf, time_new, height2)
        # no3 = median_time_forflight_1d(time2, no3aaf, time_new, height2)
        # so4 = median_time_forflight_1d(time2, so4aaf, time_new, height2)
        # nh4 = median_time_forflight_1d(time2, nh4aaf, time_new, height2)
        # chl = median_time_forflight_1d(time2, chlaaf, time_new, height2)
        
        #%% 
        # import matplotlib.pyplot as plt
        # fig = plt.figure(figsize=(8,8))
        # ax1 = fig.add_subplot(5, 1, 1)
        # ax1.plot(time2/3600, orgaaf)
        # ax1.plot(time_new/3600, org, color='r', marker='.',linewidth=2)
        # # ax1.set_ylim(-50, 5000)
        # ax2 = fig.add_subplot(5, 1, 2)
        # ax2.plot(time2/3600, no3aaf)
        # ax2.plot(time_new/3600, no3, color='r', marker='.',linewidth=2)
        # # ax2.set_ylim(-50, 5000)
        # ax3 = fig.add_subplot(5, 1, 3)
        # ax3.plot(time2/3600, so4aaf)
        # ax3.plot(time_new/3600, so4, color='r', marker='.',linewidth=2)
        # ax4 = fig.add_subplot(5, 1, 4)
        # ax4.plot(time2/3600, nh4aaf)
        # ax4.plot(time_new/3600, nh4, color='r', marker='.',linewidth=2)
        # ax5 = fig.add_subplot(5, 1, 5)
        # ax5.plot(time2/3600, chlaaf)
        # ax5.plot(time_new/3600, chl, color='r', marker='.',linewidth=2)
        
        # # h1=ax1.contourf(time/3600, d_merge, conc_merge.T, np.arange(0,160, 10))
        # # fig.colorbar(h1)
        # # ax2 = fig.add_subplot(2, 1, 2)
        # # h2=ax2.contourf(time_new/3600, d_merge, merge1.T, np.arange(0,160, 10))
        # # ax1.set_yscale('log')
        # # ax2.set_yscale('log')
        # # fig.colorbar(h2)
        
        # ax1.set_title(date)
        
        #%% output data
        outfile = prep_data_path + 'AMS_ACEENA_' + date + '.nc'
        
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
        f.input_file = lst2[0].split('/')[-1]
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
            
        #%% read in CCN data
        filename_ccna=glob.glob(ccnpath+'enaaafccn2colaF1.b1.'+date[0:8]+'*.nc')
        filename_ccnb=glob.glob(ccnpath+'enaaafccn2colbF1.b1.'+date[0:8]+'*.nc')
        
        if len(filename_ccna)==0 and len(filename_ccnb)==0:
            print('does not find any CCN data, skip this date...')
            continue
        
        if len(filename_ccna)==1:
            (timea,timeunita,ccna,qcflag,ccnunit,SSa)=read_ccn(filename_ccna[0])
            ccna=qc_mask_qcflag(ccna,qcflag)
        else:
            raise ValueError('find too many files: '+filename_ccna)
        if len(filename_ccnb)==1:
            (timeb,timeunitb,ccnb,qcflag,ccnunit,SSb)=read_ccn(filename_ccnb[0])
            ccnb=qc_mask_qcflag(ccnb,qcflag)
        else:
            raise ValueError('find too many files: '+filename_ccnb) 
    
        ccna=ccna.data
        ccnb=ccnb.data
        SSa=SSa.data
        SSb=SSb.data
        if any(timea != timeb):
            raise ValueError('time of CCNa and CCNb is inconsistent')
            
        # CCN time is slightly different with IWG, change to consistent
        if date=='20170707a':
            ccna=np.insert(ccna,5249,-9999.)
            ccnb=np.insert(ccnb,5249,-9999.)
            SSa=np.insert(SSa,5249,-9999.)
            SSb=np.insert(SSb,5249,-9999.)
            timea = np.insert(timea,5249,(timea[5248]+timea[5249])/2)
            timeb = np.insert(timeb,5249,(timeb[5248]+timeb[5249])/2)
        elif date=='20180201a':
            ccna=np.insert(ccna,3636,-9999.)
            ccnb=np.insert(ccnb,3636,-9999.)
            SSa=np.insert(SSa,3636,-9999.)
            SSb=np.insert(SSb,3636,-9999.)
            timea = np.insert(timea,3636,(timea[3635]+timea[3636])/2)
            timeb = np.insert(timeb,3636,(timeb[3635]+timeb[3636])/2)
            
        if time[-1]>timea[-1]:
            ccna = np.append(ccna, np.full(int(time[-1]-timea[-1]), -9999))
            ccnb = np.append(ccnb, np.full(int(time[-1]-timea[-1]), -9999))
            SSa = np.append(SSa, np.full(int(time[-1]-timea[-1]), -9999))
            SSb = np.append(SSb, np.full(int(time[-1]-timea[-1]), -9999))
        elif time[-1]<timea[-1]:
            ccna = ccna[0:np.where(timea==time[-1])[0][0]+1]
            ccnb = ccnb[0:np.where(timea==time[-1])[0][0]+1]
            SSa = SSa[0:np.where(timea==time[-1])[0][0]+1]
            SSb = SSb[0:np.where(timea==time[-1])[0][0]+1]
    
        if time[0]>timea[0]:
            ccna = ccna[np.where(timea==time[0])[0][0]:]
            ccnb = ccnb[np.where(timea==time[0])[0][0]:]
            SSa = SSa[np.where(timea==time[0])[0][0]:]
            SSb = SSb[np.where(timea==time[0])[0][0]:]
        elif time[0]<timea[0]:
            ccna = np.insert(ccna, np.full(int(timea[0]-time[0]),0), -9999)
            ccnb = np.insert(ccnb, np.full(int(timea[0]-time[0]),0), -9999)
            SSa = np.insert(SSa, np.full(int(timea[0]-time[0]),0), -9999)
            SSb = np.insert(SSb, np.full(int(timea[0]-time[0]),0), -9999)
            
        # quality check
        ccna=qc_remove_neg(ccna)
        SSa=qc_remove_neg(SSa)
        ccnb=qc_remove_neg(ccnb)
        SSb=qc_remove_neg(SSb)
        ccna = qc_mask_cloudflag(ccna, cldflag)
        ccnb = qc_mask_cloudflag(ccnb, cldflag)
    
        SSa_m = np.nanmean(SSa)
        SSb_m = np.nanmean(SSb)
        
        #%% re-shape the data into coarser resolution
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        SS1 = median_time_1d(time, SSa, time_new)
        SS3 = median_time_1d(time, SSb, time_new)
        ccn1 = median_time_1d(time, ccna, time_new)
        ccn3 = median_time_1d(time, ccnb, time_new)
        # ccn1 = median_time_forflight_1d(time, ccna, time_new, height)
        # ccn3 = median_time_forflight_1d(time, ccnb, time_new, height)
    
        
        #%% output data
        outfile = prep_data_path + 'CCN_ACEENA_' + date + '.nc'
        
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', len(time_new))  
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        ccn1_o = f.createVariable('CCN1', 'f8', ("time",))
        ccn3_o = f.createVariable('CCN3', 'f8', ("time",))
        ss1_o = f.createVariable('SS1', 'f8', ("time",))
        ss3_o = f.createVariable('SS3', 'f8', ("time",))
        
        
        # write data
        time_o[:] = time_new
        lon_o[:] = lon1
        lat_o[:] = lat1  
        height_o[:] = height1
        ccn1_o[:] = ccn1
        ccn3_o[:] = ccn3
        ss1_o[:] = SS1
        ss3_o[:] = SS3
        
        # attributes
        time_o.units = "seconds since " + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + " 00:00:00"
        lon_o.units = 'degree east'
        lon_o.long_name = 'Longitude'
        lat_o.units = 'degree north'
        lat_o.long_name = 'Latitude'
        height_o.units = 'm MSL'
        height_o.long_name = 'height'
        ccn1_o.units = '#/cm3'
        ccn1_o.long_name = 'CCN number for SS=SS1 (~0.1%)'
        ccn3_o.units = '#/cm3'
        ccn3_o.long_name = 'CCN number for SS=SS3 (~0.5%)'
        ss1_o.units = '%'
        ss1_o.long_name = 'measured supersaturation for CCN1'
        ss3_o.units = '%'
        ss3_o.long_name = 'measured supersaturation for CCN3'
    
        # global attributes
        f.title = "CCN number concentration"
        f.description = 'median value of each time window, in ambient condition'
        f.description2 = 'The mean SS for CCN1 is '+format(SSa_m,'.2f')+\
            '%; the mean SS for CCN3 is '+format(SSb_m,'.2f')+'%'
        f.input_file1 = [filename_ccna[0].split('/')[-1], filename_ccnb[0].split('/')[-1]]
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        lst2 = glob.glob(cpcpath + 'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')
        lst2.sort()
        
        if len(lst2)==1:
            (cpc,cpclist)=read_cpc(lst2[0])
            if date=='20180216a':
                cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1) 
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
    
        # quality check
        cpc3_o = qc_mask_cloudflag(cpc3_o, cldflag)
        cpc10_o = qc_mask_cloudflag(cpc10_o, cldflag)
    
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        cpc10 = median_time_1d(time, cpc10_o, time_new)
        cpc3 = median_time_1d(time, cpc3_o, time_new)
        # cpc10 = median_time_forflight_1d(time2, cpc10_o, time_new, height)
        # cpc3 = median_time_forflight_1d(time2, cpc3_o, time_new, height)
    
        
        #%% output data
        outfile = prep_data_path + 'CPC_ACEENA_' + date + '.nc'
        
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
        f.input_file = lst2[0].split('/')[-1]
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        lst2 = glob.glob(mergeSDpath + 'aaf.g1.aceena.mergedSD.'+date[0:8]+'*.txt')
        lst2.sort()
        
        if len(lst2)==1:
            (time2, Nd, Ndsize, dmean, dmin, dmax) = read_mergedSD(lst2[0])
        elif len(lst2)==0:
            print('does not find any data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        # quality check
        Nd = qc_remove_neg(Nd)
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
        import matplotlib.pyplot as plt
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4))
        ax1.contourf(time2, dmean, np.log(Ndsize))
        ax2.contourf(time_new, dmean, np.log(merge1.T))
        ax1.set_yscale('log')
        ax2.set_yscale('log')
        ax1.set_title(date)
        
        fig,ax=plt.subplots(figsize=(4,3))
        ax.plot(dmean,np.nanmean(Ndsize,1))
        ax.plot(dmean,np.nanmean(merge1,0),'r.')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(date)
        # e
    
        #%% output data
        outfile = prep_data_path + 'mergedSD_ACEENA_' + date + '.nc'
        
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
        nd_o.units = '#/L'
        nd_o.long_name = 'Total cloud droplet number concentration'
        merge_o.units = '#/L'
        merge_o.long_name = 'cloud droplet number concentration by bins'
    
        
        # global attributes
        f.title = "ARM MergedSD product of droplet size distribution from FCDP, 2-DS and HVPS, in ambient condition"
        f.description = 'average of each time window'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mergesize_ACEENA(fimspath, pcasppath, opcpath, iwgpath, cvipath, merged_size_path, dt=60):
    """

    Parameters
    ----------
    fimspath : str
        input path for FIMS
    pcasppath : str
        input path for PCASP
    opcpath : str
        input path for OPC
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
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'ACEENA_FIMS_bins_R0.dat')
    dlnDp_f = np.empty(len(d_fims))
    for bb in range(len(d_fims)):
        dlnDp_f[bb] = np.log(dmax_f[bb]/dmin_f[bb])
    dlnDp_f = np.mean(dlnDp_f)
        
    #%% find all data
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            # change data from #/dlnDp to number
            data2 = qc_fims_bin(data0[1:-3, :]) * dlnDp_f
            # TD mode or AMB mode. remove TD mode
            TD_AMB = data0[-1, :]
            data2[:, TD_AMB != 0] = np.nan
            T2 = data0[-2,:]
            p2 = data0[-3,:]*1000.
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time>time_fims[-1], time<time_fims[0])
            fims[:, idx] = np.nan
            # Note: FIMS measurements are in instrument temperature and pressure
            for tt in range(len(time)):
                fims[:,tt] = fims[:,tt]*((p_amb[tt]/p2[tt])*((T2[tt]+273.15)/(T_amb[tt]+273.15)))
        elif len(filename_f) == 0:
            time_fims = time
            fims = np.nan*np.empty([len(d_fims), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_f)
            
        #%% read in PCASP data
        binlen = 30
        dmax_p = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200]
        dmin_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            
        lst2 = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + '_aceena001s.ict')
        lst2.sort()
        
        if len(lst2)==1: 
            (data0, pcasplist) =read_pcasp(lst2[0])
            time_pcasp = data0[0,:]
            d_pcasp = [1000*float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
        elif len(lst2)==0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.full([len(d_pcasp), len(time)], np.nan)
            flag = np.array(time)*np.nan
        else:
            raise ValueError('find too many files')
    
        if time_pcasp.shape != time.shape:
            raise ValueError('PCASP time is inconsistent with IWG')
            
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
            
        # quality controls
        pcasp2 = qc_remove_neg(pcasp.T)
        pcasp2 = qc_mask_qcflag(pcasp2, flag)
        pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
        pcasp = pcasp2.T
    
        # CVI
        filename_c = glob.glob(cvipath + 'enaaafinletcviF1.c1.' + date[0:8] + '*.nc')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1:
            (time_c, lon_c, lat_c, alt_c, timeunit_c, cvimode, cvi_inlet, enhance_factor, dilution_factor) = read_cvi_aceena(filename_c[0])
            if date == '20180216a':
                time_c = np.insert(time_c, 1403, (time_c[1402] + time_c[1403])/2)
                cvi_inlet=np.insert(cvi_inlet, 1403, cvi_inlet[1403])
                cvimode=np.insert(cvimode, 1403, cvimode[1403])
            if not all(time_c == time):
                raise ValueError('CVI time is inconsistent with FIMS')
        elif len(filename_c) == 0:
            time_c = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvimode = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_c)
            
        # if cvi_inlet is unfunctional, assume it is isokinetic and use fims as good data
        cvi_inlet[cvi_inlet == -9] = 1  
        
        
        # read OPC
        filename_o = glob.glob(opcpath + 'OPCISO_G1_' + date[0:8] + '*.ict')
        if len(filename_o) == 1:
            (opc, dmin_o, dmax_o, d_opc, opclist) = read_opc(filename_o[0])
            time_o = opc[0, :]
            opc = opc[1:, :]
            opc = qc_remove_neg(opc)
            d_opc = d_opc*1000
            dmin_o = dmin_o*1000
            dmax_o = dmax_o*1000
        else:
            raise ValueError('can not find OPC data or find multiple files: ' + filename_o)
        if date == '20180216a':
            time_o = np.hstack((time_o[0:1403], 47873., time_o[1403:]))
            opc = np.hstack((opc[:, 0:1403], (opc[:, 1402:1403] + opc[:, 1403:1404])/2, opc[:, 1403:]))
        if any(time_o != time):
            raise ValueError('OPC time is inconsistent with FIMS')
        # fill missing timesteps if the number of missing is not too much
        if sum(np.isnan(opc[0, :])) < 0.1*len(time_o):
            for ii in range(len(d_opc)):    
                opc2 = opc[ii, :]
                opc[ii, np.isnan(opc2)] = np.interp(time[np.isnan(opc2)], time[~np.isnan(opc2)], opc[ii, ~np.isnan(opc2)])
        else:
            print('this date does not fill NaN OPC values')
    
        #%% now merge fims and pcasp
        timelen = len(time)
        nbin_merge = 67
        nbin_fims = len(d_fims)
        nbin_pcasp = len(d_pcasp)
        nbin_opc = len(d_opc)
        # low and high range of each bin
        dia_merge_l = np.full(nbin_merge, np.nan)
        dia_merge_h = np.full(nbin_merge, np.nan)
        # from bins 1-30, use FIMS bin
        for n in range(nbin_fims):
            dia_merge_l[n] = dmin_f[n]
            dia_merge_h[n] = dmax_f[n]
        # for the next bin, use upper range (0.64) of FIMS as low bound and 0.8 of PCASP as high bound
        idx = dmax_p.index(800)
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        # next bin uses 0.8 as low bound and high bound of 2nd bin (0.9) of OPC
        dia_merge_l[31] = 800
        dia_merge_h[31] = 900
        # next few bins are merged two OPC bins
        for n in range(1, 6):
            dia_merge_l[31 + n] = dmin_o[n*2]
            dia_merge_h[31 + n] = dmax_o[n*2 + 1]
        # other bins follows OPC bins
        for n in range(12, nbin_opc):
            dia_merge_l[25 + n] = dmin_o[n]
            dia_merge_h[25 + n] = dmax_o[n]
            
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration. do not treat missing as NaN. treat as -9999
        conc_merge = np.full([timelen, nbin_merge], -9999.)
        fims[np.isnan(fims)] = -9999.
        pcasp[np.isnan(pcasp)] = -9999.
        opc[np.isnan(opc)] = -9999.
        for k in range(timelen):
            # mask all data with cloud flag on
            if cldflag[k] != 0:
                continue
            # use fims data up to d_fims[24]
            for n in range(24 + 1):
                if cvi_inlet[k] == 0:     # in Jerome's code it is 0. looks like it should be 1 (CVI in cloud)
                    fims[n, k] = -9999
                conc_merge[k, n] = fims[n, k]
            # overlapping bins
            idx = dmin_p.index(300)   # start merging size. choose the index of pcasp for merging
            if fims[25, k] >=0:
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
            conc_merge[k, 25] = (fims[25, k]*ffac + pcasp[idx, k]*0.3*pfac)
            if fims[26, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.7
                    pfac = 0.3
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 26] = (fims[26, k]*ffac + (pcasp[idx, k]*0.3 + pcasp[idx + 1, k]*0.2)*pfac)
            if fims[27, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.5
                    pfac = 0.5
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 27] = (fims[27, k]*ffac + (pcasp[idx + 1, k]*0.65)*pfac)
            if fims[28, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.3
                    pfac = 0.7
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 28] = (fims[28, k]*ffac + (pcasp[idx + 1, k]*0.15 + pcasp[idx + 2, k]*0.5)*pfac)
            if fims[29, k] >=0:
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
            conc_merge[k, 29] = (fims[29, k]*ffac + (pcasp[idx + 2, k]*0.4 + pcasp[idx + 3, k]*0.2)*pfac)
            conc_merge[k, 30] = pcasp[idx + 3, k]*0.8
            if not all(pcasp[idx:idx + 4, k] >=0):
                conc_merge[k, 25:30] = fims[25:30, k]
                conc_merge[k, 30] = (conc_merge[k, 29] + opc[1, k]*1.4)/2.0
            # next merge OPC and PCASP, remove PCASP if the values is 10x larger than OPC
            pcasp2 = pcasp[18, k]*0.5
            opc2 = opc[1, k]*1.4     # the first bin of OPC contains all small-size particles. not using opc[0, k]
            if opc2<0:             
                conc_merge[k, 31] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 31] = opc2
            else:
                conc_merge[k, 31] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[18, k]*0.5 + pcasp[19, k]*0.2
            opc2 = opc[2, k] + opc[3, k]
            if opc2<0:             
                conc_merge[k, 32] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 32] = opc2
            else:
                conc_merge[k, 32] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[19, k]*0.8
            opc2 = opc[4, k] + opc[5, k]
            if opc2<0:             
                conc_merge[k, 33] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 33] = opc2
            else:
                conc_merge[k, 33] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.9
            opc2 = opc[6, k] + opc[7, k]
            if opc2<0:             
                conc_merge[k, 34] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 34] = opc2
            else:
                conc_merge[k, 34] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.1 + pcasp[21, k]
            opc2 = opc[8, k] + opc[9, k]
            if opc2<0:             
                conc_merge[k, 35] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 35] = opc2
            else:
                conc_merge[k, 35] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[22, k] + pcasp[23, k]*0.2
            opc2 = opc[10, k] + opc[11, k]
            if opc2<0:             
                conc_merge[k, 36] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 36] = opc2
            else:
                conc_merge[k, 36] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.7
            opc2 = opc[12, k]
            if opc2<0:             
                conc_merge[k, 37] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 37] = opc2
            else:
                conc_merge[k, 37] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.1 + pcasp[24, k]*0.7
            opc2 = opc[13, k]
            if opc2<0:             
                conc_merge[k, 38] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 38] = opc2
            else:
                conc_merge[k, 38] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[24, k]*0.3 + pcasp[25, k]*0.4
            opc2 = opc[14, k]
            if opc2<0:             
                conc_merge[k, 39] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 39] = opc2
            else:
                conc_merge[k, 39] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[25, k]*0.6 + pcasp[26, k]*0.3
            opc2 = opc[15, k]
            if opc2<0:             
                conc_merge[k, 40] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 40] = opc2
            else:
                conc_merge[k, 40] = opc2*0.6 + pcasp2*0.4  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[26, k]*0.7 + pcasp[27, k]*0.2
            opc2 = opc[16, k]
            if opc2<0:             
                conc_merge[k, 41] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 41] = opc2
            else:
                conc_merge[k, 41] = opc2*0.7 + pcasp2*0.3  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[27, k]*0.8 + pcasp[28, k]*0.2
            opc2 = opc[17, k]
            if opc2<0:             
                conc_merge[k, 42] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 42] = opc2
            else:
                conc_merge[k, 42] = opc2*0.8 + pcasp2*0.2  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[28, k]*0.8 + pcasp[29, k]*0.3
            opc2 = opc[18, k]
            if opc2<0:             
                conc_merge[k, 43] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 43] = opc2
            else:
                conc_merge[k, 43] = opc2*0.9 + pcasp2*0.1  # gradually reduce the weight of PCASP
            # using OPC for other bins
            for n in range(44, nbin_merge):
                conc_merge[k, n] = opc[n-25, k]
                
        #%% re-shape the data into coarser resolution
        conc_merge = qc_remove_neg(conc_merge, remove_zero='False')
        # remove all CVI inlet
        conc_merge[cvi_inlet!=1, :] = np.nan
        # remove all cloud detection
        conc_merge[cldflag!=0, :] = np.nan
        conc_merge[np.logical_or(fims[0,:]<0, pcasp[0,:]<0), :] = np.nan
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        merge1 = median_time_2d(time, conc_merge, time_new)
        # merge1 = median_time_forflight_2d(time, conc_merge, time_new, height)
        
        #%% output data
        outfile = merged_size_path + 'merged_bin_fims_pcasp_opc_ACEENA' + date + '.nc'
        
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
        f.title = "Merged size distribution from FIMS, PCASP and OPC, in ambient condition"
        f.description = 'remove all cloud flags and CVI inlet measurements. median value of each time window'
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_mergesize_withCPC_ACEENA(cpcpath, fimspath, pcasppath, opcpath, iwgpath, cvipath, merged_size_path, dt=60):
    """
    
    Parameters
    ----------
    cpcpath : str
        input path for CPC
    fimspath : str
        input path for FIMS
    pcasppath : str
        input path for PCASP
    opcpath : str
        input path for OPC
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
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'ACEENA_FIMS_bins_R0.dat')
    dlnDp_f = np.empty(len(d_fims))
    for bb in range(len(d_fims)):
        dlnDp_f[bb] = np.log(dmax_f[bb]/dmin_f[bb])
    dlnDp_f = np.mean(dlnDp_f)
        
    #%% find all data
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    for filename in lst[:]:
        
        # get date
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            # change data from #/dlnDp to number
            data2 = qc_fims_bin(data0[1:-3, :]) * dlnDp_f
            # TD mode or AMB mode. remove TD mode
            TD_AMB = data0[-1, :]
            data2[:, TD_AMB != 0] = np.nan
            T2 = data0[-2,:]
            p2 = data0[-3,:]*1000.
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time>time_fims[-1], time<time_fims[0])
            fims[:, idx] = np.nan
            # Note: FIMS measurements are in instrument temperature and pressure
            for tt in range(len(time)):
                fims[:,tt] = fims[:,tt]*((p_amb[tt]/p2[tt])*((T2[tt]+273.15)/(T_amb[tt]+273.15)))
        elif len(filename_f) == 0:
            time_fims = time
            fims = np.nan*np.empty([len(d_fims), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_f)
        
        #%% read in PCASP data
        binlen = 30
        dmax_p = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200]
        dmin_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            
        lst2 = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + '_aceena001s.ict')
        lst2.sort()
        
        if len(lst2)==1: 
            (data0, pcasplist) =read_pcasp(lst2[0])
            time_pcasp = data0[0,:]
            d_pcasp = [1000*float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
        elif len(lst2)==0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.full([len(d_pcasp), len(time)], np.nan)
            flag = np.array(time)*np.nan
        else:
            raise ValueError('find too many files')
    
        if time_pcasp.shape != time.shape:
            raise ValueError('PCASP time is inconsistent with IWG')
            
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
            
        # quality controls
        pcasp2 = qc_remove_neg(pcasp.T)
        pcasp2 = qc_mask_qcflag(pcasp2, flag)
        pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
        pcasp = pcasp2.T
    
        # CVI
        filename_c = glob.glob(cvipath + 'enaaafinletcviF1.c1.' + date[0:8] + '*.nc')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1:
            (time_c, lon_c, lat_c, alt_c, timeunit_c, cvimode, cvi_inlet, enhance_factor, dilution_factor) = read_cvi_aceena(filename_c[0])
            if date == '20180216a':
                time_c = np.insert(time_c, 1403, (time_c[1402] + time_c[1403])/2)
                cvi_inlet=np.insert(cvi_inlet, 1403, cvi_inlet[1403])
                cvimode=np.insert(cvimode, 1403, cvimode[1403])
            if not all(time_c == time):
                raise ValueError('CVI time is inconsistent with FIMS')
        elif len(filename_c) == 0:
            time_c = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvimode = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_c)
            
        # if cvi_inlet is unfunctional, assume it is isokinetic and use fims as good data
        cvi_inlet[cvi_inlet == -9] = 1  
        
        
        #%% read OPC
        filename_o = glob.glob(opcpath + 'OPCISO_G1_' + date[0:8] + '*.ict')
        if len(filename_o) == 1:
            (opc, dmin_o, dmax_o, d_opc, opclist) = read_opc(filename_o[0])
            time_o = opc[0, :]
            opc = opc[1:, :]
            opc = qc_remove_neg(opc)
            d_opc = d_opc*1000
            dmin_o = dmin_o*1000
            dmax_o = dmax_o*1000
        else:
            raise ValueError('can not find OPC data or find multiple files: ' + filename_o)
        if date == '20180216a':
            time_o = np.hstack((time_o[0:1403], 47873., time_o[1403:]))
            opc = np.hstack((opc[:, 0:1403], (opc[:, 1402:1403] + opc[:, 1403:1404])/2, opc[:, 1403:]))
        if any(time_o != time):
            raise ValueError('OPC time is inconsistent with FIMS')
        # fill missing timesteps if the number of missing is not too much
        if sum(np.isnan(opc[0, :])) < 0.1*len(time_o):
            for ii in range(len(d_opc)):    
                opc2 = opc[ii, :]
                opc[ii, np.isnan(opc2)] = np.interp(time[np.isnan(opc2)], time[~np.isnan(opc2)], opc[ii, ~np.isnan(opc2)])
        else:
            print('this date does not fill NaN OPC values')
    
        #%% read in CPC data
        lst2 = glob.glob(cpcpath + 'CPC_G1_'+date[0:8]+'*R2_ACEENA001s.ict')
        lst2.sort()
        
        if len(lst2)==1:
            (cpc,cpclist)=read_cpc(lst2[0])
            if date=='20180216a':
                cpc=np.insert(cpc,1404,(cpc[:,1403]+cpc[:,1404])/2,axis=1) 
            time2 = cpc[0,:]
            cpc10 = cpc[1,:]
            cpc3 = cpc[2,:]
        elif len(lst2)==0:
            time2 = time
            cpc10 = np.nan*np.empty([len(time)])
            cpc3 = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find too many files')
        # quality controls
        cpc10 = qc_remove_neg(cpc10)
        cpc3 = qc_mask_cloudflag(cpc3, cldflag)
                
        #%% now merge fims and pcasp
        timelen = len(time)
        nbin_merge = 67
        nbin_fims = len(d_fims)
        nbin_pcasp = len(d_pcasp)
        nbin_opc = len(d_opc)
        # low and high range of each bin
        dia_merge_l = np.full(nbin_merge, np.nan)
        dia_merge_h = np.full(nbin_merge, np.nan)
        # from bins 1-30, use FIMS bin
        for n in range(nbin_fims):
            dia_merge_l[n] = dmin_f[n]
            dia_merge_h[n] = dmax_f[n]
        # for the next bin, use upper range (0.64) of FIMS as low bound and 0.8 of PCASP as high bound
        idx = dmax_p.index(800)
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        # next bin uses 0.8 as low bound and high bound of 2nd bin (0.9) of OPC
        dia_merge_l[31] = 800
        dia_merge_h[31] = 900
        # next few bins are merged two OPC bins
        for n in range(1, 6):
            dia_merge_l[31 + n] = dmin_o[n*2]
            dia_merge_h[31 + n] = dmax_o[n*2 + 1]
        # other bins follows OPC bins
        for n in range(12, nbin_opc):
            dia_merge_l[25 + n] = dmin_o[n]
            dia_merge_h[25 + n] = dmax_o[n]
            
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration. do not treat missing as NaN. treat as -9999
        conc_merge = np.full([timelen, nbin_merge], -9999.)
        fims[np.isnan(fims)] = -9999.
        pcasp[np.isnan(pcasp)] = -9999.
        opc[np.isnan(opc)] = -9999.
        for k in range(timelen):
            # mask all data with cloud flag on
            if cldflag[k] != 0:
                continue
            # use fims data up to d_fims[24]
            for n in range(24 + 1):
                if cvi_inlet[k] == 0:     # in Jerome's code it is 0. looks like it should be 1 (CVI in cloud)
                    fims[n, k] = -9999
                conc_merge[k, n] = fims[n, k]
            # overlapping bins
            idx = dmin_p.index(300)   # start merging size. choose the index of pcasp for merging
            if fims[25, k] >=0:
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
            conc_merge[k, 25] = (fims[25, k]*ffac + pcasp[idx, k]*0.3*pfac)
            if fims[26, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.7
                    pfac = 0.3
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 26] = (fims[26, k]*ffac + (pcasp[idx, k]*0.3 + pcasp[idx + 1, k]*0.2)*pfac)
            if fims[27, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.5
                    pfac = 0.5
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 27] = (fims[27, k]*ffac + (pcasp[idx + 1, k]*0.65)*pfac)
            if fims[28, k] >=0:
                if cvi_inlet[k] == 1:
                    ffac = 0.3
                    pfac = 0.7
                elif cvi_inlet[k] == 0:
                    ffac = 0.0
                    pfac = 1.0
                else:
                    raise ValueError('cvi_inlet value is neither 0 nor 1')
            else:
                ffac = 0.0
                pfac = 1.0
            conc_merge[k, 28] = (fims[28, k]*ffac + (pcasp[idx + 1, k]*0.15 + pcasp[idx + 2, k]*0.5)*pfac)
            if fims[29, k] >=0:
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
            conc_merge[k, 29] = (fims[29, k]*ffac + (pcasp[idx + 2, k]*0.4 + pcasp[idx + 3, k]*0.2)*pfac)
            conc_merge[k, 30] = pcasp[idx + 3, k]*0.8
            if not all(pcasp[idx:idx + 4, k] >=0):
                conc_merge[k, 25:30] = fims[25:30, k]
                conc_merge[k, 30] = (conc_merge[k, 29] + opc[1, k]*1.4)/2.0
            # next merge OPC and PCASP, remove PCASP if the values is 10x larger than OPC
            pcasp2 = pcasp[18, k]*0.5
            opc2 = opc[1, k]*1.4     # the first bin of OPC contains all small-size particles. not using opc[0, k]
            if opc2<0:             
                conc_merge[k, 31] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 31] = opc2
            else:
                conc_merge[k, 31] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[18, k]*0.5 + pcasp[19, k]*0.2
            opc2 = opc[2, k] + opc[3, k]
            if opc2<0:             
                conc_merge[k, 32] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 32] = opc2
            else:
                conc_merge[k, 32] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[19, k]*0.8
            opc2 = opc[4, k] + opc[5, k]
            if opc2<0:             
                conc_merge[k, 33] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 33] = opc2
            else:
                conc_merge[k, 33] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.9
            opc2 = opc[6, k] + opc[7, k]
            if opc2<0:             
                conc_merge[k, 34] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 34] = opc2
            else:
                conc_merge[k, 34] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.1 + pcasp[21, k]
            opc2 = opc[8, k] + opc[9, k]
            if opc2<0:             
                conc_merge[k, 35] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 35] = opc2
            else:
                conc_merge[k, 35] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[22, k] + pcasp[23, k]*0.2
            opc2 = opc[10, k] + opc[11, k]
            if opc2<0:             
                conc_merge[k, 36] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 36] = opc2
            else:
                conc_merge[k, 36] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.7
            opc2 = opc[12, k]
            if opc2<0:             
                conc_merge[k, 37] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 37] = opc2
            else:
                conc_merge[k, 37] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.1 + pcasp[24, k]*0.7
            opc2 = opc[13, k]
            if opc2<0:             
                conc_merge[k, 38] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 38] = opc2
            else:
                conc_merge[k, 38] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[24, k]*0.3 + pcasp[25, k]*0.4
            opc2 = opc[14, k]
            if opc2<0:             
                conc_merge[k, 39] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 39] = opc2
            else:
                conc_merge[k, 39] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[25, k]*0.6 + pcasp[26, k]*0.3
            opc2 = opc[15, k]
            if opc2<0:             
                conc_merge[k, 40] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 40] = opc2
            else:
                conc_merge[k, 40] = opc2*0.6 + pcasp2*0.4  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[26, k]*0.7 + pcasp[27, k]*0.2
            opc2 = opc[16, k]
            if opc2<0:             
                conc_merge[k, 41] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 41] = opc2
            else:
                conc_merge[k, 41] = opc2*0.7 + pcasp2*0.3  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[27, k]*0.8 + pcasp[28, k]*0.2
            opc2 = opc[17, k]
            if opc2<0:             
                conc_merge[k, 42] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 42] = opc2
            else:
                conc_merge[k, 42] = opc2*0.8 + pcasp2*0.2  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[28, k]*0.8 + pcasp[29, k]*0.3
            opc2 = opc[18, k]
            if opc2<0:             
                conc_merge[k, 43] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 43] = opc2
            else:
                conc_merge[k, 43] = opc2*0.9 + pcasp2*0.1  # gradually reduce the weight of PCASP
            # using OPC for other bins
            for n in range(44, nbin_merge):
                conc_merge[k, n] = opc[n-25, k]
                
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
        conc_merge[np.logical_or(fims[0,:]<0, pcasp[0,:]<0), :] = np.nan
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        merge1 = median_time_2d(time, conc_merge, time_new)
        # merge1 = median_time_forflight_2d(time, conc_merge, time_new, height)
    
        #%% output data
        outfile = merged_size_path + 'merged_bin_cpc_fims_pcasp_opc_ACEENA_' + date + '.nc'
        
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
        f.title = "Merged size distribution from CPC, FIMS, PCASP and OPC, in ambient condition"
        f.description = 'remove all cloud flags and CVI inlet measurements. median value of each time window'
        f.description2 = 'size bin 3-10nm is calculated by CPC3nm - CPC10nm; size bin 10-16nm is calculated by '+\
            'CPC10nm - accumulation of CN for size>16nm'
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        lst2 = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + '_aceena001s.ict')
        lst2.sort()
        
        if len(lst2)==1: 
            (data0, pcasplist) =read_pcasp(lst2[0])
            time2 = data0[0,:]
            flag = data0[-2, :]
            pcasp100_o = data0[-5, :]  # PCASP for ACEENA starts from 100nm
        elif len(lst2)==0:
            print('does not find any PCASP data, skip this date...')
            continue
        else:
            raise ValueError('find too many files')
    
        if time2.shape != time.shape:
            raise ValueError('PCASP time is inconsistent with IWG')
            
        # quality controls
        pcasp100_o = qc_mask_qcflag(pcasp100_o, flag)
        pcasp100_o = qc_mask_cloudflag(pcasp100_o, cldflag)
        pcasp100_o = qc_remove_neg(pcasp100_o)
        
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        for tt in range(len(time)):
            pcasp100_o[tt] = pcasp100_o[tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
        
        #%% re-shape the data into coarser resolution
        
        time_new = np.arange(np.round(time[0]/dt), np.round(time[-1]/dt)) *dt
        
        lon1 = median_time_1d(time, lon, time_new)
        lat1 = median_time_1d(time, lat, time_new)
        height1 = median_time_1d(time, height, time_new)
        
        pcasp100 = median_time_1d(time, pcasp100_o, time_new)
        # pcasp100 = median_time_forflight_1d(time2, pcasp100_o, time_new, height)
    
        #%% output data
        outfile = prep_data_path + 'PCASP100_ACEENA_' + date + '.nc'
        
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
        f.input_file = lst2[0].split('/')[-1]
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        if date == '20180216a':
            iwg.insert(1403, list(iwg[1403]))
            tstr = iwg[1403][1]
            tstr = tstr[0:-1] + str(int(tstr[-1])-1)
            iwg[1403][1] = tstr
            del iwg[-1]
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
        
        if len(lst2)==1: 
            (wcm,wcmlist)=read_wcm(lst2[0])
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
        outfile = prep_data_path + 'WCM_ACEENA_' + date + '.nc'
        
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
        f.input_file = lst2[0].split('/')[-1]
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()