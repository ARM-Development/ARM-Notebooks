"""
# merge size distribution from FIMS and PCASP for Hi-Scale
# revised from size_bin_merge.pro by Jerome Fast
# Shuaiqi Tang
# 2020.10.1
"""

import os
import glob
import re
import numpy as np
from ..subroutines.read_aircraft import read_fims, read_fims_bin, read_iwg1, read_pcasp, read_cvi_hiscale as read_cvi
from ..subroutines.time_format_change import hhmmss2sec
from ..subroutines.quality_control import qc_fims_bin, qc_mask_qcflag, qc_mask_cloudflag
from netCDF4 import Dataset

def run_prep(settings):
    #%% variables from settings
    iwgpath = settings['iwgpath']
    fimspath = settings['fimspath']
    pcasppath = settings['pcasppath']
    cvipath = settings['cvipath']
    merged_size_path = settings['merged_size_path']

    #%% other settings
    
    if not os.path.exists(merged_size_path):
        os.makedirs(merged_size_path)
        
    
    # %% find all data
    # lst = glob.glob(iwgpath + 'aaf.iwg1001s.g1.hiscale.20160830*.a2.txt')
    lst = glob.glob(iwgpath + '*.a2.txt')
    lst.sort()
    
    # read in fims bin
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'HISCALE_FIMS_bins_R1.dat')
    # change unit to um
    d_fims = [x/1000 for x in d_fims]
    dmin_f = [x/1000 for x in dmin_f]
    dmax_f = [x/1000 for x in dmax_f]
    dlnDp_f = np.empty(len(d_fims))
    for bb in range(len(d_fims)):
        dlnDp_f[bb] = np.log(dmax_f[bb]/dmin_f[bb])
    dlnDp_f = np.mean(dlnDp_f)
    
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
        datestr = timestr[0]
    
        # FIMS
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*' + str(flightidx) + '_HISCALE_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            T2 = data0[-1,:]
            p2 = data0[-2,:]*1000.
            # remove some unrealistic data and change data from #/dlnDp to number
            data2 = qc_fims_bin(data0[1:-2, :]) * dlnDp_f
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time > time_fims[-1], time < time_fims[0])
            fims[:, idx] = np.nan
            #!!! FIMS measurements are in instrument temperature and pressure
            for tt in range(len(time)):
                fims[:,tt] = fims[:,tt]*((p_amb[tt]/p2[tt])*((T2[tt]+273.15)/(T_amb[tt]+273.15)))
        elif len(filename_f) == 0:
            time_fims = time
            fims = np.nan*np.empty([len(d_fims), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_f)
        fims_total = np.nansum(fims, 0)
        fims_total[fims_total <= 0] = np.nan
    
        # PCASP    
        filename_p = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + str(flightidx) + '_hiscale001s.ict.txt')
        if date[4:6] == '04' or date[4:6] == '05':
            binlen = 27
            dmax_p = [130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, 400, 500, \
                    600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        elif date[4:6] == '08' or date[4:6] == '09':
            binlen = 30
            dmax_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
            dmin_p = [90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                    400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800]
        dmin_p = [x/1000 for x in dmin_p]
        dmax_p = [x/1000 for x in dmax_p]
        # read in data
        if len(filename_p) == 1:
            (data0, pcasplist) = read_pcasp(filename_p[0])
            time_pcasp = data0[0, :]
            d_pcasp = [float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
            pcasp_total = data0[-5, :]
            # remove some questionable data
            pcasp2 = pcasp.T
            pcasp2 = qc_mask_qcflag(pcasp2, flag)
            pcasp2 = qc_mask_cloudflag(pcasp2, cldflag)
            pcasp = pcasp2.T
            if not all(time_pcasp == time):
                raise ValueError('PCASP time is inconsistent with FIMS')
        elif len(filename_p) == 0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.nan*np.empty([len(d_pcasp), len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_p)
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
            
            
        # CVI
        filename_c = glob.glob(cvipath + 'CVI_G1_' + date[0:8] + '*R4_HISCALE_001s.ict.txt')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1 or len(filename_c) == 2:
            (cvi, cvilist) = read_cvi(filename_c[flightidx-1])
            time_cvi = cvi[0, :]
            cvi_inlet = cvi[-1, :]
            enhance_factor = cvi[2, :]
            enhance_factor[enhance_factor < -9000] = np.nan
            dilution_factor = cvi[3, :]
            dilution_factor[dilution_factor < -9000] = np.nan
            cvi_mode = cvi[4, :]
            cvi_qc = cvi[5, :]
            if not all(time_cvi == time):
                raise ValueError('CVI time is inconsistent with FIMS')
        elif len(filename_c) == 0:
            time_cvi = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvi_mode = np.nan*np.empty([len(time)])
            dilution_factor = np.nan*np.empty([len(time)])
            enhance_factor = np.nan*np.empty([len(time)])
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
        idx = dmax_p.index(0.5)
        # use upper range (0.425) of FIMS as low bound and 0.5 of PCASP as high bound
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        for n in range(idx + 1, nbin_pcasp):
            dia_merge_l[nbin_fims + n-idx] = dmin_p[n]
            dia_merge_h[nbin_fims + n-idx] = dmax_p[n]
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration
        conc_merge = np.empty([timelen, nbin_merge])
        fims[np.isnan(fims)] = -9999.   # do not treat missing as NaN. treat -9999
        for k in range(timelen):
            # use fims data up to d_fims[23] (~0.19 um)
            for n in range(23 + 1):
                if cvi_inlet[k] == 0:     # in Jerome's code it is 0. looks like it should be 1 (CVI in cloud)
                    fims[n, k] = -9999
                conc_merge[k, n] = fims[n, k]
            # overlapping bins
            idx = dmin_p.index(0.2)   # start merging size. corresponding to 10 in IOP2
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
            
        #%% output data
        if not os.path.exists(merged_size_path):
            os.mkdir(merged_size_path)
        outfile = merged_size_path + 'merged_bin_fims_pcasp_HISCALE_' + date + '.nc'
        # define filename
        f = Dataset(outfile, 'w', format = 'NETCDF4')
        
        # define dimensions
        t = f.createDimension('time', None)  # unlimited
        s = f.createDimension('size', nbin_merge)  # unlimited
        
        # create variable list
        time_o = f.createVariable("time", "f8", ("time", ))
        size_o = f.createVariable("size", "f8", ("size", ))
        sizeh_o = f.createVariable("size_high", "f8", ("size", ))
        sizel_o = f.createVariable("size_low", "f8", ("size", ))
        lon_o = f.createVariable("lon", 'f8', ("time", ))
        lat_o = f.createVariable("lat", 'f8', ("time", ))
        height_o = f.createVariable("height", 'f8', ("time", ))
        cflag_o = f.createVariable('cld_flag', 'i8', ("time", ))
        legnum_o = f.createVariable('leg_number', 'i8', ("time", ))
        cvi_o = f.createVariable('CVI_inlet', 'i8', ("time", ))
        cvim_o = f.createVariable('CVI_mode', 'i8', ("time", ))
        df_o = f.createVariable('CVI_Dilution_Factor', 'f8', ("time", ))
        ef_o = f.createVariable('CVI_Enhancement_Factor', 'f8', ("time", ))
        merge_o = f.createVariable('size_distribution_merged', 'f8', ("time", "size"))
        fims_o = f.createVariable('totalnum_fims', 'f8', ("time", ))
        pcasp_o = f.createVariable('totalnum_pcasp', 'f8', ("time", ))
        
        # write data
        time_o[:] = time
        size_o[:] = d_merge
        sizeh_o[:] = dia_merge_h
        sizel_o[:] = dia_merge_l
        lon_o[:] = lon
        lat_o[:] = lat    
        height_o[:] = height
        cflag_o[:] = cldflag
        legnum_o[:] = legnum
        cvi_o[:] = cvi_inlet
        cvim_o[:] = np.array(cvi_mode)
        dilution_factor[np.isnan(dilution_factor)] = -9999.
        df_o[:] = dilution_factor
        enhance_factor[np.isnan(enhance_factor)] = -9999.
        ef_o[:] = enhance_factor
        conc_merge[np.isnan(conc_merge)] = -9999.
        conc_merge[conc_merge < 0] = -9999.
        merge_o[:, :] = conc_merge
        fims_total[np.isnan(fims_total)] = -9999.
        fims_total[fims_total < 0] = -9999.
        fims_o[:] = fims_total
        pcasp_total[np.isnan(pcasp_total)] = -9999.
        pcasp_total[pcasp_total < 0] = -9999.
        pcasp_o[:] = pcasp_total
        
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
        cflag_o.units = 'N/A'
        cflag_o.long_name = 'cloud flag'
        cflag_o.description = '1-cloud; 0-no cloud'
        legnum_o.units = 'N/A'
        legnum_o.long_name = 'leg number'
        cvi_o.units = 'N/A'
        cvi_o.long_name = 'CVI inlet status'
        cvi_o.description = '0-CVI inlet on; 1-Isokinetic inlet on'
        cvim_o.units = 'N/A'
        cvim_o.long_name = 'CVI mode flag'
        cvim_o.description = '0: CVI mode; 1: under-kinetic; -1: transition'
        df_o.units = 'N/A'
        df_o.long_name = 'CVI Dilution Factor'
        df_o.description = 'Dilution Factor after under-kinetic mode. Some measurements such as AMS, need to divide by this number'
        ef_o.units = 'N/A'
        ef_o.long_name = 'CVI Enhancement Factor'
        ef_o.description = 'Enhancement Factor after CVI mode. Some measurements such as AMS, need to divide by this number'
        merge_o.units = '#/cm3'
        merge_o.long_name = 'merged size distribution'
        fims_o.units = '#/cm3'
        fims_o.long_name = 'total aerosol concentration from FIMS'
        pcasp_o.units = '#/cm3'
        pcasp_o.long_name = 'total aerosol concentration from PCASP'
        
        # global attributes
        import time as ttt
        f.description = "Merged size distribution from FIMS and PCASP, in ambient condition"
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
