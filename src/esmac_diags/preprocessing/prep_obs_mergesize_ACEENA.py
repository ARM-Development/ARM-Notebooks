"""
# merge size distribution from FIMS, PCASP and OPC for ACE-ENA
# revised from size_bin_merge.pro by Jerome Fast
# Shuaiqi Tang
# 2020.10.1
"""

import os
import glob
import re
import numpy as np
from ..subroutines.read_aircraft import read_fims, read_fims_bin, read_iwg1, read_pcasp, read_opc
from ..subroutines.read_ARMdata import read_cvi_aceena as read_cvi
from ..subroutines.time_format_change import hhmmss2sec
from netCDF4 import Dataset

def run_prep(settings):
    #%% variables from settings
    iwgpath = settings['iwgpath']
    fimspath = settings['fimspath']
    pcasppath = settings['pcasppath']
    opcpath = settings['opcpath']
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
    (d_fims, dmin_f, dmax_f) = read_fims_bin(fimspath + 'ACEENA_FIMS_bins_R0.dat')
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
        fname = re.split('aceena.|.a2', filename)
        date = fname[-2]
        print(date)
        
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
        datestr = timestr[0]
    
        # FIMS
        filename_f = glob.glob(fimspath + 'FIMS_G1_' + date[0:8] + '*_001s.ict')
        # read in data
        if len(filename_f) == 1:
            (data0, fimslist) = read_fims(filename_f[0])
            time_fims = data0[0, :]
            # change data from #/dlnDp to number
            data2 = data0[1:-3, :] * dlnDp_f
            # TD mode or AMB mode. remove TD mode
            TD_AMB = data0[-1, :]
            data2[:, TD_AMB != 0] = -9999.            
            T2 = data0[-2,:]
            p2 = data0[-3,:]*1000.
            T2 = np.interp(time,time_fims,T2)
            p2 = np.interp(time,time_fims,p2)
            fims = np.empty([30, len(time)])
            for ii in range(30):
                fims[ii, :] = np.interp(time, time_fims, data2[ii, :])
            idx = np.logical_or(time>time_fims[-1], time<time_fims[0])
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
        filename_p = glob.glob(pcasppath + 'pcasp_g1_' + date[0:8] + '*' + '_aceena001s.ict')
    
        binlen = 30
        dmax_p = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200]
        dmin_p = [100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 280, 300, \
                400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
        dmin_p = [x/1000 for x in dmin_p]
        dmax_p = [x/1000 for x in dmax_p]
        # read in data
        if len(filename_p) == 1:
            (data0, pcasplist) = read_pcasp(filename_p[0])
            pcasp2 = data0[1:-5, :]
            time_pcasp = data0[0, :]
            d_pcasp = [float(i) for i in pcasplist[1:-5]]
            pcasp = data0[1:-5, :]
            flag = data0[-2, :]
            pcasp_total = data0[-5, :]
            # remove some questionable data
            # pcasp[np.isnan(pcasp)] = -9999
            # pcasp[np.logical_or(pcasp <=0, pcasp>1e6)] = np.nan
            # pcasp[:, flag != 0] = np.nan
            if not all(time_pcasp == time):
                raise ValueError('PCASP time is inconsistent with FIMS')
        elif len(filename_p) == 0:
            time_pcasp = time
            d_pcasp = [(dmin_p[x] + dmax_p[x])/2 for x in range(len(dmin_p))]
            pcasp = np.full([len(d_pcasp), len(time)], -9999.)
            pcasp_total = np.full(len(time), -9999.)
        else:
            raise ValueError('find more than one file: ' + filename_p)
    
        # !! PCASP data is for standard T and p (Conc = Conc_orig*[(1013.25/Pamb)*(Tamb/293.15)]), change to ambient T/p
        pcasp2 = np.array(pcasp)
        for tt in range(len(time)):
            pcasp[:, tt] = pcasp[:, tt]/((1013.25/p_amb[tt])*((T_amb[tt] + 273.15)/293.15))
        
        # CVI
        filename_c = glob.glob(cvipath + 'enaaafinletcviF1.c1.' + date[0:8] + '*.nc')
        filename_c.sort()
        # read in data
        if len(filename_c) == 1:
            (time_c, lon_c, lat_c, alt_c, timeunit_c, cvimode, cvi_inlet, enhance_factor, dilution_factor) = read_cvi(filename_c[0])
            if date == '20180216a':
                time_c = np.insert(time_c, 1403, (time_c[1402] + time_c[1403])/2)
                cvi_inlet=np.insert(cvi_inlet, 1403, cvi_inlet[1403])
                cvimode=np.insert(cvimode, 1403, cvimode[1403])
                enhance_factor = np.insert(enhance_factor, 1403, enhance_factor[1403])
                dilution_factor = np.insert(dilution_factor, 1403, dilution_factor[1403])
            enhance_factor[enhance_factor<-9000] = np.nan
            dilution_factor[dilution_factor<-9000] = np.nan
            if not all(time_c == time):
                raise ValueError('CVI time is inconsistent with FIMS')
        elif len(filename_c) == 0:
            time_c = time
            cvi_inlet = np.nan*np.empty([len(time)])
            cvimode = np.nan*np.empty([len(time)])
            enhance_factor = np.nan*np.empty([len(time)])
            dilution_factor = np.nan*np.empty([len(time)])
        else:
            raise ValueError('find more than one file: ' + filename_c)
    
        cvi_inlet[cvi_inlet == -9] = 1  # if cvi_inlet is unfunctional, assume it is isokinetic and use fims as good data
        
        
        # read OPC
        filename_o = glob.glob(opcpath + 'OPCISO_G1_' + date[0:8] + '*.ict')
        if len(filename_o) == 1:
            (opc, dmin_o, dmax_o, d_opc, opclist) = read_opc(filename_o[0])
            time_o = opc[0, :]
            opc = opc[1:, :]
            opc[opc<0] = np.nan
        else:
            raise ValueError('can not find OPC data or find multiple files: ' + filename_o)
        if date == '20180216a':
            time_o = np.hstack((time_o[0:1403], 47873., time_o[1403:]))
            opc = np.hstack((opc[:, 0:1403], (opc[:, 1402:1403] + opc[:, 1403:1404])/2, opc[:, 1403:]))
        if any(time_o != time):
            raise ValueError('OPC time is inconsistent with FIMS')
        if sum(np.isnan(opc[0, :]))<0.1*len(time_o):
            for ii in range(len(d_opc)):    # fill missing timesteps
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
        idx = dmax_p.index(0.8)
        dia_merge_l[nbin_fims] = dmax_f[-1]
        dia_merge_h[nbin_fims] = dmax_p[idx]
        # next bin uses 0.8 as low bound and high bound of 2nd bin (0.9) of OPC
        dia_merge_l[31] = 0.8
        dia_merge_h[31] = 0.9
        # next few bins are merged two OPC bins
        for n in range(1, 6):
            dia_merge_l[31 + n] = dmin_o[n*2]
            dia_merge_h[31 + n] = dmax_o[n*2 + 1]
        # other bins follows OPC bins
        for n in range(12, nbin_opc):
            dia_merge_l[25 + n] = dmin_o[n]
            dia_merge_h[25 + n] = dmax_o[n]
            
        d_merge = (dia_merge_h + dia_merge_l)/2
        
        # merged concentration
        conc_merge = np.full([timelen, nbin_merge], -9999.)
        fims[np.isnan(fims)] = -9999.   # do not treat missing as NaN. treat -9999
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
            idx = dmin_p.index(0.3)   # start merging size. choose the index of pcasp for merging
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
            if np.isnan(opc2):             
                conc_merge[k, 31] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 31] = opc2
            else:
                conc_merge[k, 31] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[18, k]*0.5 + pcasp[19, k]*0.2
            opc2 = opc[2, k] + opc[3, k]
            if np.isnan(opc2):             
                conc_merge[k, 32] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 32] = opc2
            else:
                conc_merge[k, 32] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[19, k]*0.8
            opc2 = opc[4, k] + opc[5, k]
            if np.isnan(opc2):             
                conc_merge[k, 33] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 33] = opc2
            else:
                conc_merge[k, 33] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.9
            opc2 = opc[6, k] + opc[7, k]
            if np.isnan(opc2):             
                conc_merge[k, 34] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 34] = opc2
            else:
                conc_merge[k, 34] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[20, k]*0.1 + pcasp[21, k]
            opc2 = opc[8, k] + opc[9, k]
            if np.isnan(opc2):             
                conc_merge[k, 35] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 35] = opc2
            else:
                conc_merge[k, 35] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[22, k] + pcasp[23, k]*0.2
            opc2 = opc[10, k] + opc[11, k]
            if np.isnan(opc2):             
                conc_merge[k, 36] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 36] = opc2
            else:
                conc_merge[k, 36] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.7
            opc2 = opc[12, k]
            if np.isnan(opc2):             
                conc_merge[k, 37] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 37] = opc2
            else:
                conc_merge[k, 37] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[23, k]*0.1 + pcasp[24, k]*0.7
            opc2 = opc[13, k]
            if np.isnan(opc2):             
                conc_merge[k, 38] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 38] = opc2
            else:
                conc_merge[k, 38] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[24, k]*0.3 + pcasp[25, k]*0.4
            opc2 = opc[14, k]
            if np.isnan(opc2):             
                conc_merge[k, 39] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 39] = opc2
            else:
                conc_merge[k, 39] = (opc2 + pcasp2)/2.0    # assume equal weight
            pcasp2 = pcasp[25, k]*0.6 + pcasp[26, k]*0.3
            opc2 = opc[15, k]
            if np.isnan(opc2):             
                conc_merge[k, 40] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 40] = opc2
            else:
                conc_merge[k, 40] = opc2*0.6 + pcasp2*0.4  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[26, k]*0.7 + pcasp[27, k]*0.2
            opc2 = opc[16, k]
            if np.isnan(opc2):             
                conc_merge[k, 41] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 41] = opc2
            else:
                conc_merge[k, 41] = opc2*0.7 + pcasp2*0.3  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[27, k]*0.8 + pcasp[28, k]*0.2
            opc2 = opc[17, k]
            if np.isnan(opc2):             
                conc_merge[k, 42] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 42] = opc2
            else:
                conc_merge[k, 42] = opc2*0.8 + pcasp2*0.2  # gradually reduce the weight of PCASP
            pcasp2 = pcasp[28, k]*0.8 + pcasp[29, k]*0.3
            opc2 = opc[18, k]
            if np.isnan(opc2):             
                conc_merge[k, 43] = pcasp2
            elif pcasp2>10*opc2 or pcasp2<0:
                conc_merge[k, 43] = opc2
            else:
                conc_merge[k, 43] = opc2*0.9 + pcasp2*0.1  # gradually reduce the weight of PCASP
            # using OPC for other bins
            for n in range(44, nbin_merge):
                conc_merge[k, n] = opc[n-25, k]
            
        
        #%% output data
        if not os.path.exists(merged_size_path):
            os.mkdir(merged_size_path)
        outfile = merged_size_path + 'merged_bin_fims_pcasp_opc_ACEENA_' + date + '.nc'
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
        cvim_o[:] = np.array(cvimode)
        dilution_factor[np.isnan(dilution_factor)] = -9999.
        df_o[:] = dilution_factor
        enhance_factor[np.isnan(enhance_factor)] = -9999.
        ef_o[:] = enhance_factor
        conc_merge[np.isnan(conc_merge)] = -9999.
        conc_merge[conc_merge<0] = -9999.
        merge_o[:, :] = conc_merge
        fims_total[np.isnan(fims_total)] = -9999.
        fims_total[fims_total<0] = -9999.
        fims_o[:] = fims_total
        pcasp_total[np.isnan(pcasp_total)] = -9999.
        pcasp_total[pcasp_total<0] = -9999.
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
        f.description = "Merged size distribution from FIMS, PCASP and OPC"
        f.create_time = ttt.ctime(ttt.time())
        
        f.close()
        
    
