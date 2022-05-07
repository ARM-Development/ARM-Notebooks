"""
all quality controls for data in ESMAC_Diags
examples of quality controls:
1. qc flags in ARM data
2. minimum and maximum cutoffs
3. systematic data corrections
4. additional data masking
"""

import numpy as np


#%% 
def qc_mask_cloudflag(data,cflag):
    """
    mask data with cloud flag, typically remove cloud_flag=1 

    Parameters
    ----------
    data : numpy array
        input data
    cflag : numpy array, int format
        cloud flag data. remove when cflag==1

    Returns
    -------
    data : output data

    """
    if len(data.shape)==1:
        data[cflag==1]=np.nan
    elif len(data.shape)==2:
        data[cflag==1,:]=np.nan
    else: 
        raise ValueError("dimension of input data and flag data are inconsistent")
    return(data)

#%% 
def qc_mask_qcflag(data,qc):
    """
    mask data with quality-control flags, typically remove all qc != 0

    Parameters
    ----------
    data : numpy array
        input data
    qc : numpy array, int format
        quality control flag data. remove when qc!=0

    Returns
    -------
    data : output data

    """
    if len(data.shape)==1:
        data[qc!=0]=np.nan
    elif len(data.shape)==2:
        data[qc!=0,:]=np.nan
    else: 
        raise ValueError("dimension of input data and flag data are inconsistent")
    return(data)

#%% 
def qc_mask_qcflag_cpc(data,qc):
    """
    mask data with quality-control flags, only for specific bit.
    
    this is only for CPC since QC in some bits may be good data and should be retained.
    according to Singh, Ashish <asingh@bnl.gov>, only the following qc_bit is critical 
     	C1/CPCf (I guess this is sgpaoscpcC1*) – bit_value 7, 8 
     	S1/CPCF – bit value 11, 12, 13, 14 
     	S1/CPCu or CPCuf- bit value 4, 5, 6, 15, 16  

    Parameters
    ----------
    data : numpy array
        input data
    qc : numpy array, int format
        quality control flag data.

    Returns
    -------
    data : output data

    """
    data_out = np.array(data)
    qc = qc.astype(np.int)
    for tt in range(len(qc)):
        if qc[tt] != 0:
            qc_bin = bin(qc[tt])
            qc_bit = []
            for ii in range(len(qc_bin)-2):
                if qc_bin[-1-ii]=='1':
                    qc_bit.append(ii+1)
            if any([x in qc_bit for x in [7,8,11,12,13,14]]):
                data_out[tt] = np.nan
    return(data_out)

def qc_mask_qcflag_cpcu(data,qc):
    data_out = np.array(data)
    qc = qc.astype(np.int)
    for tt in range(len(qc)):
        if qc[tt] != 0:
            qc_bin = bin(qc[tt])
            qc_bit = []
            for ii in range(len(qc_bin)-2):
                if qc_bin[-1-ii]=='1':
                    qc_bit.append(ii+1)
            if any([x in qc_bit for x in [4, 5, 6, 15, 16]]):
                data_out[tt] = np.nan
    return(data_out)

#%% 
def qc_mask_takeoff_landing(time,data):
    """
    mask aircraft takeoff/landing time
    mask time is set as 30min to exclude possible land contamination for ocean measurements

    Parameters
    ----------
    time : measurement time in seconds
    data : input data
    time is only 1-dimensional while data can be 1- or 2-dimensional

    Returns
    -------
    masked data

    """
    # time should be in unit of seconds
    idx=np.logical_or(time<(time[0]+1800), time>(time[-1]-1800))
    if len(data.shape)==1:
        data[idx]=np.nan
    elif len(data.shape)==2:
        data[:,idx]=np.nan
    else:
        raise ValueError("input data can only be 1-d or 2-d array")
    return(data)

#%% 
def qc_remove_neg(data, remove_zero='False'):
    """
    remove negative values
    options of keep or remove zero value
    """
    if remove_zero == 'False' or remove_zero == 'false':
        data[data<0]=np.nan
    elif remove_zero == 'True' or remove_zero == 'true':
        data[data<=0]=np.nan
    else:
        raise ValueError("remove_zero can only be true or false")
    return(data)
    
#%%
def qc_ccn_max(ccn, SS):
    """
    set a max value for CCN measurement
    different maximum threshold for different supersaturations
    SS can be a fixed value of changable for different ccn values
    
    Parameters
    ----------
    ccn : ccn number concentration
    SS : supersaturation
    ccn and SS should be the same size

    Returns
    -------
    ccn with quality controls

    """
    ccn[np.logical_and(SS<0.2, ccn>2000)] = np.nan
    ccn[np.logical_and(SS<0.6, ccn>4000)] = np.nan
    ccn[ccn>8000] = np.nan
    return(ccn)
    
#%%
def qc_cn_max(data, size):
    """
    set a max value for CN measurements
    different threshold for different minimal cufoff size
    
    Parameters
    ----------
    data : input aerosol number concentration
    size : minimum detected aerosol size (nm) in CPC

    Returns
    -------
    data with value greater than maximum removed

    """
    if size == 3:
        data[data>5e4] = np.nan
    elif size == 10:
        data[data>2.5e4] = np.nan
    elif size == 100:
        data[data>0.5e4] = np.nan
    return(data)

#%%
def qc_cpc_air(cpc3, cpc10):
    """
    set a minimum threshold for ARM G1 flight CPC measurements
    """
    cpc3[cpc3<20] = np.nan
    cpc10[cpc10<10] = np.nan
    return(cpc3, cpc10)

#%%
def qc_correction_nanosmps(data):
    """
    nanoSMPS at SGP is known to overcount particle number. make a correction to ensure smooth transition with SMPS
    the fraction of 3.8 is chosen so that nanoSMPS and SMPS is smoothly transit at ~18 nm
    """
    data = data/3.8
    return(data)

#%%
def qc_fims_bin(data):
    """
    remove some suspicious value of FIMS measurements
    """
    data2 = np.array(data)
    # data2[:, data2[0, :] > 3e4] = np.nan
    data2[np.logical_or(data2 < 0, data2 > 3e4)] = np.nan
    return(data2)

#%%
def qc_acsm_org_max(data):
    """
    set max value for ACSM measured total organic matters
    """
    data[data>10] = np.nan # max value set as 10 ug/m3
    return(data)

#%%
def qc_uhsas_RF_NCAR(data):
    """
    set a maximum value for NCAR research flight UHSAS measurements
    """
    
    data[data>500] = np.nan
    return(data)