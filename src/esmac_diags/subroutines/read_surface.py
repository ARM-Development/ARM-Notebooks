"""
functions of reading surface data, typically in ascii format
"""
import numpy as np
#%%
# filename='../data/arm-ccn/HS_SGP_Nccn_data.dat'
def read_CCN_hiscale_IOP1(filepath):
    """
    read in surface CCN data for HISCALE field campaign in IOP1
    data provided by Jerome Fast and GK

    Parameters
    ----------
    filepath : string
        directory of the CCN file

    Returns
    -------
    time : time of data
    ccn : ccn data
    SS : super saturation for measured CCN
    timeunit : unit of time
    """
    filename = filepath+'HS_SGP_Nccn_data.dat'
        
    # nccn=12174  # the length of measurements. for checking purpose
    
    time = list()
    SS = list()
    ccn = list()
    f = open(filename, 'r')
    
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split()
        time.append(float(columns[0])/86400. + 127.65014)
        SS.append(float(columns[1]))
        ccn.append(float(columns[2]))
    f.close()
    timeunit = 'calendar day'
    return(time, ccn, SS, timeunit)

# filename = '../data/arm-ccn/N_CCN_corrected_IOP2.dat'
def read_CCN_hiscale_IOP2(filepath):
    """
    read in surface CCN data for HISCALE field campaign in IOP2
    data provided by Jerome Fast and GK

    Parameters
    ----------
    filepath : string
        directory of the CCN file

    Returns
    -------
    time : time of data
    ccn : ccn data
    SS : super saturation for measured CCN
    timeunit : unit of time
    """
    filename = filepath+'N_CCN_corrected_IOP2.dat'
        
    # nccn=38880  # the length of measurements. for checking purpose
    
    time = list()
    SS = list()
    ccn = list()
    f = open(filename, 'r')
    
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split()
        time.append(float(columns[0])/1440. + 241.)
        SS.append(float(columns[1]))
        ccn.append(float(columns[2]))
    f.close()
    timeunit = 'calendar day'
    return(time, ccn, SS, timeunit)

#%% PBLH from doppler lidar
# filename = '../data/dl-pblh/sgpdlC1_mlh_0.08.txt'
def read_dl_pblh(filename):
    """
    read in PBLH from doppler lidar
    data provided by Jerome Fast

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    data2 : PBL height in meter
    """
    f = open(filename, 'r')
    h = f.readline()
    h = f.readline()
    h = f.readline()
    h = h.strip()
    data = []
    if 'data2' in locals():
        del(data2)
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split()
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2 = np.asarray(source)
        else:
            data2 = np.column_stack((data2, source))
    f.close()
    return(data2)


#%% pnnl smps data

# for total information
# filename = '../data/pnnl-smps/HiScaleSMPSa_SGP_20160827_R0.ict.txt'
def read_smpsa_pnnl(filename):
    """
    read in surface SMPS data for HISCALE campaign

    Parameters
    ----------
    filename : string
        input filename. variables include:
        
        # dat_SMPS_UTC, second
        # Number, cm-3, particle number concentration
        # SurfaceArea, nm2cm-3, particle surface area
        # Volume, nm3cm-3, particle volume concentration
        # flag, NA, 0=good data, 1=zero period, 2=bad data

    Returns
    -------
    data2: SMPS aerosol number size distribution
    """
    f = open(filename, 'r')
    # read in data:
    h = 'aaa'
    while h[0:20] != 'dat_SMPS_UTC, Number':
        h = f.readline()
    h = h.strip()
    varlist = h.split(',')
    # data = []
    if 'data2' in locals():
        del(data2)
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2 = np.asarray(source)
        else:
            data2 = np.column_stack((data2, source))
    f.close()
    data2[data2 < -9990] = np.nan
    return(data2)

# filename = '../data/pnnl-smps/HiScaleSMPSb_SGP_20160827_R1.ict'
def read_smpsb_pnnl(filename):
    """
    read in surface SMPS data for HISCALE campaign, 
    format slight different with read_smpsa_pnnl

    Parameters
    ----------
    filename : string
        input filename. variables include:
        
        # dat_SMPS_UTC, second
        # Number, cm-3, particle number concentration
        # SurfaceArea, nm2cm-3, particle surface area
        # Volume, nm3cm-3, particle volume concentration
        # flag, NA, 0=good data, 1=zero period, 2=bad data

    Returns
    -------
    data2: SMPS aerosol number size distribution
    """
    f = open(filename, 'r')
    # read in data:
    h = 'aaa'
    while h[0:18] != 'dat_SMPS_UTC, NSD1':
        h = f.readline()
    h = h.strip()
    varlist = h.split(',')
    # data = []
    if 'data2' in locals():
        del(data2)
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2 = np.asarray(source)
        else:
            data2 = np.column_stack((data2, source))
    f.close()
    data2[data2 < -9990] = np.nan
    return(data2)

#%%
# filename = '../data/pnnl-smps/NSD_column_size_chart.txt'
def read_smps_bin(filename):
    """
    read in bin information of surface SMPS data for HISCALE campaign

    Parameters
    ----------
    filename : string
        input filename. 
        
    Returns
    -------
    dmean: mean diameter for each bin of SMPS aerosol number size distribution
    """
    f = open(filename, 'r')
    h = f.readline()
    h = h.strip()
    dmean = []
    for line in f:
        line = line.strip()  # remove \n
        columns = line.split('\t')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        dmean.append(source[1])
    f.close()
    return(dmean)

