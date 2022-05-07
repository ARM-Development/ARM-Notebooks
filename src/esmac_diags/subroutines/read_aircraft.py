"""
functions of reading aircraft data, mostly ascii format
"""
import numpy as np
from netCDF4 import Dataset
from .time_format_change import hhmmss2sec

#%% 
# filename='../../data/HiScale/obs/aircraft/shilling-ams\\HiScaleAMS_G1_20160425_R0.ict'
def read_ams(filename):    
    """
    READ AMS composition data

    Parameters
    ----------
    filename : str
        input filename

    Returns
    -------
    data2: ams measured variables
    varlist: list of variable names
    """
        
    f = open(filename, 'r')
    
    # read in data:
    
    h = 'aaa'
    for ii in range(40):
        h = f.readline()
    while(h[0:7]!='dat_ams' and h[0:7]!='dat_utc'):
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2, varlist)

#%% 
def read_beasd(filename):    
    """
    read best estimate aerosol size distribution data

    Parameters
    ----------
    filename : str
        input filename

    Returns
    -------
    data2: beasd variables
    varlist: list of variable names
    size_h : upper bound of size bin
    size_l : lower bound of size bin
    size_m : middle value of size bin

    """    
    
    f=open(filename,'r')
    
    # read in data:
    h='aaa'
    while h[0:24]!='UPPER_BIN_SIZE_nanometer':
        h=f.readline()
    h=h.strip()
    size_h=h.split(',')
    h=f.readline()
    h=h.strip()
    size_l=h.split(',')
    h=f.readline()
    h=h.strip()
    size_m=h.split(',')
    
    size_h = np.array([np.float32(x.split(' ')[-1]) for x in size_h])
    size_m = np.array([np.float32(x.split(' ')[-1]) for x in size_m])
    size_l = np.array([np.float32(x.split(' ')[-1]) for x in size_l])
    
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0,len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if 'data2' not in locals():
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))    
    f.close()
    
    return(data2,varlist, size_h, size_l, size_m)


#%% 
# filename='../data/mei-cpc/CPC_G1_20160830143515_R2_HiScale001s.ict.txt'
def read_cpc(filename): 
    """
    READ CPC data (from mei-cpc)

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """    
        
    f=open(filename, 'r')
    
    # read in data:
    
    h='aaa'
    while h[0:19]!='Start_UTC,CPC_Conc_':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    data2[data2<-9990]=np.nan
    return(data2, varlist)

#%%
# filename='../data/mei-ccn/CCN_G1_20160518170015_R2_HiScale001s.ict'
def read_ccn_hiscale(filename):
    """
    READ CCN data (from mei-ccn)

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """    
        
    f=open(filename, 'r')
    
    # read in data:
    
    h='aaa'
    while h[0:14]!='Start_UTC,DT_A':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2, varlist)

#%%
def read_ccn_socrates(filename):
    """
    READ CPC data for SOCRATES

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """    
        
    f=open(filename, 'r')
    
    # read in data:
    
    h='aaa'
    # for CCNscanning or CCNspectra data
    while h[0:14] not in ['Start_UTC, CCN', 'Start_UTC, Sto']:
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2, varlist)

#%%
# filename='../data/pekour-cvi/CVI_G1_20160518170015_R4_HISCALE_001s.ict.txt'
def read_cvi_hiscale(filename):  
    """
    READ in CVI data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """      
    f=open(filename, 'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2, varlist)

#%%
def read_fims(filename):
    """
    READ in wang-fims data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """
    
    f=open(filename, 'r', errors='replace')
    
    # read in data:
    
    h='aaa'
    while h[0:10]!='UTC,n_Dp_1':
        h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    
    # data2[data2<-9990]=np.nan
    
    # multiply dlnDp since FIMS data are dN/dlnDp 
    # data2[1:31, :]=data2[1:31, :]*0.1272
    return(data2, varlist)


#%%
# filename='../data/wang-fims/HISCALE_FIMS_bins_R1.dat'
def read_fims_bin(filename):
    """
    READ in bin information of wang-fims data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    dmean: mean aerosol diameter for each bin
    dmin: lower bound of diameter
    dmax: upper bound of diameter
    """
    f=open(filename, 'r')
    h=f.readline()
    h=h.strip()
    dmean=[]
    dmin=[]
    dmax=[]
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        dmin.append(source[0])
        dmean.append(source[1])
        dmax.append(source[2])
    f.close()

    return(dmean, dmin, dmax)

#%%
# filename='../data/mei-iwg1/aaf.iwg1001s.g1.hiscale.20160511a.a2.txt'
def read_iwg1(filename):
    """
    READ in mei-iwg data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data: all measured variables
    a: list of variable names and units
    """
    
    
    f=open(filename, 'r')
    varname=f.readline()
    varunit=f.readline()
    varname=varname.strip()
    varname=varname.split(',')
    varunit=varunit.strip()
    varunit=varunit.split(',')
    a = np.column_stack((np.asarray(varname), np.asarray(varunit)))
    data=[]
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(columns[i])
        data.append(source)
    
    f.close()
    # data[data<-9990]=np.nan
    return(data, a)
    
# except:
#     print('error')

#%%
# filename='../data/Kappa/FIMS_kappa_IOP2_part2/20160830aAir_data_table_FIMS_kappa_col_A.dat'
def read_kappa(filename):    
    """
    READ in kappa data. 
    this data is calculated from FIMS and CCN measurements.
    currently it is not used in the ESMAC Diags

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varname: list of variable names
    """
    f=open(filename, 'r')
    varname=f.readline()
    # varunit=f.readline()
    varname=varname.strip()
    varname=varname.split()
    # data=[]
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split()
        if columns==[]:
            continue
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    return(data2, varname)

#%% read mergedSD data
# filename = 'C:/Users/tang357/OneDrive - PNNL/EAGLES/python_diag_pkg/ESMAC_Diags_Tool/'+\
#     'data/HISCALE/obs/aircraft/mergedSD/aaf.g1.hiscale.mergedSD.20160903a.txt'
def read_mergedSD(filename):   
    """
    READ in mergedSD data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    time2: time in seconds
    n_total: total droplet number concentration
    n_bins: measured cloud droplet number concentration in each size bin
    d_min: lower bound of each size bins
    d_max: upper bound of each size bins
    d_mean: center diameter for each size bin
    """ 

    f=open(filename, 'r')
    varname=f.readline()
    # varunit=f.readline()
    varname=varname.strip()
    varname=varname.split(',')
    varname=varname[1:]
    # get droplet size information
    dmin=[]
    dmax=[]
    dmean=[]
    for i in range(1,len(varname)):
        d_str = varname[i].split('_')
        d1 = np.float64(d_str[1].replace('p','.'))
        d2 = np.float64(d_str[2].replace('p','.'))
        dmean.append((d1+d2)/2)
        dmin.append(d1)
        dmax.append(d2)
    dmin = np.array(dmin)
    dmax = np.array(dmax)
    dmean = np.array(dmean)
    # read in time and data
    date=[]
    time2=np.array([])
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        if line==[]:
            continue
        columns = line.split(',')
        # print(columns[0])
        time = columns[0].split()
        date.append(time[0])
        time2 = np.hstack((time2, hhmmss2sec(time[1])))
        source = []
        for i in range(1, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    
    f.close()
    n_total = data2[0,:]
    n_bins = data2[1:, :]
    return(time2, n_total, n_bins, dmean, dmin, dmax)


#%% read OPC
# filename='../data/opciso/OPCISO_G1_20170707103233_R3_ACEENA_001s.ict'
def read_opc(filename):   
    """
    READ in OPC data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: measured aerosol number concentration in each size bin
    d_min: lower bound of aerosol size bins
    d_max: upper bound of aerosol size bins
    d_center: center diameter for each size bin
    varlist: list of variable names in data2
    """ 
    f=open(filename, 'r')
    h='aaa'
    while h[0:26]!='UPPER_BIN_SIZE_micrometer:':
        h=f.readline()
    h=h.strip()
    d_max = h[26:].split(',')
    d_max=[float(x) for x in d_max]
    h=f.readline()
    h=h.strip()
    d_min = h[26:].split(',')
    d_min=[float(x) for x in d_min]
    h=f.readline()
    h=h.strip()
    d_center = h[24:].split(',')
    d_center=[float(x) for x in d_center]
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    # data2[data2<-9990]=np.nan
    return(data2, np.array(d_min), np.array(d_max), np.array(d_center), varlist)


#%%
# filename='../data/tomlinson-pcasp/pcasp_g1_20160511163038_R2_L1_hiscale001s.ict.txt'
def read_pcasp(filename):    
    """
    READ in PCASP aerosol size distribution data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """
    f=open(filename, 'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    data2[data2<-9990]=np.nan
    return(data2, varlist)
    
#%% read research flight data from NCAR
def read_RF_NCAR(filename, varname):
    """
    READ in NCAR research flight data

    Parameters
    ----------
    filename : input filename
    varname : choose variables to read

    Returns
    -------
    time: time of measurements
    data: measured variables
    timeunit: description of time
    dataunit: unit of data
    long_name: long_name or description of the data
    cellsize: if read in size distribution, return size information of each cell
    cellunit: unit of cell size
    
    cellsize and cellunit are 'N/A' if variable is not size distribution
    """
    f = Dataset(filename, 'r')
    
    # read in variables
    t_id = f.variables['Time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    try:
        cellsize = d_id.CellSizes
    except: # pylint: disable=bare-except
        cellsize = 'N/A'
    try:
        cellunit = d_id.CellSizeUnits
    except: # pylint: disable=bare-except
        cellunit = 'N/A'

    
    f.close()
    
    return(time, data, timeunit, dataunit, long_name, cellsize, cellunit)

#%%
# filename='../data/springston-tracegases/SO2/SO2_G1_20160510_R2.ict'
def read_so2(filename):
    """
    READ in SO2 data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """
    f=open(filename, 'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    data2[data2<-9990]=np.nan
    return(data2, varlist)
    
#%%
# filename='../data/tomlinson-uhsas/uhsasa_g1_20160425155810_R1_L1_hiscale001s.ict.txt'
def read_uhsas(filename):  
    """
    READ in UHSAS data

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """  
    f=open(filename, 'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    data2[data2<-9990]=np.nan
    return(data2, varlist)

#%% read LWC from WCM
# filename='../data/matthews-wcm/WCM_G1_20160909150045_R2_HISCALE001s.ict.txt'
def read_wcm(filename):
    """
    READ in LWC data from WCM measurements

    Parameters
    ----------
    filename : input filename

    Returns
    -------
    data2: all measured variables
    varlist: list of variable names
    """
    f=open(filename, 'r')
    h='aaa'
    while h[0:3]!='R0:':
        h=f.readline()
    h=f.readline()
    h=h.strip()
    varlist=h.split(',')
    if 'data2' in locals():
        del(data2)
    for line in f:
        line=line.strip()  # remove \n
        columns = line.split(',')
        source = []
        for i in range(0, len(columns)):
            source.append(float(columns[i]))
        # data.append(source)
        if not('data2' in locals()):
            data2=np.asarray(source)
        else:
            data2=np.column_stack((data2, source))
    f.close()
    data2[data2<-9990]=np.nan
    return(data2, varlist)


# In[test read in data]
# filename = '../data/wang-fims/FIMS_G1_20160830_R1_L1_HISCALE_001s.ict'

# (x, y)=read_fims(filename)

# filename='../data/mei-iwg1/aaf.iwg1001s.g1.hiscale.20160511a.a2.txt'
# (x, y)=read_iwg1(filename)