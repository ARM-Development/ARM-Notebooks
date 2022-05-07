"""
functions of reading ARM standard NETCDF files
"""
import numpy as np
from netCDF4 import Dataset
#%%
# filename='../data/arm-cpcu/iop1/sgpaoscpcuS01.b1.20160422.000000.nc'
def read_acsm(filename, varname):
    """
    read in ACSM data

    Parameters
    ----------
    filename : string
        input filename
    varname : string
        variable name to be read

    Returns
    -------
    time : time of measurments
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    """
    
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    # flag = f.variables['qc_'+varname][:]
    # data[flag!=0]=np.nan
    f.close()
    return(time, data, timeunit, dataunit)

#%%
# filename='../../data/ACEENA/obs/profile/enaarmbecldrad/enaarmbecldradC1.c1.20170101.003000.nc'
# varname='cld_frac'
def read_armbe(filename, varname):
    """
    read in ARMBE data

    Parameters
    ----------
    filename : string
        input filename
    varname : string
        variable name to be read

    Returns
    -------
    time : time of measurments
    height : vertical dimension in meter
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    """
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    height = f.variables['height'][:]
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    # flag = f.variables['qc_'+varname][:]
    # data[flag!=0]=np.nan
    f.close()
    return(time, height, data, timeunit, dataunit)

#%%
# filename='../data/ccn_aaf/enaaafccn2colaF1.b1.20180121.094335.nc'
def read_ccn(filename):
    """
    read in CCN data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    timeunit : unit of time
    ccn : ccn data
    qc_ccn : quality control flag for ccn data
    dataunit : unit of ccn
    SS : super saturation for measured CCN
    """ 

    f = Dataset(filename, 'r')
    # read in variables
    # lon = f.variables['lon'][:]
    # lat = f.variables['lat'][:]
    # alt = f.variables['alt'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    SS = f.variables['supersaturation_calculated'][:]
    d_id = f.variables['N_CCN']
    qc_ccn = f.variables['qc_N_CCN'][:]
    ccn = d_id[:]
    dataunit = d_id.units
    f.close()
    # ccn[qc_ccn!=0] = -9999.
    return(time, timeunit, ccn, qc_ccn, dataunit, SS)

# MAGIC
def read_ccn_magic(filename):
    """
    read in CCN data for MAGIC field campaign, no qc_ccn

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    timeunit : unit of time
    ccn : ccn data
    dataunit : unit of ccn
    SS : super saturation for measured CCN
    """ 

    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    SS = f.variables['CCN_ss_set'][:]
    d_id = f.variables['N_CCN']
    ccn = d_id[:]
    dataunit = d_id.units
    f.close()
    return(time, timeunit, ccn, dataunit, SS)

#%%
# filename='../data/arm-cpcu/iop1/sgpaoscpcuS01.b1.20160422.000000.nc'
def read_cpc(filename):
    """
    read in CPC data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    data : data of the variable
    qc_flag : quality control flag for data
    timeunit : unit of time
    dataunit : unit of data
    """
    
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables['concentration']
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    try:
        qc_flag = f.variables['qc_concentration'][:]
    except: # pylint: disable=bare-except
        qc_flag = np.full(len(data), 0.0)
    f.close()
    return(time, data, qc_flag, timeunit, dataunit)

#%%
# filename='../data/inletcvi/enaaafinletcviF1.c1.20170718.083145.nc'
def read_cvi_aceena(filename):  
    """
    read in CVI data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    lon : longitude
    lat : latitude
    alt : altitude
    timeunit : unit of time
    cvimode : CVI mode in the data
    cviinlet : CVI inlet status
    enhance_factor : enhancement factor
    dilution_factor : dilution factor
    
    see ARM CVI handbook for more information
    """  

    f = Dataset(filename, 'r')
    # read in variables
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    alt = f.variables['alt'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables['cvi_mode']
    cvimode = d_id[:]
    d_id = f.variables['inlet_selector']
    cviinlet = d_id[:]
    d_id = f.variables['enhancement_factor']
    enhance_factor = d_id[:]
    d_id = f.variables['inlet_dilution_factor']
    dilution_factor = d_id[:]
    f.close()
    return(time, lon, lat, alt, timeunit, cvimode, cviinlet, enhance_factor, dilution_factor)

#%%
# filename='../data/arm-met/iop1/sgpmetE13.b1.20160502.000000.cdf'
# varname='wdir_vec_mean'
def read_met(filename, varname):
    """
    read in met data

    Parameters
    ----------
    filename : string
        input filename
    varname : string
        variable name to be read

    Returns
    -------
    time : time of measurments
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : longname of data
    """
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    f.close()
    return(time, data, timeunit, dataunit, long_name)

#%% mwrret 
def read_mwr(filename, varname):
    """
    read in microwave radiometer (MWR) data

    Parameters
    ----------
    filename : string
        input filename
    varname : string
        variable name to be read

    Returns
    -------
    time : time of measurments
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    flag : quality control flag
    """
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    # long_name = d_id.long_name
    try:
        flag = f.variables['qc_'+varname][:]
    except: # pylint: disable=bare-except
        flag = data*0
    f.close()
    return(time, data, timeunit, dataunit, flag)

#%%
# filename='../data/arm-pblh/sgppblhtsonde1mcfarlC1.c1.20160501.052900.cdf'
def read_pblhtsonde1(filename):
    """
    read in PBLH from radiosonde data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    timeunit : unit of time
    height2 : height of measurements
    p : pressure
    T : temperature
    rh : relative humidity
    wspd : wind speed
    pblh : PBL height
    pbl_type : different types to estimate PBLH
    long_name : long name of PBLH for different types
    """
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    height = f.variables['height_ss'][:]
    T = f.variables['air_temp'][:]
    p = f.variables['atm_pres'][:]
    p2 = f.variables['pressure_gridded'][:]
    rh = f.variables['rh'][:]
    wspd = f.variables['wspd'][:]
    pbl_type = ['heffter', 'liu_liang', 'bulk_richardson_pt25', 'bulk_richardson_pt5']
    pblh = np.nan*np.zeros(len(pbl_type))
    long_name=list()
    for i in range(len(pbl_type)):
        d_id = f.variables['pbl_height_'+pbl_type[i]]
        qc = f.variables['qc_pbl_height_'+pbl_type[i]][:]
        if qc==0:
            pblh[i] = d_id[:]
        long_name.append(d_id.long_name)
    f.close()
    
    height2=np.interp(p, np.flip(p2), np.flip(height))
    return(time, timeunit, height2, p, T, rh, wspd, pblh, pbl_type, long_name)

#%%
# filename='../data/arm-pblh/sgppblhtmpl1sawyerliC1.c1.20160830.000002.nc'
def read_pblhtmpl1(filename):
    """
    read in PBLH from MPL data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    timeunit : unit of time
    pblh : PBL height
    qc : quality control flag
    """
    f = Dataset(filename, 'r')
    # read in variables
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    # height = f.variables['height'][:]
    d_id = f.variables['annealing_pbl_height_sawyer_li']
    qc = f.variables['qc_annealing_pbl_height_sawyer_li'][:]
    pblh = d_id[:]
    # pblh[qc!=0] = np.nan
    f.close()
    return(time, timeunit, pblh, qc)

#%%
# filename='../data/bnl-smps/sgpaossmpsS01.a1.20160513.000000.nc'
# varname='total_concentration'

def read_smps_bnl(filename, varname):  
    """
    read in SMPS data processed by BNL

    Parameters
    ----------
    filename : string
        input filename
    varname : string
        variable name to be read

    Returns
    -------
    time : time of measurments
    diameter : diameter for aerosol size bins
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : data long names
    """  

    f = Dataset(filename, 'r')
    # read in variables
    diameter = f.variables['diameter_midpoint'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    f.close()
    return(time, diameter, data, timeunit, dataunit, long_name)



#%%
# filename='../data/arm-uhsas/sgpaosuhsasS01.a1.20160501.000003.nc'
def read_uhsas(filename):    
    """
    read in UHSAS data

    Parameters
    ----------
    filename : string
        input filename

    Returns
    -------
    time : time of measurments
    dmin : lower bound of aerosol size bins
    dmax : upper bound of aerosol size bins
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : data long names
    """

    f = Dataset(filename, 'r')
    # read in variables
    dmin = f.variables['lower_size_limit'][:]
    dmax = f.variables['upper_size_limit'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    try:        # direct read concentration
        d_id = f.variables['concentration']
        data = d_id[:]
        dataunit = d_id.units
        long_name = d_id.long_name
    except: # pylint: disable=bare-except
         # calculate from raw count
        raw_count = f.variables['size_distribution'][:]
        flow_rate = f.variables['sampling_volume'][:]/60.    # cc/min to cc/s
        sample_time = 10    # sample interval is 10s
        data=np.full(raw_count.shape, np.nan)
        for bb in range(data.shape[1]):
            data[:, bb] = raw_count[:, bb] /flow_rate /sample_time
        dataunit='1/cc'
        long_name='size distribution'
    f.close()
    return(time, dmin, dmax, data, timeunit, dataunit, long_name)
