"""
functions of reading some NETCDF files 
"""
from netCDF4 import Dataset

#%%
# filename='../data/MAM/cutoff_number_MAM5.nc'
# varname='T'
def read_E3SM_z(filename, varname):
    """
    read variables from E3SM output with height variable

    Parameters
    ----------
    filename : string
        filename
    varname : string
        variable name

    Returns
    -------
    time : time of data
    height : vertical dimension in meter
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : longname of data
    """
    
    
    f = Dataset(filename, 'r')
    
    # read in variables
    height = f.variables['height'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    
    return(time, height, data, timeunit, dataunit, long_name)


#%%
def read_E3SM(filename, varname):
    """
    read variables from E3SM output

    Parameters
    ----------
    filename : string
        filename
    varname : string
        variable name

    Returns
    -------
    time : time of data
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
    if type(varname) is str:
        d_id = f.variables[varname]
        data = d_id[:]
        try:
            dataunit = d_id.units
        except:
            dataunit='N/A'    
        long_name = d_id.long_name
    elif type(varname) is list:
        data=list()
        dataunit=list()
        long_name=list()
        for vv in range(len(varname)):
            d_id=f.variables[varname[vv]]
            data.append(d_id[:])
            try:
                dataunit.append(d_id.units)
            except:
                dataunit.append('N/A')    
            long_name.append(d_id.long_name)
    
    f.close()
    
    return(time, data, timeunit, dataunit, long_name)

#%%
# filename='../data/MAM/NCN_MAM5_20160907_R1_L2.nc'
# varname='NCN'
def read_extractflight(filename, varname):
    """
    read variables from model-extracted files along flight tracks

    Parameters
    ----------
    filename : string
        filename
    varname : string
        variable name

    Returns
    -------
    time : time of data
    height : vertical dimension in meter
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : longname of data
    """
    f = Dataset(filename, 'r')
    # read in variables
    height = f.variables['height'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    return(time, height, data, timeunit, dataunit, long_name)


#%% read merged size distribution data
def read_merged_size(filename, varname):
    """
    read merged aerosol size distribution data

    Parameters
    ----------
    filename : string
        filename
    varname : string
        variable name

    Returns
    -------
    time : time of data
    size : center value of diameter bins
    data : data of the variable
    timeunit : unit of time
    dataunit : unit of data
    long_name : longname of data
    """
    f = Dataset(filename, 'r')
    # read in variables
    size = f.variables['size'][:]
    t_id = f.variables['time']
    time = t_id[:]
    timeunit = t_id.units
    d_id = f.variables[varname]
    data = d_id[:]
    dataunit = d_id.units
    long_name = d_id.long_name
    
    f.close()
    return(time, size, data, timeunit, dataunit, long_name)

#%% read exhaust-free ship data data
def read_ship_exhaustfree(filename, varname):
    """
    read variables from exhaust-free ship data

    Parameters
    ----------
    filename : string
        filename
    varname : string
        variable name

    Returns
    -------
    time : time of data
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