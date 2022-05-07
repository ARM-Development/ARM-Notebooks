'''
functions of change time resolution (e.g., average into a coarser resolution).
'''

import numpy as np
#%%
def avg_time_1d(time0, data0, time):
    """
    average 1d data into coarser time resolution

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if data0.shape[0] != len(time0):
        raise ValueError("Arrays must have the same size")
    data = np.full((len(time)), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt] = np.nanmean(data0[idx], axis = 0)
    return(data)

#%%
def avg_time_2d(time0, data0, time):
    """
    average 2d data into coarser time resolution

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if data0.shape[0] != len(time0):
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time), data0.shape[1]), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt, :] = np.nanmean(data0[idx, :], axis = 0)
    return(data)

#%%
def weightmean_time_1d(time0,data0,weight,time):
    """
    average 1d data into coarser time resolution with weight

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    weight : numpy array
        weighting the input data for averaging
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if data0.shape[0]!=len(time0) or data0.shape!=weight.shape:
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time)),np.nan)
    dt=(time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0>=time[tt]-dt,time0<=time[tt]+dt)
        wtsum = np.nansum(weight[idx])
        dataall = data0[idx]*weight[idx]/wtsum
        if ~np.isnan(dataall).all():  # skip hours with all missing
            data[tt]=np.nansum(dataall) 
    return(data)

#%%
def median_time_1d(time0, data0, time):
    """
    rescale 1d data into coarser time resolution
    get median value in each coarser time window

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if len(data0) != len(time0):
        raise ValueError("Arrays must have the same size")
    data = np.full((len(time)), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt] = np.nanmedian(data0[idx], axis = 0)
    return(data)

#%%
def median_time_2d(time0, data0, time):
    """
    rescale 2d data into coarser time resolution
    get median value in each coarser time window

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if data0.shape[0] != len(time0):
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time), data0.shape[1]), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt,:] = np.nanmedian(data0[idx,:], axis = 0)
        # for dd in range(data0.shape[1]):
        #     if sum(~np.isnan(data0[idx,dd]))/sum(idx) > 0.4:
        #         data[tt,dd] = np.nanmedian(data0[idx,dd], axis = 0)
    return(data)

#%%
def median_time_forflight_1d(time0, data0, time, height, hdiff=50.):
    """
    Note: this is for flight measurements only !!!
    rescale 1d data into coarser time resolution, with a selection on flight height variance
    in each coarser time window, if height changes too much (aircraft is ascending or decending),
    either exclude the entire window (if flight height varies too much for the entire window)
    or only count the time when flight height does not vary (if flight height are mostly stable)

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data
    height: numpy array
        flight height
    hdiff: float
        height difference threshold (unit: m) to determine if data is chosen:
            if 25% and 75% percentile heights differ > hdiff, exclude this time window
            else, get median value of data with height between median +/- hdff of the time window

    Returns
    -------
    data : output data

    """
    if len(data0) != len(time0) or len(height) != len(time0):
        raise ValueError("Input arrays must have the same size")
    data = np.full((len(time)), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        if any(idx)==False:
            data[tt] = np.nan
            continue
        h = height[idx]
        h1 = np.nanpercentile(h, 25)
        h2 = np.nanpercentile(h, 75)
        if (h2-h1) > hdiff:
            data[tt] = np.nan
        else:
            data1 = data0[idx]
            h0 = np.nanmedian(h)
            # keep data in flight heights between median +/- hdiff
            idx2 = np.logical_and(h>(h0-hdiff), h<(h0+hdiff))
            data[tt] = np.nanmedian(data1[idx2], axis = 0)
    return(data)

#%%
def median_time_forflight_2d(time0, data0, time, height, hdiff=50.):
    """
    Note: this is for flight measurements only !!!
    rescale 2d data (time, dim2) into coarser time resolution, with a selection on flight height variance
    in each coarser time window, if height changes too much (aircraft is ascending or decending),
    either exclude the entire window (if flight height varies too much for the entire window)
    or only count the time when flight height does not vary (if flight height are mostly stable)

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data
    height: numpy array
        flight height
    hdiff: float
        height difference threshold (unit: m) to determine if data is chosen:
            if 25% and 75% percentile heights differ > hdiff, exclude this time window
            else, get median value of data with height between median +/- hdff of the time window

    Returns
    -------
    data : output data

    """
    if len(height) != len(time0):
        raise ValueError("Input arrays must have the same size")
    if data0.shape[0] != len(time0):
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time), data0.shape[1]), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        if any(idx)==False:
            data[tt, :] = np.nan
            continue
        h = height[idx]
        h1 = np.nanpercentile(h, 25)
        h2 = np.nanpercentile(h, 75)
        if (h2-h1) > hdiff:
            data[tt, :] = np.nan
        else:
            data1 = data0[idx, :]
            h0 = np.nanmedian(h)
            # keep data in flight heights between median +/- hdiff
            idx2 = np.logical_and(h>(h0-hdiff), h<(h0+hdiff))
            data[tt, :] = np.nanmedian(data1[idx2, :], axis = 0)
    return(data)
#%%
def median_time_2d(time0, data0, time):
    """
    rescale 2d data (time, dim2) into coarser time resolution
    get median value in each coarser time window

    Parameters
    ----------
    time0 : numpy array
        time dimension for input data
    data0 : numpy array
        input data
    time : numpy array
        time dimension for output data

    Returns
    -------
    data : output data

    """
    if data0.shape[0] != len(time0):
        raise ValueError("the first dimension of input data must have the same size with time")
    data = np.full((len(time), data0.shape[1]), np.nan)
    dt = (time[1]-time[0])/2
    for tt in range(len(time)):
        idx = np.logical_and(time0 >= time[tt]-dt, time0 <= time[tt] + dt)
        data[tt, :] = np.nanmedian(data0[idx, :], axis = 0)
    return(data)

