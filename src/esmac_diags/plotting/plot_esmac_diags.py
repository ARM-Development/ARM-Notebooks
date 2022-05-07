"""
# functions of plot diagnostics
"""
import matplotlib.pyplot as plt
import numpy as np

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def timeseries(time, data, figsize=(10,4), xlimit=None, ylimit=None, 
                 xlabel=None, ylabel=None, title=None, legend=None, 
                 color=['k','b','g','c','r','orange','gray'], **kwargs):
    """
    plot timeseries

    Parameters
    ----------
    time : list of 1-d xarrays
        input of time coordinate
    data : list of 1-d xarrays
        input timeseries
        other plotting parameters (marker, etc.) should also be list
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if legend is None:
        legend = [None for mm in range(ndata)]
        
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(figsize=figsize)
    for nn in range(ndata):
        ax.plot(time[nn],data[nn],color=color[nn],label=legend[nn], **kwargs)
    ax.legend()
    ax.set_ylim(ylimit)
    ax.set_xlim(xlimit)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title, fontsize=20)
    plt.tight_layout()

    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def timeseries_2d(time, height, data, figsize=None, levellist=None, title=None,
                 xlimit=None, ylimit=None, xticks=None, yticks=None, xlabel='Time', ylabel=None,
                 legend=None, **kwargs):
    """
    plot aerosol size timeseries

    Parameters
    ----------
    time : list of 1-d xarrays
        input of time coordinate
    height : list of 1-d xarrays
        input of vertical coordinate
    data : list of 2-d xarrays
        input timeseries of 2d data
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if figsize is None:
        figsize=(10, 2*(ndata))
    if legend is None:
        legend = [None for mm in range(ndata)]
    if levellist is None:
        levellist = np.linspace(np.nanmin(data[0]), np.nanmax(data[0]), 10)
    
    fig = plt.figure(figsize=figsize)
    plt.rcParams.update({'font.size': 16})
    ax = []
    for nn in range(ndata):
        ax0 = fig.add_subplot(ndata, 1, nn+1)
        ax.append(ax0)
        h = ax0.contourf(time[nn],height[nn],data[nn],levellist, **kwargs)
        ax0.set_ylim(ylimit)
        ax0.set_ylabel(ylabel)
        if xticks is not None:
            ax0.set_xticks(xticks)
        if yticks is not None:
            ax0.set_yticks(yticks)
        ax0.grid()
        ax0.text(0.013, 0.94, legend[nn], transform=ax0.transAxes,bbox=dict(facecolor='w',edgecolor='none'), verticalalignment='top')
        if nn<ndata-1:
            ax0.set_xticklabels([])
    ax[0].set_title(title, fontsize=18)
    ax[-1].set_xlabel(xlabel)
    # colorbar
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    cbar=fig.colorbar(h, cax=cax)
    cbar.ax.set_yticks(levellist)
    plt.subplots_adjust(right=0.9,bottom=0.1)
    
    return(fig, ax)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def timeseries_size(time, size, data, figsize=None, leveltick=[0.1,1,10,100,1000,10000, 100000],
                 xlimit=None, ylimit=(3,3000), xlabel='Time', ylabel=None, title=None,
                 legend=None, **kwargs):
    """
    plot aerosol size timeseries

    Parameters
    ----------
    time : list of 1-d xarrays
        input of time coordinate
    size : list of 1-d xarrays
        input of size coordinate
    data : list of 2-d xarrays
        input timeseries of aerosol size distribution
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if figsize is None:
        figsize=(10, 2*(ndata))
    if legend is None:
        legend = [None for mm in range(ndata)]
        
    levellist=np.arange(np.log10(leveltick[0]),np.log10(leveltick[-1])+0.1,.2)
    
    # for plotting in log scale, remove zeros by assigning small values to the lowest levelstick
    for nn in range(ndata):
        data[nn][data[nn]<leveltick[0]] = leveltick[0]
    
    fig = plt.figure(figsize=figsize)
    plt.rcParams.update({'font.size': 16})
    ax = []
    for nn in range(ndata):
        ax0 = fig.add_subplot(ndata, 1, nn+1)
        ax.append(ax0)
        h = ax0.contourf(time[nn],size[nn],np.log10(data[nn]),levellist,cmap=plt.get_cmap('jet'))
        ax0.set_yscale('log')
        ax0.set_ylim(ylimit)
        ax0.set_ylabel(ylabel)
        ax0.text(0.013, 0.94, legend[nn], transform=ax0.transAxes,bbox=dict(facecolor='w',edgecolor='none'), verticalalignment='top')
        if nn<ndata-1:
            ax0.set_xticklabels([])
    ax[0].set_title(title, fontsize=18)
    ax[-1].set_xlabel(xlabel)
    # colorbar
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    cbar=fig.colorbar(h, cax=cax, ticks=np.log10(leveltick))
    cbar.ax.set_yticklabels(leveltick)
    plt.subplots_adjust(right=0.9,bottom=0.1)
    
    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def mean_size(size, data, figsize=(8,6), xlimit=None, ylimit=None, xscale='log',yscale='log',
                 xlabel=None, ylabel=None, title=None, legend=None, linestyles=None,
                 marker=None, color=['k','b','g','c','r','orange','gray']):
    """
    plot timeseries

    Parameters
    ----------
    size : list of 1-d xarrays
        input of size coordinate
    data : list of 1-d xarrays
        input timeseries
        other plotting parameters (marker, etc.) should also be list
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if marker is None:
        marker = [None for mm in range(ndata)]
    if linestyles is None:
        linestyles = ['-' for mm in range(ndata)]
    if legend is None:
        legend = [None for mm in range(ndata)]
    
    if xlimit is None and xscale=='log':
        xlimit = (np.min([np.min(x) for x in size]), np.max([np.max(x) for x in size]))
    if ylimit is None and yscale=='log':
        ylimit = [np.nanmin([np.min(y) for y in data]+[1e-2]), np.nanmax([np.max(y) for y in data]+[1e4])]
        ylimit[0] = np.max([ylimit[0], 1e-4])
        
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(figsize=figsize)
    for nn in range(ndata):
        ax.plot(size[nn],data[nn],marker=marker[nn],linestyle=linestyles[nn],color=color[nn],label=legend[nn])
    ax.legend()
    ax.set_ylim(ylimit)
    ax.set_xlim(xlimit)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid()
    ax.set_title(title, fontsize=20)
    plt.tight_layout()

    return(fig, ax)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def diurnalcycle(data, nozero_percentile=False, figsize=(8,6), 
                 xlimit=(-0.5,24), xticks=np.arange(0,25,3), ylimit=None, 
                 xlabel='Hours', ylabel='Y-label', title=None, legend=None, 
                 fmt=None, marker=None, color=['k','b','g','c','r','orange','gray']):
    """
    plot diurnal cycle

    Parameters
    ----------
    data : list of 1-d xarrays
        input timeseries to plot diurnal cycle
        other plotting parameters (fmt, marker, etc.) should also be list
    nozero_percentile : True or False
        option of whether include zeros to calculate percentiles. 
        usually set as True for precipitation, LWP etc
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if marker is None:
        marker = ['.' for mm in range(ndata)]
    if fmt is None:
        fmt = ['-' for mm in range(ndata)]
    if legend is None:
        legend = [None for mm in range(ndata)]
        
    # errorbar with median/mean value
    # # observation diurnal cycle
    # data_dc = []
    # data_dc_error = []
    # data_dcmn = []
    # for nn in range(ndata):
    #     #interquartile range ("error"), median, and mean values
    #     # remove zero-precip samples
    #     if nozero_percentile==True:
    #         obs_nozero = data[nn]
    #         obs_nozero[obs_nozero==0] = np.nan
    #         data_dc1 = obs_nozero.groupby('time.hour').quantile([0.25, 0.5, 0.75])
    #     else:
    #         data_dc1 = data[nn].groupby('time.hour').quantile([0.25, 0.5, 0.75]) 
    #     data_dc_error1 = [data_dc1[:,1]-data_dc1[:,0], data_dc1[:,2]-data_dc1[:,1]]
    #     data_dcmn1 = data[nn].groupby('time.hour').mean()
    #     data_dc.append(data_dc1)
    #     data_dc_error.append(data_dc_error1)
    #     data_dcmn.append(data_dcmn1)
    # # make plot
    # # set position shift so that models and obs are not overlapped
    # p_shift = np.arange(ndata)
    # p_shift = (p_shift - p_shift.mean())*0.2
    # plt.rcParams.update({'font.size': 16})
    # fig,ax = plt.subplots(figsize=figsize)
    # for nn in range(ndata):
    #     ax.errorbar(np.arange(24)+p_shift[nn], data_dc[nn][:,1], yerr=data_dc_error[nn], fmt=fmt[nn],color=color[nn])
    #     ax.plot(np.arange(24)+p_shift[nn], data_dcmn[nn], linestyle='none',marker=marker[nn],color=color[nn],label=legend[nn])
            
    
    # percentile box plot
    # make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(ndata)
    p_shift = (p_shift - p_shift.mean())*0.2
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(figsize=figsize)
    for nn in range(ndata):
        c = color[nn]
        if nozero_percentile==True:
            data_nozero = data[nn]
            data_nozero[data_nozero==0] = np.nan
            dataa = data_nozero.groupby('time.hour')
        else:
            dataa = data[nn].groupby('time.hour')
        datab = [dataa[i].data for i in range(24)]
        datac = [d[~np.isnan(d)] for d in datab]
        ax.boxplot(datac,whis=(10,90),showmeans=False,showfliers=False,
                positions=np.arange(24)+p_shift[nn],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax.plot([],c=c, label=legend[nn])      
    ax.legend()
    
    ax.set_ylim(ylimit)
    ax.set_xlim(xlimit)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title, fontsize=20)
    plt.tight_layout()

    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def diurnalcycle_2d(data, y=None, figsize=None, x=np.arange(24), 
                 xlimit=(0,23), xticks=np.arange(0,24,3), ylimit=None, 
                 levellist=None,
                 xlabel='Hours', ylabel=None, title=None, **kwargs):
    """
    plot 2d (t-z) diurnal cycle

    Parameters
    ----------
    data : list of 2-d xarrays
        input timeseries to plot diurnal cycle
    y : list of numpy array, optional
        coordinate of the other dimension of the 2d data
    levellist: numpy array, optional
        determine the values of contour lines
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.contourf
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    # calculate mean diurnal cycle
    
    data_dc = []
    for nn in range(ndata):
        data_dc.append(data[nn].groupby('time.hour').mean())
        
    if y is None:
        y = [np.arange(data[mm].shape[0]) for mm in range(ndata)]
    if levellist is None:
        levellist = np.linspace(np.nanmin(data_dc[0]), np.nanmax(data_dc[0]), 10)
    if title is None:
        title = [None for mm in range(ndata)]
    if figsize is None:
        figsize=(6*(ndata),6)
    
    # make plot
    plt.rcParams.update({'font.size': 16})
    fig = plt.figure(figsize=figsize)
    for nn in range(ndata):
        ax1 = fig.add_subplot(1,ndata,nn+1)
        h0=ax1.contourf(x,y[nn],data_dc[nn],levellist, **kwargs)
        ax1.set_xlim(xlimit)
        ax1.set_ylim(ylimit)
        ax1.set_xticks(xticks)
        ax1.set_xlabel(xlabel)
        ax1.set_title(title[nn], fontsize=18)
        ax1.grid()
        if nn==0:
            ax=[ax1]
        else:
            ax.append(ax1)
    ax[0].set_ylabel(ylabel)
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    fig.colorbar(h0, cax=cax)    
    return(fig, ax)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def seasonalcycle(data,figsize=(9,6), xlimit=(0.5,12.5), xticks=np.arange(1,13,1), 
                 ylimit=None, nozero_percentile=False,
                 xlabel='Months', ylabel='Y-label', title=None, legend=None, 
                 fmt=None, marker=None, color=['k','b','g','c','r','orange','gray']):
    """
    plot seasonal cycle

    Parameters
    ----------
    data : list of 1-d xarrays
        input timeseries to plot seasonal cycle
        other plotting parameters (fmt, marker, etc.) should also be list
    nozero_percentile : True or False
        option of whether include zeros to calculate percentiles. 
        usually set as True for precipitation, LWP etc
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if marker is None:
        marker = ['.' for mm in range(ndata)]
    if fmt is None:
        fmt = ['-' for mm in range(ndata)]
    if legend is None:
        legend = [None for mm in range(ndata)]
    
    # # observation seasonal cycle
    # data_sc = []
    # data_sc_error = []
    # data_scmn = []
    # for nn in range(ndata):
    #     #interquartile range ("error"), median, and mean values
    #     # remove zero-precip samples
    #     if nozero_percentile==True:
    #         data_nozero = data[nn][data[nn]>0.00]
    #         data_sc1 = data_nozero.groupby('time.month').quantile([0.25, 0.5, 0.75])
    #     else:
    #         data_sc1 = data[nn].groupby('time.month').quantile([0.25, 0.5, 0.75]) 
    #     data_sc_error1 = [data_sc1[:,1]-data_sc1[:,0], data_sc1[:,2]-data_sc1[:,1]]
    #     data_scmn1 = data[nn].groupby('time.month').mean()
    #     data_sc.append(data_sc1)
    #     data_sc_error.append(data_sc_error1)
    #     data_scmn.append(data_scmn1)
    # # make plot
    # # set position shift so that models and obs are not overlapped
    # p_shift = np.arange(ndata)
    # p_shift = (p_shift - p_shift.mean())*0.2
    # plt.rcParams.update({'font.size': 16})
    # fig,ax = plt.subplots(figsize=figsize)
    # for nn in range(ndata):
    #     ax.errorbar(np.arange(1,13)+p_shift[nn], data_sc[nn][:,1], yerr=data_sc_error[nn], fmt=fmt[nn],color=color[nn])
    #     ax.plot(np.arange(1,13)+p_shift[nn], data_scmn[nn], linestyle='none',marker=marker[nn],color=color[nn],label=legend[nn])
    # ax.legend()
    
    # percentile box plot
    # make plot
    # set position shift so that models and obs are not overlapped
    p_shift = np.arange(ndata)
    p_shift = (p_shift - p_shift.mean())*0.2
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(figsize=figsize)
    for nn in range(ndata):
        c = color[nn]
        if nozero_percentile==True:
            data_nozero = data[nn]
            data_nozero[data_nozero==0] = np.nan
            dataa = data_nozero.groupby('time.month')
        else:
            dataa = data[nn].groupby('time.month')
        datab = [dataa[i].data for i in range(12)]
        datac = [d[~np.isnan(d)] for d in datab]
        ax.boxplot(datac,whis=(10,90),showmeans=False,showfliers=False,
                positions=np.arange(1,13)+p_shift[nn],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax.plot([],c=c, label=legend[nn])      
    ax.legend()
    
    ax.set_ylim(ylimit)
    ax.set_xlim(xlimit)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title, fontsize=20)
    plt.tight_layout()

    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def seasonalcycle_2d(data, figsize=None, x=np.arange(12)+1, y=None, 
                xlimit=(1,12), xticks=np.arange(1,13,1), 
                ylimit=None, levellist=None,
                xlabel='Months', ylabel=None, title=None, **kwargs):
    """
    plot 2d (t-z) seasonal cycle

    Parameters
    ----------
    data : list of 2-d xarrays
        input timeseries to plot seasonal cycle
    y : numpy array, optional
        coordinate of the other dimension of the 2d data
    levellist: numpy array, optional
        determine the values of contour lines
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.contourf
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    
    data_sc = []
    for nn in range(ndata):
        data_sc.append(data[nn].groupby('time.month').mean())
        
    if y is None:
        y=np.arange(data[0].shape[0])
    if levellist is None:
        levellist = np.linspace(np.nanmin(data_sc[0]), np.nanmax(data_sc[0]), 10)
    if title is None:
        title = [None for mm in range(ndata)]
    if figsize is None:
        figsize=(6*(ndata),6)
    
    # make plot
    plt.rcParams.update({'font.size': 16})
    fig = plt.figure(figsize=figsize)
    for nn in range(ndata):
        ax1 = fig.add_subplot(1,ndata,nn+1)
        h0=ax1.contourf(x,y,data_sc[nn],levellist, **kwargs)
        ax1.set_xlim(xlimit)
        ax1.set_ylim(ylimit)
        ax1.set_xticks(xticks)
        ax1.set_xlabel(xlabel)
        ax1.set_title(title[nn], fontsize=18)
        ax1.grid()
        if nn==0:
            ax=[ax1]
        else:
            ax.append(ax1)
    ax[0].set_ylabel(ylabel)
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    fig.colorbar(h0, cax=cax)    
    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def hist(data, figsize=(8, 5), yscale='linear',xlimit=None, ylimit=None, 
              xlabel=None, ylabel=None, title=None, legend=None, **kwargs):
    """
    plot histgram or 1d data

    Parameters
    ----------
    data : list of 1-d data
        input timeseries to plot histogram
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.hist
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    if ylabel is None:
        ylabel = 'Sample Number'

    # make plot
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(data,  **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yscale(yscale)
    ax.set_xlim(xlimit)
    ax.set_ylim(ylimit)
    ax.set_title(title)
    if legend is not None:
        ax.legend(legend)
    ax.grid()
    plt.tight_layout()
    return(fig, ax)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def jointhist(xdata, ydata, figsize=None, xedges=None, yedges=None, weight=None,
              normalize_x=False, xlimit=None, ylimit=None, 
              xlabel=None, ylabel=None, title=None, **kwargs):
    """
    plot joint histogram of two variables

    Parameters
    ----------
    xdata, ydata : list of array(s)
        input 1-d timeseries for joint histogram
    xedges, yedges : numpy array, optional
        edge value of data bins. 
    weight: list of array(s), optional
        weights of the data, same dimension of xdata and ydata
    normalize_x : True or False, optional
        options of whether histogram is normalized in x bins
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.pcolormesh
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(xdata)
    if xedges is None:
        xedges = np.arange(21)*np.nanmax(xdata[0])/20.
    if yedges is None:
        yedges = np.arange(21)*np.nanmax(ydata[0])/20.
    if title is None:
        title = [None for i in range(ndata)]
    if figsize is None:
        figsize=(6*(ndata),6)
    
    # calculate joint histogram    
    H_norm = []
    pdf_x = []
    for mm in range(ndata):
        H_count, x, y = np.histogram2d(xdata[mm], ydata[mm], bins=(xedges, yedges))
        # apply weight if available
        if weight is None:
            # normalize by sample number in each ccn bin
            if normalize_x == True:
                n_x, x2 = np.histogram(xdata[mm][~np.isnan(ydata[mm])], bins=xedges)
                n_x = np.float16(n_x)
                n_x[n_x==0] = np.nan
                H_count = H_count.T/n_x
                pdf_x.append(n_x)
            else:
                H_count = H_count.T/np.count_nonzero(~np.isnan(xdata[mm]*ydata[mm]))
        else:
            H_count = H_count.T * weight[mm]
        H_count[H_count == 0] = np.nan
        H_norm.append(H_count)

    # plot joint histogram
    X, Y = np.meshgrid(xedges, yedges)
    plt.rcParams.update({'font.size': 16})    
    # not normalize x
    if normalize_x == True:
        fig,ax = plt.subplots(2,ndata,figsize=figsize, 
                              sharex=True, gridspec_kw={'height_ratios':[5,1]})   # figsize in inches
        fig.subplots_adjust(hspace=0)
        h1=[]
        for mm in range(ndata):
            h = ax[0,mm].pcolormesh(X, Y, H_norm[mm], **kwargs)
            h1.append(h)
            ax[1,mm].plot(0.5*(X[0,0:-1]+X[0,1:]), pdf_x[mm],color='k')
            ax[1,mm].set_yticks([])
            if type(xlabel) is list:
                ax[1,mm].set_xlabel(xlabel[mm])
            else:
                ax[1,mm].set_xlabel(xlabel)
            ax[0,mm].set_title(title[mm], fontsize=18)
            ax[0,mm].grid()
            ax[1,mm].grid()
            if xlimit is not None:
                ax[0,mm].set_xlim(xlimit)
                ax[1,mm].set_xlim(xlimit)
            if ylimit is not None:
                ax[0,mm].set_ylim(ylimit)
            if type(ylabel) is list:
                ax[0,mm].set_ylabel(ylabel[mm])
            else:
                ax[0,0].set_ylabel(ylabel)
        ax[1,0].set_ylabel('PDF')
        ax[1,0].yaxis.set_label_coords(-0.13,0.6)
    else:
        fig = plt.figure(figsize=figsize)
        h1 = []
        ax = []
        for mm in range(ndata):
            ax1 = fig.add_subplot(1,ndata,mm+1)
            h = ax1.pcolormesh(X, Y, H_norm[mm], **kwargs)
            ax.append(ax1)
            h1.append(h)
            if xlimit is not None:
                ax1.set_xlim(xlimit)
            if ylimit is not None:
                ax1.set_ylim(ylimit)
            ax1.grid()
            ax1.set_title(title[mm], fontsize=18)
            if type(xlabel) is list:
                ax1.set_xlabel(xlabel[mm])
            else:
                ax1.set_xlabel(xlabel)
            if type(ylabel) is list:
                ax1.set_ylabel(ylabel[mm])
            else:
                ax[0].set_ylabel(ylabel)
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    fig.colorbar(h1[0], cax=cax) 
    
    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def heatmap(xdata, ydata, zdata, figsize=None, xedges=None, yedges=None,
              xlabel=None, ylabel=None, zlabel=None, title=None, **kwargs):
    """
    plot heatmaps (*median* value of zdata in each x-y bin) 
    

    Parameters
    ----------
    xdata, ydata, zdata : list of array(s)
        input 1-d timeseries for heatmaps
    xedges, yedges : numpy array, optional
        edge value of data bins. 
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.imshow
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(xdata)
    if xedges is None:
        xedges = np.arange(21)*np.nanmax(xdata[0])/20.
    if yedges is None:
        yedges = np.arange(21)*np.nanmax(ydata[0])/20.
    if title is None:
        title = [None for i in range(ndata)]
    if figsize is None:
        figsize=(6*(ndata),6)
    
    heatmaps = []
    samplenum = []
    for mm in range(ndata):
        heatmap_tmp = np.full((len(xedges)-1, len(yedges)-1), np.nan)
        sample_tmp = np.full((len(xedges)-1, len(yedges)-1), 0)
        for j in range(len(yedges)-1):
            for i in range(len(xedges)-1):
                mask = np.logical_and(np.logical_and(xdata[mm]>xedges[i], xdata[mm]<xedges[i+1]), \
                                      np.logical_and(ydata[mm]>yedges[j], ydata[mm]<yedges[j+1])).data
                heatmap_tmp[i, j] = np.nanmedian(zdata[mm][mask])
                # heatmap_tmp[i, j] = np.nanmean(zdata[mm][mask])
                sample_tmp[i, j] = np.size(np.where(mask == True))
        heatmaps.append(heatmap_tmp)
        samplenum.append(sample_tmp)
        
    xticklabels = (xedges[0:-1] + xedges[1:])/2
    yticklabels = (yedges[0:-1] + yedges[1:])/2
    xticks = np.arange(0, xticklabels.size, 3)
    yticks = np.arange(0, yticklabels.size, 2)
    plt.rcParams.update({'font.size': 16})
    fig = plt.figure(figsize=figsize)
    h = []
    ax = []
    for mm in range(ndata):
        ax1 = fig.add_subplot(1,ndata,mm+1)
        ax.append(ax1)
        h0 = ax[mm].imshow(heatmaps[mm].T, origin='lower', **kwargs)
        h.append(h0)
        ax[mm].grid()
        ax[mm].set_xticks(xticks)
        ax[mm].set_xticklabels(xticklabels[xticks].astype(np.int))
        ax[mm].set_yticks(yticks)
        ax[mm].set_yticklabels(yticklabels[yticks].astype(np.int))
        ax[mm].set_title(title[mm], fontsize=18)
        if type(xlabel) is list:
            ax[mm].set_xlabel(xlabel[mm])
        else:
            ax[mm].set_xlabel(xlabel)
        if type(ylabel) is list:
            ax[mm].set_ylabel(ylabel[mm])
        else:
            ax[0].set_ylabel(ylabel)
    cax = plt.axes([0.92, 0.2, 0.02, 0.6])
    fig.colorbar(h[0], cax=cax) 
    plt.text(1, 0, zlabel)
    
    return(fig, ax)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def scatter(xdata, ydata, figsize=None, xlimit=None, ylimit=None, 
                 xlabel='xdata', ylabel='ydata', title=None, 
                 linear_fit=False, intercept=False, **kwargs):
    """
    plot scatter plot of two variables

    Parameters
    ----------
    xdata, ydata : list of arrays
        input 1-d timeseries for scatter plot
    linear_fit : True of False, optional
        options of whether calculating and plotting linear fit
    intercept : True or False, optional
        if linear_fit is True, options of whether the linear fit has intercept
        intercept == True : y=ax+b
        intercept == False : y=ax
    
    Other Parameters (**kwargs))
    ----------
    refer to matplotlib.pyplot.scatter
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
        
    ndata = len(xdata)
    if title is None:
        title = [None for i in range(ndata)]
    if figsize is None:
        figsize=(6*(ndata),6)
    if xlimit is None:
        xlimit = (min([np.nanmin(dd) for dd in xdata]), max([np.nanmax(dd) for dd in xdata]))
    if ylimit is None:
        ylimit = (min([np.nanmin(dd) for dd in ydata]), max([np.nanmax(dd) for dd in ydata]))
    
    plt.rcParams.update({'font.size': 16})
    fig = plt.figure(figsize=figsize)
    ax = []
    for mm in range(ndata):
        ax1 = fig.add_subplot(1,ndata,mm+1)
        ax.append(ax1)
        ax1.scatter(xdata[mm],ydata[mm], **kwargs)
        ax1.set_xlabel(xlabel)
        ax1.set_xlim(xlimit)
        ax1.set_ylim(ylimit)
        ax1.set_title(title[mm])
        # plot 1:1 line
        minval = min(xlimit[0], ylimit[0])
        maxval = max(xlimit[1], ylimit[1])
        ax1.plot([minval, maxval],[minval, maxval],color='silver',linestyle=':')
        ax1.grid()
        
        if linear_fit==True:
            idx = np.logical_and(~np.isnan(xdata[mm]), ~np.isnan(ydata[mm]))
            if intercept==True:   # y=ax+b
                coef = np.polyfit(xdata[mm][idx],ydata[mm][idx],1)
                poly1d_fn = np.poly1d(coef)
                ax1.plot([minval, maxval],poly1d_fn([minval, maxval]),color='k')
                if coef[1]>=0:
                    fit_line = 'y = '+format(coef[0],'2.2f')+'x + '+format(coef[1],'2.2f')
                else:
                    fit_line = 'y = '+format(coef[0],'2.2f')+'x - '+format(-coef[1],'2.2f')
                print('linear fit: ' + fit_line)
                ax1.text(0.45*xlimit[1], 0.9*ylimit[1], fit_line)
            else:   # y=ax
                x = xdata[mm][idx,np.newaxis]
                a,_,_,_ = np.linalg.lstsq(x, ydata[mm][idx])
                # # correlation coefficient and p value of linear fit
                # import scipy.stats
                # r, p = scipy.stats.pearsonr(xdata[idx],ydata[idx])
                y=a*x
                ax1.plot(x,y,color='k')
                fit_line = 'y = '+format(a[0],'2.2f')+'x'
                print('linear fit: ' + fit_line)
                # ax1.text(x[-1], y[-1], fit_line)
                ax1.text(0.7*xlimit[1], 0.9*ylimit[1], fit_line)
        else:
            fit_line = 'None'
    
    # plt.tight_layout()
    ax[0].set_ylabel(ylabel)
    return(fig, ax)
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def bar(dataall, figsize=None, datalabel=None, varlabel=None, xlabel=None, ylabel=None, 
                 ylimit=None,title=None, 
                 colorall=['k','b','g','c','r','limegreen','lightblue','orange','silver']):
    """
    plot bar plot for aerosol compositions

    Parameters
    ----------
    dataall : list of 1-d arrays
        input data for all data source [datagroup1, datagroup2, datagroup3, ...]
        each data group includes several variables, missing variable are placehold with []:
            datagroup1 = [var1, var2, var3, ...]
            datagroup2 = [var1, [], var3, ...]
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(dataall)
    nvars = len(dataall[0])
    if datalabel is None:
        datalabel = ['data'+format(nn) for nn in range(ndata)]
    if varlabel is None:
        varlabel = ['var'+format(nn) for nn in range(nvars)]
    if figsize is None:
        figsize = (ndata*4,4)
            
    plt.rcParams.update({'font.size': 16})
    fig,ax = plt.subplots(1,1,figsize=figsize)   # figsize in inches
    
    for nn in range(ndata):
        if len(dataall[nn]) != nvars:
            raise ValueError('each data group should have same number of elements with label')
        barsize = np.array([np.nanmean(vv) for vv in dataall[nn]])
        barsize[np.isnan(barsize)] = 0
        for ii in range(nvars):
            if nn==0:
                ax.bar(nn,barsize[ii],bottom=sum(barsize[0:ii]),width=0.6,color=colorall[ii],label=varlabel[ii])
            else:
                ax.bar(nn,barsize[ii],bottom=sum(barsize[0:ii]),width=0.6,color=colorall[ii])
        
    ax.set_xticks(np.arange(ndata))
    ax.set_xticklabels(datalabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='right', shadow=False, bbox_to_anchor=(1.02+0.6/ndata, 0.5))

    return(fig, ax)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def percentile_lat(data, lat, latbin, figsize=(8,2), xlimit=None, ylimit=None, 
                 xlabel='Latitude', ylabel=None, title=None, legend=None, 
                 color=['k','b','g','c','r']):
    """
    plot vertical percentiles for aircraft measurements

    Parameters
    ----------
    data : list of 1-d arrays
        input timeseries in different latitude
    lat : list of 1-d arrays
        input timeseries of latitude, corresponding to the data
    latbin : 1-d array
        center of lat bins to calculate percentiles
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata=len(data)
    if legend is None:
        legend = [None for mm in range(ndata)]
    
    # put data in different height bins
    dlat = latbin[1]-latbin[0]
    latmin = latbin-dlat/2
    latmax = latbin+dlat/2
    latlen = len(latbin)
    data_lat = []
    for mm in range(ndata):
        data_lat.append([])
        data2 = data[mm].data
        idx0 = ~np.isnan(data2)
        data2 = data2[idx0]
        lat2 = lat[mm][idx0]
        for zz in range(latlen):
            idx = np.logical_and(lat2>=latmin[zz], lat2<latmax[zz])
            data_lat[mm].append(data2[idx])
    
    # make plot
    # set position shift so that they are not overlapped
    p_shift = np.arange(ndata)
    p_shift = (p_shift - p_shift.mean())*0.2
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(figsize=figsize)
    for mm in range(ndata):
        c = color[mm]
        ax.boxplot(data_lat[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(latlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=True, patch_artist=True)    # need patch_artist to fill color in box
        ax.plot([],c=c, label=legend[mm])
    ax.legend(loc='upper right', fontsize='medium')
    ax.set_xlim(-1,latlen)
    ax.set_xticks(np.arange(-0.5*dlat,latlen-1,2))
    ax.set_xticklabels([int(np.floor(a)) for a in latbin[0::2]])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlimit)
    ax.set_ylim(ylimit)
    ax.set_title(title)

    return(fig, ax)    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def percentile_z(data, height, height_bin, figsize=(3,8), xlimit=None, ylimit=None, 
                 xlabel=None, ylabel=None, title=None, legend=None, 
                 color=['k','b','g','c','r']):
    """
    plot vertical percentiles for aircraft measurements

    Parameters
    ----------
    data : list of 1-d arrays
        input timeseries in different height
    height : list of 1-d arrays
        input timeseries of height, corresponding to the data
    height_bin : 1-d array
        center of height bins to calculate percentiles
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata=len(data)
    if legend is None:
        legend = [None for mm in range(ndata)]
    
    # put data in different height bins
    z = height_bin
    dz = z[1]-z[0]
    zmin=z-np.insert((z[1:]-z[0:-1])/2,0,dz)
    zmax=z+np.append((z[1:]-z[0:-1])/2,dz)
    zlen=len(z)
    data_z = []
    for mm in range(ndata):
        data_z.append([])
        data2 = data[mm]
        idx0 = ~np.isnan(data2)
        data2 = data2[idx0]
        h2 = height[mm][idx0]
        for zz in range(zlen):
            idx = np.logical_and(h2>=zmin[zz], h2<zmax[zz])
            data_z[mm].append(data2[idx])
    
    # make plot
    # set position shift so that they are not overlapped
    p_shift = np.arange(ndata)
    p_shift = (p_shift - p_shift.mean())*0.2
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(figsize=figsize)
    for mm in range(ndata):
        c = color[mm]
        ax.boxplot(data_z[mm],whis=(5,95),showmeans=False,showfliers=False,
                positions=np.array(range(zlen))+p_shift[mm],widths=0.15,
                boxprops=dict(facecolor=c, color=c),whiskerprops=dict(color=c),
                medianprops=dict(color='lightyellow',linewidth=1),capprops=dict(color=c),
                vert=False, patch_artist=True)    # need patch_artist to fill color in box
        ax.plot([],c=c, label=legend[mm])
    ax.legend(loc='upper right', fontsize='medium')
    ax.set_ylim(-1,zlen)
    ax.set_yticks(range(zlen))
    ax.set_yticklabels(z)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlimit)
    ax.set_ylim(ylimit)
    ax.set_title(title)

    return(fig, ax)    


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def vertical(data, z, figsize=(4,6), xlimit=None, ylimit=None, errorevery=1,
                 xlabel=None, ylabel=None, title=None, legend=None, pct=[0.25,0.75],
                 fmt=None, marker=None, color=['k','b','g','c','r','orange','gray']):
    """
    plot vertical profiles and percentiles

    Parameters
    ----------
    data : list of arrays
        input 2-d timeseries
    z : 1-d array
        coordinate of vertical dimension
    errorevery : int, optional
        options of plotting errorbar every "errorevery" data points
    pct : [float, float]
        lower and upper percentile levels
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    ndata = len(data)
    if marker is None:
        marker = ['.' for mm in range(ndata)]
    if fmt is None:
        fmt = ['-' for mm in range(ndata)]
    if legend is None:
        legend = [None for mm in range(ndata)]
        
        
    if data[0].shape[1]!=z.shape[0]:
        raise ValueError('the second dimension of data should be consistent with z dimension')
    ndata = len(data)
    data_pct = []
    data_pct_error = []
    for mm in range(ndata):
        data_pct_tmp = np.zeros((np.size(z),3))
        for ii in np.arange(np.size(z)):
            data_pct_tmp[ii,:] = data[mm][:,ii].quantile([pct[0], 0.5, pct[1]])
        data_pct.append(data_pct_tmp)
        data_pct_error.append([data_pct_tmp[:,1]-data_pct_tmp[:,0], data_pct_tmp[:,2]-data_pct_tmp[:,1]])
    
    plt.rcParams.update({'font.size': 16})
    fig, ax = plt.subplots(figsize=figsize)
    for mm in range(ndata):
        ax.errorbar(data_pct[mm][:,1], z, xerr=data_pct_error[mm], errorevery=errorevery,\
                      fmt=fmt[mm], marker=marker[mm],color=color[mm],label=legend[mm])
    ax.legend()
    ax.set_xlim(xlimit)
    ax.set_ylim(ylimit)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=18)
    ax.grid()
    plt.tight_layout()  
    return(fig, ax)