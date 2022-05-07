"""
# functions of calculate statistics of data
"""
import numpy as np
import scipy.stats

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def mean_std_percentiles(data,legend=None,outfile=None):
    """
    calculate mean value, standard deviation, [5,25,50,75,95] percentiles

    Parameters
    ----------
    data : list of 1-d xarrays
        input timeseries
    legend : list of str
        legend of data
    outfile : str
        filename of output statistics
    
    Returns
    -------
    mean
    std
    percentiles

    """
    ndata = len(data)
    if legend is None:
        legend=['' for i in range(ndata)]
    
    # calculate statistics
    meanall=[]
    stdall=[]
    pctall=[]
    for i in range(ndata):
        mean = np.nanmean(data[i])
        std = np.nanstd(data[i])
        pct = [5,25,50,75,95]
        percentiles = np.nanpercentile(data[i],pct)
        
        meanall.append(mean)
        stdall.append(std)
        pctall.append(percentiles)
        
        print('-------- statistics '+legend[i]+' ---------')
        print('mean    std.dev.    [5/25/50/75/95] percentiles')
        print(mean,std,percentiles)
    
    # output as txt file
    # print('write statistics to file '+outfile)
    if outfile is not None:
        with open(outfile, 'w') as f:
            for i in range(ndata):
                f.write('\n ************************************************ \n')
                f.write('mean statistics for '+legend[i]+'. sample size '+
                        format(np.int64(sum(~np.isnan(data[i]))))+'\n')
                # write mean
                f.write('\n mean: \t')
                f.write(format(meanall[i],'10.2f'))
                # write std
                f.write('\n std. dev.: ')
                f.write(format(stdall[i],'10.2f'))
                # write percentiles
                f.write('\n\n percentile: ')
                for j in range(len(pct)):
                    f.write('\n '+format(pct[j])+'%: ')
                    f.write(format(pctall[i][j],'10.2f'))
                f.write('\n')
    
    return(meanall,stdall,pctall)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def bias_corrcoef_RMSE(data1,data2,outfile,label1='data',label2='reference'):
    """
    calculate mean bias, correlation coefficient (include p-value) and RMSE of data1 related to data2

    Parameters
    ----------
    data1 : list of 1-d xarrays
        input timeseries
    data2 : list of 1-d xarrays
        input timeseries of reference data
    outfile : str
        filename of output statistics
    label1 : str
        label/description of data1 (analysis data)
    label2 : str
        label/description of data2 (reference data)
    
    Returns
    -------
    bias
    corrcoef
    RMSE

    """
    idx = ~np.logical_or(np.isnan(data1), np.isnan(data2))
    dataa=data1[idx]
    datab=data2[idx]
    
    bias = np.nanmean(dataa)-np.nanmean(datab)
    corr = scipy.stats.pearsonr(dataa,datab)
    corrcoef = [corr[0],corr[1]]
    rmse = np.sqrt(np.nanmean((dataa-datab)**2))
    
    print('bias    [corr, pvalue]    RMSE')
    print(bias,corrcoef,rmse)
    
    # print('write statistics to file '+outfile)
    with open(outfile, 'w') as f:
        f.write('mean statistics of '+label1+' related to '+label2+'. sample size '+format(np.int64(sum(idx)))+'\n')
        # write mean bias
        f.write('\n mean bias: ')
        f.write(format(bias,'10.2f'))
        # write corr coef
        f.write('\n corrcoef: ')
        f.write(format(corrcoef[0],'10.2f'))
        f.write('\n P-value: ')
        f.write(format(corrcoef[1],'10.4f'))
        # write percentiles
        f.write('\n RMSE: ')
        f.write(format(rmse,'10.2f'))
            
    return(bias,corrcoef,rmse)
