"""
functions of calculating time or changing time format
"""

#%%
def yyyymmdd2cday(datestr,calendar='leap'):
    """
    change from date string to calendar day

    Parameters
    ----------
    datestr : string
        yyyymmdd format or yyyy-mm-dd format
        e.g., '20190312', '2019-03-12'
    calendar : string, optional
        option if use leap year or nonleap year. The default is 'leap'.

    Returns
    -------
    cday : calendar day in this year

    """
    # must be yyyy-mm-dd format or yyyymmdd format
    if len(datestr)==10:
        dstr=datestr.split('-')
        yr=int(dstr[0])
        mon=int(dstr[1])
        day=int(dstr[2])
    elif len(datestr)==8:
        yr=int(datestr[0:4])
        mon=int(datestr[4:6])
        day=int(datestr[6:8])
    else:
        raise ValueError("input date format must be yyyymmdd or yyyy-mm-dd")
        
    
    if calendar=='leap' and yr%4 == 0:
        mmm=[31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        mmm=[31,28,31,30,31,30,31,31,30,31,30,31]
    
    
    cday=sum(mmm[0:mon-1])+day
    # if mon==1:
    #     cday=day
    # else:
    #     cday=sum(mmm[0:mon-1])+day
        
    return(cday)

#%%
def hhmmss2sec(timestr):
    """
    change time from hhmmss format to seconds of the day

    Parameters
    ----------
    timestr : string
        time, in hh:mm:ss format.
        e.g., 12:23:11

    Returns
    -------
    tsec : time in seconds of the day

    """
    #must be hh:mm:ss format
    tstr=timestr.split(':')
    hh=int(tstr[0])
    mm=int(tstr[1])
    ss=int(tstr[2])
    
    tsec=hh*3600+mm*60+ss
    
    return(tsec)


#%%
def cday2hhmmss(cday):
    """
    find the time from calendar day in hour-minute-second format

    Parameters
    ----------
    cday : float
        calendar day

    Returns
    -------
    hhmmss : string of time of the day

    """
    
    time=86400*(cday-int(cday))
    hh=int(time/3600)
    mm=int(time/60-hh*60)
    ss=int(time - hh*3600 - mm*60)
    
    hhmmss=str(hh).zfill(2) +':'+ str(mm).zfill(2) +':'+ str(ss).zfill(2)
    
    return(hhmmss)

#%% 
def cday2mmdd(cday,calendar='noleap'):
    """
    find the date from calendar day in month and day format

    Parameters
    ----------
    cday : float
        calendar day
    calendar : string, optional
        option if use leap year or noleap year. The default is 'noleap'.

    Returns
    -------
    mmdd : string of day of the year

    """
    
    if calendar=='leap':
        mmm=[31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        mmm=[31,28,31,30,31,30,31,31,30,31,30,31]
    for ii in range(13):
        if sum(mmm[0:ii])>=int(cday):
            mm = ii
            break
    dd = int(cday)-sum(mmm[0:mm-1])
    
    mmdd = str(mm).zfill(2) +str(dd).zfill(2)
    return(mmdd)

#%%
def datetime2cday(time):
    """
    change the time format from datetime64 to float of calendar day

    Parameters
    ----------
    time : numpy array of DateTime64 format
        input time array

    Returns
    -------
    cday0 : calendar day of yyyy-mm-dd hh:mm:ss in float format

    """
    
    # change time to calendar day
    sec_from_Jan1 = time.astype('datetime64[s]') - time.astype('datetime64[Y]')
    calday = sec_from_Jan1.astype('float64')/86400 + 1
    
    return(calday)

#%%

def timeunit2cday(timeunit,calendar='leap'):
    """
    find the calendar day of the initial time in timeunit

    Parameters
    ----------
    timeunit : string
        timeunit, in format of "*** since yyyy-mm-dd hh:mm:ss"
        *** may be days, hours, minutes or seconds
    calendar : string, optional
        option if use leap year or noleap year. The default is 'leap'.

    Returns
    -------
    cday0 : calendar day of yyyy-mm-dd hh:mm:ss in float format

    """
    tstr=timeunit.split(' ')
    cday0 = yyyymmdd2cday(tstr[2],calendar) + hhmmss2sec(tstr[3])/86400
    return(cday0)

# print(yyyymmdd2cday(datestr))
# print(hhmmss2sec(timestr))
