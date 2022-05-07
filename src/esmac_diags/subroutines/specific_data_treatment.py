"""
function of some specific data treatment
"""

import numpy as np

#%% 
def calc_cdnc_ARM(lwp,cod,H):
    """
    calculate cloud droplet number concentration using ARM surface-based retrievals
    references: https://www.arm.gov/publications/tech_reports/doe-sc-arm-tr-140.pdf 

    Parameters
    ----------
    lwp : numpy array
        liquid water path, unit: mm (kg/m2)
    cod : numpy array
        cloud optical depth, unit: N/A
    H : numpy array
        cloud depth, unit: m

    Returns
    -------
    Nd : output data
        column-integrated layer-mean cloud droplet number concentration. 
        unit: cm-3

    """
    # constants and parameters
    C1 = 0.05789
    k = 0.74
    rho_liq = 1000.   # unit: kg/m3
    
    Nd = 1e-6*(C1/k)*(2.**0.5)*(rho_liq**2)*(cod**3)/(lwp**2)/H 
    Nd[Nd == 0] = np.nan
    return(Nd)

#%% 
def calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=0.8):
    """
    calculate cloud droplet number concentration using VISST satellite retrievals
    options for different adiabaticity
    default value is 80% (0.8) adiabaticity following the reference
    references: follows Bennartz (2007, JGR)

    Parameters
    ----------
    lwp : numpy array
        liquid water path, unit: mm (kg/m2)
    ctt : numpy array
        average cloud temperature, for liquid only, unit: K
    cod : numpy array
        cloud optical depth, unit: N/A

    Returns
    -------
    Nd : output data
        column-integrated layer-mean cloud droplet number concentration. 
        unit: cm-3

    """
    # constants and parameters
    G = 9.8
    Cp = 1005.7
    Rd = 287.
    Rv = 461.
    lv = 2.477e6 #at 10 C
    epsilon = Rd/Rv
    pres_const = 85000. #Pa (used by Bennartz, 2007, JGR), could use cloud top pressure, but shouldn't alter estimates much
    Q = 2.
    k = 0.74 #Bennartz uses 0.8 +/- 0.1 but can be 0.5-0.9
    rho_liq = 1000.
        
    rho_air = pres_const/(Rd*ctt)
    es = 611.2*np.exp(17.62*(ctt-273.15)/(243.12 + ctt - 273.15))
    ws = epsilon*es/(pres_const - es)
    gamma_w = G*((1 + lv*ws/(Rd*ctt))/(Cp + lv**2*ws*epsilon/(Rd*ctt**2)))
    gamma_ad = (((epsilon + ws)*ws*lv*gamma_w)/(Rd*ctt**2) - (G*ws*pres_const/(Rd*ctt*(pres_const - es))))*rho_air
    Nd = 1e-6*(cod**3/k)*((2*(1e-3*lwp))**(-2.5))*((0.6*np.pi*Q)**(-3))*((3./(4.*np.pi*rho_liq))**(-2))*((adiabaticity*gamma_ad)**0.5) 
    return(Nd)

#%% 
def calc_clouddepth_VISST(lwp, ctt, adiabaticity=0.8):
    """
    calculate cloud depth using VISST satellite retrievals
    options for different adiabaticity
    default value is 80% (0.8) adiabaticity following the reference
    references: follows Bennartz (2007, JGR)

    Parameters
    ----------
    lwp : numpy array
        liquid water path, unit: mm (kg/m2)
    ctt : numpy array
        average cloud temperature, for liquid only, unit: K

    Returns
    -------
    Nd : output data
        column-integrated layer-mean cloud droplet number concentration. 
        unit: cm-3

    """
    # constants and parameters
    G = 9.8
    Cp = 1005.7
    Rd = 287.
    Rv = 461.
    lv = 2.477e6 #at 10 C
    epsilon = Rd/Rv
    pres_const = 85000. #Pa (used by Bennartz, 2007, JGR), could use cloud top pressure, but shouldn't alter estimates much
        
    rho_air = pres_const/(Rd*ctt)
    es = 611.2*np.exp(17.62*(ctt-273.15)/(243.12 + ctt - 273.15))
    ws = epsilon*es/(pres_const - es)
    gamma_w = G*((1 + lv*ws/(Rd*ctt))/(Cp + lv**2*ws*epsilon/(Rd*ctt**2)))
    gamma_ad = (((epsilon + ws)*ws*lv*gamma_w)/(Rd*ctt**2) - (G*ws*pres_const/(Rd*ctt*(pres_const - es))))*rho_air
    H = (2.*1e-3*lwp/(adiabaticity*gamma_ad))**0.5
    return(H)

#%% 
def find_nearest(xall, yall, x, y):
    """
    find the index of nearest point at 2-d space
    for E3SM spectrum element core output, so xall and yall are longitude and latitude
    but in 1-dimensional

    Parameters
    ----------
    xall : 1-d numpy array
        x dimension value of the 2d space
    yall : 1-d numpy array
        y dimension value of the 2d space.
    x : float or int
        x value of the given point.
    y : float or int
        y value of the given point.

    Returns
    -------
    None.

    """
    distance = np.square(xall-x) + np.square(yall-y)
    idx = distance.argmin()
    return(idx)

#%%
def insolation(time, lon, lat, leap_year='noleap'):
    """
    calculate insolation from given time and location

    Parameters
    ----------
    time : 1-d numpy array
        time in calendar day
    lon : float
        longitude
    lat : float
        latitude
    leap_year : str
        leap year has 366 days for a year

    Returns
    -------
    ins : 1-d numpy array
        insolation (W/m2) for the given time and lat/lon

    """
    
    days_in_year = 365.0
    # check if this is for the leap year
    if leap_year=='leap':
        days_in_year = 366.0
    
    thepi = 3.14159265
    tw=2.0*thepi*(time-1.)/days_in_year
    
    # eccentricity
    ecc= 1.000110+0.034221*np.cos(tw)+0.001280*np.sin(tw)+ 0.000719*np.cos(2.0*tw)+0.000077*np.sin(2.0*tw)
    
    pif=thepi/180.0
    
    delta=0.006918 - 0.399912*np.cos(tw) + 0.070257*np.sin(tw) - 0.006758*np.cos(2.*tw) + \
        0.000907*np.sin(2.*tw) - 0.002697*np.cos(3.*tw) + 0.001480*np.sin(3.*tw)
    
    cr1=pif*279.367 + 0.985647*days_in_year/360.0*tw
    
    dt=-105.4*np.sin(cr1)+596.2*np.sin(2*cr1)+4.3*np.sin(3*cr1)-12.7*np.sin(4*cr1)- \
        429.2*np.cos(cr1)-2.1*np.cos(2*cr1)+19.3*np.cos(3*cr1) # in seconds
    
    tt=dt/86400.0*2.0*thepi/days_in_year + tw
    
    cosz = np.sin(lat*pif)*np.sin(delta) - np.cos(lat*pif)*np.cos(delta)*np.cos(tt*days_in_year+lon*pif)
    
    ins=1368.2*cosz*ecc
    #ins=1378.95*cosz*ecc
    
    ins[ins<0]=0.0
    
    return(ins)

#%%
def lwc2cflag(lwc, lwcunit):
    """
    estimate cloud flag based on LWC

    Parameters
    ----------
    lwc : numpy array
        liquid water content data
    lwcunit : string
        unit of lwc

    Returns
    -------
    cldflag : estimated cloud flag

    """
    if lwcunit == 'kg/m3':
        lwc = lwc*0.001
        lwcunit = 'g/m3'
    elif lwcunit == 'g m-3':
        lwcunit = 'g/m3'
    if lwcunit not in ['g/m3', 'gram/m3']:
        print('unit of LWC: ' + lwcunit)
        raise ValueError("unit of LWC should be gram/m3 or kg/m3")

    cldflag = 0*np.array(lwc)
    
    # set threshold of LWC to identify cloud
    cldflag[lwc > 0.02] = 1
    return(cldflag)

#%%
def mask_model_ps(timem, psm, legnum, campaign, shipmetpath):
    """
    set model masks if the difference of Ps with observation is too large

    Parameters
    ----------
    timem : numpy array
        time in model
    psm : numpy array
        surface pressure in model
    legnum : string
        leg number,  or trip number
    campaign : string
        name of field campaign
    shipmetpath : string
        file path of shipmet data

    Returns
    -------
    datamask : mask flag of large Ps difference

    """
    import glob
    from ..subroutines.read_ship import read_marmet
    from ..subroutines.read_ARMdata import read_met
    from ..subroutines.time_format_change import yyyymmdd2cday,  cday2mmdd
    
    if campaign == 'MAGIC':
        filenameo = shipmetpath + 'marmet' + legnum + '.txt'
        (shipdata, shipvarlist) = read_marmet(filenameo)
        year = [a[1] for a in shipdata]
        month = [a[2] for a in shipdata]
        day = [a[3] for a in shipdata]
        hh = [int(a[4]) for a in shipdata]
        mm = [int(a[5]) for a in shipdata]
        ss = [int(a[6]) for a in shipdata]
        yyyymmdd = [year[i] + month[i] + day[i] for i in range(len(year))]   # yyyymmdd
        # get time in calendar day
        time = np.array(hh)/24. + np.array(mm)/1440. + np.array(ss)/86400. 
        time = np.array([time[i] + yyyymmdd2cday(yyyymmdd[i], 'noleap') for i in range(len(time))])
        if time[-1] < time[0]:
            time[time <= time[-1]] = time[time <= time[-1]] + 365
        # get variables
        ps = np.array([float(a[shipvarlist.index('bp')]) for a in shipdata])    
        ps[ps == -999] = np.nan

    elif campaign == 'MARCUS':
        if legnum == '1':
            startdate = '2017-10-30'
            enddate = '2017-12-02'
        elif legnum == '2':
            startdate = '2017-12-13'
            enddate = '2018-01-11'
        elif legnum == '3':
            startdate = '2018-01-16'
            enddate = '2018-03-04'
        elif legnum == '4':
            startdate = '2018-03-09'
            enddate = '2018-03-22'
        cday1 = yyyymmdd2cday(startdate, 'noleap')
        cday2 = yyyymmdd2cday(enddate,  'noleap')
        if startdate[0:4] != enddate[0:4]:
            cday2 = cday2 + 365  # cover two years
        time = np.empty(0)
        ps = np.empty(0)
        for cc in range(cday1, cday2 + 1):
            if cc <= 365:
                yyyymmdd = startdate[0:4] + cday2mmdd(cc)
            else:
                yyyymmdd = enddate[0:4] + cday2mmdd(cc-365)
            lst0 = glob.glob(shipmetpath + 'maraadmetX1.b1.' + yyyymmdd + '*')
            (time0, ps0, timeunit, psunit, ps_long_name) = read_met(lst0[0], 'atmospheric_pressure')
            time = np.hstack((time, time0/86400. + cc))
            ps = np.hstack((ps, ps0))
        ps[ps <= -999] = np.nan

    if len(timem) != len(time):
        raise ValueError("model and obs have inconsistent size")
        
    datamask = (ps-psm) > 10
    return(datamask)