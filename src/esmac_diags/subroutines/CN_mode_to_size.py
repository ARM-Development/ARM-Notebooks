'''
  calculate aerosol number concentration with 1nm 
  increment from 1nm to 3000nm (dimension [size, level, time])

  rewritten from Kai Zhang and Jian Sun's NCL code
  cutoff_CTRL_regions_0_3000nm.ncl
  
  Shuaiqi Tang
  2021.4.30

'''

import numpy as np
import scipy.special

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def func_cutoff_NN(num,dn,lnsg,dmax,dmin):
    """
    function of calculate total number concentration between two cutoff size
    in a log-normal distribution
    
    Parameters
    ----------
    num : total number concentration of the log-normal distribution
    dn : median size of the distribution
    lnsg: log value of geometric standard deviation of the distribution
    dmax: maximum value of cutoff size
    dmin: minimum value of cutoff size

    Returns
    -------
    nc: total number between dmin and dmax
    
    """
    logdr = np.log(dmax/dn)
    logdl = np.log(dmin/dn)
    # erf is the error function
    erfr = scipy.special.erf(logdr/(np.sqrt(2)*lnsg))
    erfl = scipy.special.erf(logdl/(np.sqrt(2)*lnsg))
    nc = 0.5*num*(erfr-erfl)
    return(nc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
def calc_CNsize_cutoff_0_3000nm(dnall,numall,T,P):
    """
    main function to calculate CN size of 1nm bins from 1-3000nm in
    MAM4 or MAM5 (with nucleation mode) aerosol modes in CESM or E3SM
    
    Parameters
    ----------
    numall : [num_a1, num_a2, num_a3, num_a4]
            number concentration in each mode in #/kg (num_a*)
    dnall : [dgnd_a01, dgnd_a02, dgnd_a03, dgnd_a04]
            number median diameter of each mode in m (dgnd_a*)
    T: temperature, in K
    P: pressure, in Pa

    Returns
    -------
    ncut: numpy array of [3000, size_of_num_a1]
            aerosol number in 1-3000nm bins with 1nm increment
    
    """
    
    nmode = len(dnall)
    if nmode!=4 and nmode!=5:
        raise ValueError("Error: Currently only apply for MAM4 and MAM5 with 4 or 5 modes")
        
    #%% set constant and other variables
    
    # SHR_CONST_STEBOL  = 5.67e-8      # Stefan-Boltzmann constant ~ W/m^2/K^4
    SHR_CONST_BOLTZ   = 1.38065e-23  # Boltzmann's constant ~ J/K/molecule
    SHR_CONST_AVOGAD  = 6.02214e26   # Avogadro's number ~ molecules/kmole
    SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ       # Universal gas constant ~ J/K/kmole
    SHR_CONST_MWDAIR  = 28.966     # molecular weight dry air ~ kg/kmole
    # SHR_CONST_MWWV    = 18.016     # molecular weight water vapor
    SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR        # Dry air gas constant     ~ J/K/kg
    # pi = 3.1415926
    
    # calculate rho
    Rho = P / (SHR_CONST_RDAIR * T)   # kg/m3 
    
    # set geometric standard deviation of each mode
    Sg = [1.8, 1.6, 1.8, 1.6, 1.6] # [Acc, Aik, Coa, PrC, Nuc]
    lnSg = np.log(Sg)
    
    # set cutoff size with 1nm increment from 0 to 3000nm
    cnsize = np.arange(1,3001)*1e-9
    dminx = cnsize - 1e-9# lower bound
    dminx[0] = 1e-33      # set as non-zero small value
    dmaxx = cnsize       # upper bound
    ns = len(cnsize)  # dimension length
    
    
    #%% dry diameter and number concentration in each mode
    dn1=dnall[0]
    dn2=dnall[1]
    dn3=dnall[2]
    dn4=dnall[3]
    num1=numall[0]
    num2=numall[1]
    num3=numall[2]
    num4=numall[3]
    if nmode==5:
        dn5=dnall[4]
        num5=numall[4]
    
    # unit convert to #/m3
    num1 = num1*Rho
    num2 = num2*Rho
    num3 = num3*Rho
    num4 = num4*Rho
    if nmode==5:
        num5=num5*Rho
    
    #%% calculate aerosol size in each size bin
    ncut = np.full(tuple([ns]+list(dn1.shape)),np.nan)
    for ss in range(ns):
        ncut1 = func_cutoff_NN(num1,dn1,lnSg[0],dmaxx[ss],dminx[ss])
        ncut2 = func_cutoff_NN(num2,dn2,lnSg[1],dmaxx[ss],dminx[ss])
        ncut3 = func_cutoff_NN(num3,dn3,lnSg[2],dmaxx[ss],dminx[ss])
        ncut4 = func_cutoff_NN(num4,dn4,lnSg[3],dmaxx[ss],dminx[ss])
        if len(ncut.shape)==1:
            ncut[ss]=ncut1+ncut2+ncut3+ncut4
        else:
            ncut[ss,:]=ncut1+ncut2+ncut3+ncut4
    if nmode==5:
        for ss in range(ns):
            ncut5 = func_cutoff_NN(num5,dn5,lnSg[4],dmaxx[ss],dminx[ss])
            if len(ncut.shape)==1:
                ncut[ss]=ncut[ss]+ncut5
            else:
                ncut[ss,:]=ncut[ss,:]+ncut5
    
    
    # calculate total aerosol number and relative error
    if nmode==4:
        NN = num1+num2+num3+num4
    elif nmode==5:
        NN = num1+num2+num3+num4+num5
    
    NN2=np.sum(ncut,0)
    diff = np.abs(NN2-NN)/NN
    
    if np.nanmax(diff)>0.01:
        print(['maximum absolute relative error: '+str(100*np.nanmax(diff))+'%'])

    
    return(ncut)
