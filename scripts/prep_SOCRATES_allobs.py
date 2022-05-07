# -*- coding: utf-8 -*-

"""
prepare all observational data from SOCRATES
inlcude:
    prep_SOCRATES_flight
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_SOCRATES_flight as air

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# flight data path
iwgpath = '/global/cscratch1/sd/sqtang/EAGLES/ESMAC_Diags_v2/raw_data/obs/SOCRATES/aircraft/mei-iwg1/'
RFpath = '../raw_data/obs/SOCRATES/aircraft/aircraft_lowrate/'
ccnpath = '../raw_data/obs/SOCRATES/aircraft/CCN/'

# output data path
prep_data_path = '../prep_data/SOCRATES/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare flight data; output time in 1min (dt=60s) resolution
print('prepare aircraft data:')
air.prep_CCN(ccnpath, RFpath, prep_data_path+'flight/', dt=60)        # CCN number concentration
air.prep_CNsize(RFpath, prep_data_path+'flight/', dt=60)    # aerosol size distribution
air.prep_CN(RFpath, prep_data_path+'flight/', dt=60)        # CN number concentration
air.prep_LWC(RFpath, prep_data_path+'flight/', dt=60)       # cloud liquid water content
air.prep_Nd(RFpath, prep_data_path+'flight/', dt=60)        # cloud size distribution