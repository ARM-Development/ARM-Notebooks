# -*- coding: utf-8 -*-

"""
prepare all observational data from CSET
inlcude:
    prep_CSET_flight
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_CSET_flight as air

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# flight data path
RFpath = '/global/cscratch1/sd/sqtang/EAGLES/ESMAC_Diags_v2/raw_data/obs/CSET/aircraft/aircraft_lowrate/'

# output data path
prep_data_path = '/global/cscratch1/sd/sqtang/EAGLES/ESMAC_Diags_v2/prep_data/CSET/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare flight data; output time in 1min (dt=60s) resolution
print('prepare aircraft data:')
air.prep_CNsize(RFpath, prep_data_path+'flight/', dt=60)    # aerosol size distribution
air.prep_CN(RFpath, prep_data_path+'flight/', dt=60)        # CN number concentration
air.prep_LWC(RFpath, prep_data_path+'flight/', dt=60)       # cloud liquid water content
air.prep_Nd(RFpath, prep_data_path+'flight/', dt=60)        # cloud size distribution