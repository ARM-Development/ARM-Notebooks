# -*- coding: utf-8 -*-

"""
prepare all observational data from MAGIC
inlcude:
    prep_MAGIC_ship
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_MAGIC_ship as ship

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# ship data path
shipmetpath = '../raw_data/obs/MAGIC/ship/magmarinemetM1.b1/'
mwrpath = '../raw_data/obs/MAGIC/ship/magmwrret1liljclouM1.s2/'
cpcpath = '../raw_data/obs/MAGIC/ship/magaoscpcfM1.a1/'
ccnpath = '../raw_data/obs/MAGIC/ship/magaosccn100M1.a1/'
uhsaspath = '../raw_data/obs/MAGIC/ship/magaosuhsasM1.a1/'
Ndpath = '../raw_data/obs/MAGIC/ship/Cloud_Micro_Retrieval/'

# output data path
prep_data_path = '../prep_data/MAGIC/'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare ship data. output time in 1hr (dt=3600s) resolution
print('prepare ship data:')
ship.prep_CCN(shipmetpath, ccnpath, prep_data_path+'ship/', dt=3600)              # CCN number concentration
ship.prep_CN(shipmetpath, cpcpath, uhsaspath, prep_data_path+'ship/', dt=3600)    # aerosol number concentration (>3 or 10nm)
ship.prep_CNsize(shipmetpath, uhsaspath, prep_data_path+'ship/', dt=3600)   # aerosol size distribution from UHSAS
ship.prep_MWR(shipmetpath, mwrpath, prep_data_path+'ship/', dt=3600) # cloud liquid water path
ship.prep_Nd_Wu_etal(Ndpath, prep_data_path+'ship/', dt=3600)          # cloud droplet number retrieval from Wu et al.