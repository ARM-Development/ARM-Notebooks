"""
prepare all observational data from HISCALE
inlcude:
    prep_HISCALE_flight
    prep_HISCALE_sfc
    prep_HISCALE_satellite
"""

import numpy as np
import esmac_diags
import esmac_diags.preprocessing.prep_HISCALE_flight as air
import esmac_diags.preprocessing.prep_HISCALE_sfc as sfc
import esmac_diags.preprocessing.prep_HISCALE_satellite as sat


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set input data path
# flight data path
iwgpath = '../raw_data/obs/HISCALE/aircraft/mei-iwg1/'
amspath = '../raw_data/obs/HISCALE/aircraft/shilling-ams/'
beasdpath = '../raw_data/obs/HISCALE/aircraft/pekour-aafbe/'
ccnairpath = '../raw_data/obs/HISCALE/aircraft/mei-ccn/'
wcmpath = '../raw_data/obs/HISCALE/aircraft/matthews-wcm/'
fimspath = '../raw_data/obs/HISCALE/aircraft/wang-fims/'
pcasppath = '../raw_data/obs/HISCALE/aircraft/tomlinson-pcasp/'
cvipath = '../raw_data/obs/HISCALE/aircraft/pekour-cvi/'
cpcairpath = '../raw_data/obs/HISCALE/aircraft/mei-cpc/'
mergeSDpath = '../raw_data/obs/HISCALE/aircraft/mergedSD/'
# surface data path
acsmpath = '../raw_data/obs/HISCALE/surface/arm_acsm/'
armbepath = '../raw_data/obs/HISCALE/profile/armbe/'
arsclpath = '../raw_data/obs/HISCALE/profile/arscl/'
ccnsfcpath = '../raw_data/obs/HISCALE/surface/arm-ccn/'
cpcpath = '../raw_data/obs/HISCALE/surface/arm-cpc/'
cpcupath = '../raw_data/obs/HISCALE/surface/arm-cpcu/'
uhsaspath = '../raw_data/obs/HISCALE/surface/arm-uhsas/'
smpspath = '../raw_data/obs/HISCALE/surface/bnl-smps/'
nanosmpspath = '../raw_data/obs/HISCALE/surface/bnl-nanosmps/'
smps_pnnl_path = '../raw_data/obs/HISCALE/surface/pnnl-smps/'
mfrsrpath = '../raw_data/obs/HISCALE/surface/arm_mfrsr/'
ndroppath = '../raw_data/obs/HISCALE/surface/sgpndrop/'
# satellite data path
visstgridpath = '../raw_data/obs/HISCALE/satellite/visst/grid/'
visstpixpath = '../raw_data/obs/HISCALE/satellite/visst/pix_3x3/'

# output data path
prep_data_path = '../prep_data/HISCALE/'

# other settings
height_out = np.array([0.,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,\
                1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,\
                3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,\
                10000,10500,11000,11500,12000,12500,13000,14000,15000,16000,17000,18000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prepare flight data; output time in 1min (dt=60s) resolution
print('prepare aircraft data:')
air.prep_AMS(amspath, iwgpath, prep_data_path+'flight/', dt=60)       # aerosol composition
air.prep_beasd(beasdpath,iwgpath, prep_data_path+'flight/', dt=60)    # best estimate aerosol size distribution
air.prep_CCN(ccnairpath, iwgpath, prep_data_path+'flight/', dt=60)       # CCN number concentration
air.prep_CPC(cpcairpath, iwgpath, prep_data_path+'flight/', dt=60)       # aerosol number concentration (>3 or 10nm)
air.prep_PCASP100(pcasppath, iwgpath, prep_data_path+'flight/', dt=60)# aerosol number concentration (>100nm)
air.prep_mergeSD(mergeSDpath, iwgpath, prep_data_path+'flight/', dt=60) # merged cloud size distribution
air.prep_WCM(wcmpath, iwgpath, prep_data_path+'flight/', dt=60)       # cloud liquid water content
# prepare surface data. output time in 1hr (dt=3600s) resolution
print('prepare surface data:')
sfc.prep_ACSM(acsmpath, prep_data_path+'surface/', dt=3600)            # aerosol composition
sfc.prep_ccn(ccnsfcpath, prep_data_path+'surface/', dt=3600)              # CCN number concentration
sfc.prep_CPC(cpcpath, cpcupath, prep_data_path+'surface/', dt=3600)    # aerosol number concentration (>3 or 10nm)
sfc.prep_CNsize_UHSAS(uhsaspath, prep_data_path+'surface/', dt=3600)   # aerosol size distribution from UHSAS
sfc.prep_CNsize_SMPS_IOP1(smpspath, nanosmpspath, prep_data_path+'surface/', dt=3600)# aerosol size distribution from SMPS for IOP1
sfc.prep_CNsize_SMPS_IOP2(smps_pnnl_path, prep_data_path+'surface/', dt=3600)# aerosol size distribution from SMPS for IOP2
sfc.prep_cloud_2d(armbepath, prep_data_path+'surface/', height_out, dt=3600)   # 2D cloud fraction
sfc.prep_cloudheight_ARSCL(arsclpath, prep_data_path+'surface/', dt=3600)   # cloud height 
sfc.prep_totcld(armbepath, prep_data_path+'surface/', dt=3600)         # cloud fraction. from ARSCL, TSI and satellite sources
sfc.prep_LWP(armbepath, mfrsrpath, prep_data_path+'surface/', dt=3600) # cloud liquid water path
sfc.prep_Ndrop(ndroppath, prep_data_path+'surface/', dt=3600)          # cloud droplet number retrieval from ARM Ndrop VAP
sfc.prep_mfrsr_cod(mfrsrpath,  prep_data_path+'surface/', dt=3600)     # cloud optical depth from MFRSR
sfc.prep_mfrsr_Reff(mfrsrpath,  prep_data_path+'surface/', dt=3600)    # cloud effective radius from MFRSR
sfc.prep_radiation(armbepath, prep_data_path+'surface/', dt=3600)      # surface radiation
sfc.prep_LTS(armbepath, prep_data_path+'surface/', dt=3600)            # lower tropospheric stability
sfc.prep_precip(armbepath, prep_data_path+'surface/', dt=3600)         # surface precipitation
# prepare satellite data. output time in 1hr (dt=3600s) resolution
print('prepare satellite data:')
sat.prep_VISST_grid(visstgridpath, prep_data_path+'satellite/', dt=3600)     # VISST 0.5x0.5 degree gridded data
sat.prep_VISST_pixel(visstpixpath, prep_data_path+'satellite/', dt=3600)     # VISST 4km pixel-level data