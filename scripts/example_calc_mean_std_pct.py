"""
script to calculate mean, standard deviation and percentiles

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.calc_statistics as calc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
# set output path for plots
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_obs_path + 'sfc_CPC_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cpc = obsdata['time'].load()
cpc10 = obsdata['cpc10'].load()
cpc3 = obsdata['cpc3'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
ncn10_e3sm2 = modeldata['NCN10'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

# trim for the same time period
IOP = 'IOP1'
time1 = np.datetime64('2016-04-25')
time2 = np.datetime64('2016-05-22')
time = pd.date_range(start='2016-04-25', end='2016-05-22', freq="H")
# IOP = 'IOP2'
# time1 = np.datetime64('2016-08-28')
# time2 = np.datetime64('2016-09-23')
# time = pd.date_range(start='2016-08-28', end='2016-09-23', freq="H")

cpc10 = cpc10[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
cpc3 = cpc3[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
ncn10_e3sm2 = ncn10_e3sm2[np.logical_and(time_m2>=time1, time_m2<=time2)]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# diagnose statistics
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
calc.mean_std_percentiles([cpc10,ncn10_e3sm2],['Obs','E3SMv2'],outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
