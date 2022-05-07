"""
example to generate 1-d histogram plot

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
prep_obs_path_2 = '../prep_data/'+site+'/satellite/'
# set output path for plots
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_obs_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_arm = obsdata['time'].load()
nd_arm = obsdata['cdnc'].load()
obsdata.close()

filename = prep_obs_path_2 + 'Nd_VISSTgrid_'+site+'.nc'
satdata = xr.open_dataset(filename)
time_sat = satdata['time'].load()
nd_sat = satdata['Nd'].load()
satdata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
nd_e3sm2 = modeldata['Nd_mean'].load()
modeldata.close()

nd_e3sm2[nd_e3sm2>5000] = np.nan
nd_e3sm2[nd_e3sm2<1] = np.nan

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

nd_arm = nd_arm[np.logical_and(time_arm>=time1, time_arm<=time2)]
nd_sat = nd_sat[np.logical_and(time_sat>=time1, time_sat<=time2)]
nd_e3sm2 = nd_e3sm2[np.logical_and(time_m2>=time1, time_m2<=time2)]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)

w1 = np.ones_like(nd_arm)/sum(~np.isnan(nd_arm.data))
w2 = np.ones_like(nd_sat)/sum(~np.isnan(nd_sat.data))
w3 = np.ones_like(nd_e3sm2)/sum(~np.isnan(nd_e3sm2.data))
fig,ax = plot.hist([nd_arm, nd_sat, nd_e3sm2], bins=np.arange(10,500,20), weights=[w1, w2, w3], \
                   legend=['ARM Ndrop','VISST', 'E3SMv2'], color=['k','b','r'],
                   xlabel='cm$^{-3}$', ylabel='Fraction', title=None)
fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   