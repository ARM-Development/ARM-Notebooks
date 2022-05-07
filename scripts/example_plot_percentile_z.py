"""
example to generate percentile plots in vertical coordinate
for aircraft data

"""
import os
import glob
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
prep_obs_path = '../prep_data/'+site+'/flight/'
# set output path for plots
figpath= '../figures/'+site+'/flight/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'CPC_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_cpc = obsdata['time'].load()
height = obsdata['height'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()

lst2 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst2)
time_m = modeldata['time'].load()
ncn10_m = modeldata['NCN10'].load()
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
height = height[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
ncn10_m = ncn10_m[np.logical_and(time_m>=time1, time_m<=time2)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)
    
height_bin = np.arange(300,3200,200)
fig,ax = plot.percentile_z([cpc10,ncn10_m], [height,height], 
                      height_bin, figsize=(3,8), title='CN (>10nm)',
                      xlabel='cm$^{-3}$', ylabel='height (m)', legend = ['Obs','E3SMv2'], )
ax.set_xticks([0,6000,12000])
fig.savefig(figpath+'percentile_z_CN10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
