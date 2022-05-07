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
# set output path for plots
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data

# trim for the same time period
IOP = 'IOP1'
time1 = np.datetime64('2016-04-25')
time2 = np.datetime64('2016-05-22')
time = pd.date_range(start='2016-04-25', end='2016-05-22', freq="H")
# IOP = 'IOP2'
# time1 = np.datetime64('2016-08-28')
# time2 = np.datetime64('2016-09-23')
# time = pd.date_range(start='2016-08-28', end='2016-09-23', freq="H")

filename = prep_obs_path + 'sfc_CCN_'+site+'_'+IOP+'.nc'
obsdata = xr.open_dataset(filename)
time_ccn = obsdata['time'].load()
ccn5 = obsdata['CCN5'].load()
ccn2 = obsdata['CCN2'].load()
obsdata.close()
filename = prep_obs_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_nd = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
ccn2_m = modeldata['CCN4'].load()
ccn5_m = modeldata['CCN5'].load()
nd_m = modeldata['Nd_mean'].load()
modeldata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment


ccn2 = ccn2[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
ccn5 = ccn5[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
ndrop = ndrop[np.logical_and(time_nd>=time1, time_nd<=time2)]
ccn2_m = ccn2_m[np.logical_and(time_m>=time1, time_m<=time2)]
ccn5_m = ccn5_m[np.logical_and(time_m>=time1, time_m<=time2)]
nd_m = nd_m[np.logical_and(time_m>=time1, time_m<=time2)]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.jointhist([ccn2,ccn2_m], [ndrop,nd_m], 
                    xedges=np.arange(0,1200,100),yedges=np.arange(0,800,100),
                    xlabel='0.2%CCN (cm-3)', ylabel='Nd (cm-3)',
                   title=['Obs','E3SMv2'])
# ax.set_xscale('log')
fig.savefig(figpath+'jointhist_CCN2_Nd_sfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   

