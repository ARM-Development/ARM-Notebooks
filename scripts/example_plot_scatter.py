"""
example to generate scatter plot

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
from esmac_diags.subroutines.time_resolution_change import avg_time_1d
import matplotlib.dates as mdates


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_obs_path = '../prep_data/'+site+'/flight/'
prep_obs_path_2 = '../prep_data/'+site+'/surface/'
# set output path for plots
figpath= '../figures/'+site+'/flight/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'mergedSD_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst,concat_dim='time',combine='nested')
time_air = obsdata['time'].load()
nd_air = obsdata['Nd'].load()
obsdata.close()
nd_air = nd_air*0.001   # #/L to #/cm3

obsdata = xr.open_dataset(prep_obs_path_2 + 'Ndrop_'+site+'.nc')
time_sfc = obsdata['time'].load()
nd_sfc = obsdata['cdnc'].load()
obsdata.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

nd_air[nd_air<1]=np.nan
nd_sfc[nd_sfc<1]=np.nan

nd_sfc_1min = np.interp(time_air, time_sfc, nd_sfc)
nd_air_1hr = avg_time_1d(time_air, nd_air, time_sfc)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.scatter(nd_sfc_1min, nd_air.data, xlimit=(0,600), ylimit=(0,600),
                    xlabel='Ndrop', ylabel='flight', title='Nd (cm$^{-3}$) '+site,
                    linear_fit=True, intercept=False)
fig,ax = plot.scatter(nd_sfc.data, nd_air_1hr, xlimit=(0,600), ylimit=(0,600),
                    xlabel='Ndrop', ylabel='flight', title='Nd (cm$^{-3}$) '+site,
                linear_fit=True, intercept=False)
fig.savefig(figpath+'Nd_compare_flight_sfc_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
