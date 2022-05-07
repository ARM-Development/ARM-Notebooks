"""
example to generate timeseries plot

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import matplotlib.dates as mdates

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_obs_path = '../prep_data/'+site+'/surface/'
# set output path for plots
figpath= 'C:/Users/tang357/Downloads/figures/'+site+'/'
figpath= '../figures/'+site+'/surface/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_obs_path + 'sfc_ACSM_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_acsm = obsdata['time'].load()
org = obsdata['org'].load()
so4 = obsdata['so4'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
pom_m = modeldata['pom'].load()
mom_m = modeldata['mom'].load()
so4_m = modeldata['so4'].load()
soa_m = modeldata['soa'].load()
modeldata.close()
org_e3sm2 = pom_m + mom_m + soa_m
so4_e3sm2 = so4_m

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

org = org[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
so4 = so4[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
org_e3sm2 = org_e3sm2[np.logical_and(time_m2>=time1, time_m2<=time2)]
so4_e3sm2 = so4_e3sm2[np.logical_and(time_m2>=time1, time_m2<=time2)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.timeseries([time,time,], [org,org_e3sm2], 
                         legend=['ACSM','E3SMv2'], color=['k','r'],
                         xlabel='Time', title='Organic aerosol ($\mu$g/m$^3$) '+site+' '+IOP)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%M-%D'))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_org_sfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time,time,], [so4,so4_e3sm2], 
                         legend=['ACSM','E3SMv2'], color=['k','r'],
                         xlabel='Time', title='Sulfate aerosol ($\mu$g/m$^3$) '+site+' '+IOP)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%M-%D'))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
fig.savefig(figpath+'timeseries_so4_sfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   