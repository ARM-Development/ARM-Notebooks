"""
script to run a test case
compare the figures generated at testcase/figures/ with testcase/figures_verify
to makesure testcase works as expected
"""

import os
import xarray as xr
import esmac_diags.plotting.plot_esmac_diags as plot
import matplotlib.dates as mdates

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../testcase/data/model/'
prep_obs_path = '../testcase/data/obs/'
# set output path for plots
figpath= '../testcase/figures/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_obs_path + 'AMS_HISCALE_20160920b.nc'
obsdata = xr.open_dataset(filename)
time_o = obsdata['time'].load()
org = obsdata['ORG'].load()
so4 = obsdata['SO4'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv2_HISCALE_flight_20160920b.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
pom_m = modeldata['pom'].load()
mom_m = modeldata['mom'].load()
so4_m = modeldata['so4'].load()
soa_m = modeldata['soa'].load()
modeldata.close()
org_e3sm2 = pom_m + mom_m + soa_m
so4_e3sm2 = so4_m

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# output plot
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.timeseries([time_o,time_m,], [org,org_e3sm2], 
                         legend=['AMS','E3SMv2'], color=['k','r'],
                         xlabel='Time (UTC) in 2016-09-20', title='Organic aerosol ($\mu$g/m$^3$) '+site)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
fig.savefig(figpath+'timeseries_organic_'+site+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

