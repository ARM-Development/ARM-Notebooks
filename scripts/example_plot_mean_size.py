"""
example to generate mean size distribution

"""
import os
import glob
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
prep_obs_path = '../prep_data/'+site+'/flight/'
# set output path for plots
figpath= '../figures/'+site+'/flight/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
lst = glob.glob(prep_obs_path + 'beasd_'+site+'_*.nc')
obsdata = xr.open_mfdataset(lst)
time_beasd = obsdata['time'].load()
size_beasd = obsdata['size'].load()
sizeh_beasd = obsdata['size_high'].load()
sizel_beasd = obsdata['size_low'].load()
beasd = obsdata['size_distribution_merged'].load()
obsdata.close()

lst2 = glob.glob(prep_model_path + 'E3SMv2_'+site+'_flight_*.nc')
modeldata = xr.open_mfdataset(lst2)
time_m2 = modeldata['time'].load()
ncn_e3sm2 = modeldata['NCNall'].load()
modeldata.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

# size distribution /dlogDp
dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
ncn_e3sm2 = ncn_e3sm2.T/dlogDp_e3sm
dlogDp_beasd = np.log10(sizeh_beasd/sizel_beasd)
beasd = beasd/dlogDp_beasd

# trim for the same time period
IOP = 'IOP1'
time1 = np.datetime64('2016-04-25')
time2 = np.datetime64('2016-05-22')
time = pd.date_range(start='2016-04-25', end='2016-05-22', freq="H")
# IOP = 'IOP2'
# time1 = np.datetime64('2016-08-28')
# time2 = np.datetime64('2016-09-23')
# time = pd.date_range(start='2016-08-28', end='2016-09-23', freq="H")

beasd = beasd[np.logical_and(time_beasd>=time1, time_beasd<=time2),:]
ncn_e3sm2 = ncn_e3sm2[np.logical_and(time_m2>=time1, time_m2<=time2),:]

pdf_beasd = np.nanmean(beasd, axis=0)
pdf_e3sm2 = np.nanmean(ncn_e3sm2, axis=0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot 
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.mean_size([size_beasd, np.arange(1,3001)], [pdf_beasd, pdf_e3sm2], legend = ['Obs','E3SMv2'], 
                  xlimit=(10, 1000), ylimit=(1e-2,1e4), 
                  xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
                 color=['k','r'], title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
fig.savefig(figpath+'AerosolSize_mean_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
