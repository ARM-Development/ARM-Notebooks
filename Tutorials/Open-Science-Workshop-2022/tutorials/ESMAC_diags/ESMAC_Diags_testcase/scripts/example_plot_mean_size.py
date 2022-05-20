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

filename = prep_obs_path + 'sfc_SMPS_'+site+'_'+IOP+'.nc'
obsdata = xr.open_dataset(filename)
time_smps = obsdata['time'].load()
smpsall = obsdata['dN_dlogDp'].load()
size_smps = obsdata['size'].load()
obsdata.close()

filename = prep_model_path + 'E3SMv1_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m1 = modeldata['time'].load()
ncn_e3sm1 = modeldata['NCNall'].load()
modeldata.close()

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
ncn_e3sm2 = modeldata['NCNall'].load()
modeldata.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatment

# convert to dN/dlog10Dp
# SMPS is already dlogDp
dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
ncn_e3sm1 = ncn_e3sm1.T/dlogDp_e3sm
ncn_e3sm2 = ncn_e3sm2.T/dlogDp_e3sm

pdf_smps = np.nanmean(smpsall, axis=0)
pdf_smps[pdf_smps<1e-4] = np.nan
pdf_e3sm1 = np.nanmean(ncn_e3sm1, axis=0)
pdf_e3sm2 = np.nanmean(ncn_e3sm2, axis=0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot 
if not os.path.exists(figpath):
    os.makedirs(figpath)

fig,ax = plot.mean_size([size_smps, np.arange(1,3001), np.arange(1,3001)], [pdf_smps, pdf_e3sm1, pdf_e3sm2], 
                  legend = ['SMPS','E3SMv1','E3SMv2'],  color=['k','r','b'], 
                  xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', xlimit=(2, 3000), ylimit=(1e-1,5e4), 
                 title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
fig.savefig(figpath+'AerosolSize_mean_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# show figures in interactive commandline screen
import matplotlib.pyplot as plt
plt.show()   
