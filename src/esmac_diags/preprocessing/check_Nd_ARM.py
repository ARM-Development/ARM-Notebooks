"""
check Nd_ARM retrieval and compare to Ndrop data
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import median_time_1d, median_time_2d,\
                avg_time_1d, avg_time_2d
from esmac_diags.subroutines.quality_control import qc_remove_neg, qc_mask_qcflag, \
                qc_mask_qcflag_cpc
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_ARM


site='ACEENA'
arsclpath = '../../../data/ACEENA/obs/profile/arscl/'
mfrsrpath = '../../../data/ACEENA/obs/surface/arm_mfrsr/'
ndroppath = '../../../data/ACEENA/obs/surface/enandrop/'
date='20180211'

# site='HISCALE'
# date='20160831'
# arsclpath = '../../../data/HISCALE/obs/profile/arscl/'
# mfrsrpath = '../../../data/HISCALE/obs/surface/arm_mfrsr/'
# ndroppath = '../../../data/HISCALE/obs/surface/sgpndrop/'


#%% read in data
lst = glob.glob(os.path.join(ndroppath, '*ndropmfrsrC1.c1.'+date+'*.nc'))
lst.sort()
obsdata = xr.open_mfdataset(lst, combine='by_coords')
time0 = obsdata['time'].load()
ndrop = obsdata['drop_number_conc'].load()
qc_nd = obsdata['qc_drop_number_conc'].load()
obsdata.close()
# quality control
ndrop = qc_mask_qcflag(ndrop,qc_nd)
ndrop = ndrop*1e-6   # m-3 to cm-3


#%% read in data
lst1 = glob.glob(os.path.join(arsclpath, '*arsclkazrbnd1kolliasC1.c0.'+date+'*.nc'))
lst1.sort()
arscldata = xr.open_mfdataset(lst1, combine='by_coords')
arscltime = arscldata['time']
cbh = arscldata['cloud_base_best_estimate']
cbhs = arscldata['cloud_layer_base_height']
cths = arscldata['cloud_layer_top_height']
arscldata.close()

lst2 = glob.glob(os.path.join(mfrsrpath, '*mfrsrcldod1minC1.c1.'+date+'*.cdf'))
lst2.sort()
# first data
mfrsrdata = xr.open_dataset(lst2[0])
mfrsrtime = mfrsrdata['time']
lwp = mfrsrdata['lwp']
qc_lwp = mfrsrdata['qc_lwp']
cod = mfrsrdata['optical_depth_instantaneous']
qc_cod = mfrsrdata['qc_optical_depth_instantaneous']
mfrsrdata.close()
for file in lst2[1:]:
    mfrsrdata = xr.open_dataset(file)
    mfrsrtime = xr.concat([mfrsrtime, mfrsrdata['time']], dim="time")
    lwp = xr.concat([lwp, mfrsrdata['lwp']], dim="time")
    qc_lwp = xr.concat([qc_lwp, mfrsrdata['qc_lwp']], dim="time")
    cod = xr.concat([cod, mfrsrdata['optical_depth_instantaneous']], dim="time")
    qc_cod = xr.concat([qc_cod, mfrsrdata['qc_optical_depth_instantaneous']], dim="time")
    mfrsrdata.close()

lwp.load()
qc_lwp.load()
cod.load()
qc_cod.load()
cbh.load()
cbhs.load()
cths.load()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# calculate CDNC
# calculate cloud depth for the lowest cloud layer in ARSCL
cth = cths[:,0]
H = cths[:,0] - cbh
H.load()
H = qc_remove_neg(H, remove_zero='true')
# H[cbhs[:,1] > 0] = np.nan

# filter some data samples
# H[cbhs[:,0] > 5000] = np.nan
# H[cth<100] = np.nan
# H[cth>5000.] = np.nan   # remove deep clouds with cloud top >5km
# lwp[np.logical_or(lwp<0.02, lwp>0.3)] = np.nan
# cod[np.logical_or(cod<2, cod>60)] = np.nan

# calculate CDNC first then average into 1hr
time = mfrsrtime.data
H_tmp = np.interp(np.int64(time), np.int64(arscltime), H)
nd = calc_cdnc_ARM(lwp, cod, H_tmp)

#%%
# nd[nd>500]=np.nan
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
fig,ax = plt.subplots(1,1,figsize=(8,2))
ax.plot(time0,ndrop,'k',label='Ndrop')
ax.plot(time,nd,'r.',linewidth=2,label='Nd_ARM')
ax.legend(loc='upper right', fontsize='large')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax.set_xlabel('time in '+date)
ax.set_ylim(-10,1000)
x = ax.get_xlim()

fig,ax = plt.subplots(1,1,figsize=(8,2))
# ax.plot(time,lwp*100,'k',label='LWP*100')
# ax.plot(time,cod,'b',label='cod')
# ax.plot(time,H_tmp*0.001,'r',label='H_cld')
# ax.legend(loc='upper right')
# ax.plot(arscltime,cbhs[:,0]*0.001,'m')
ax.plot(arscltime,cbhs[:,0]*0.001,'.',label='1')
ax.plot(arscltime,cbhs[:,1]*0.001,'.',label='2')
ax.plot(arscltime,cbhs[:,2]*0.001,'.',label='3')
ax.plot(arscltime,cbhs[:,3]*0.001,'.',label='4')
ax.plot(arscltime,cbhs[:,4]*0.001,'.',label='5')
ax.legend(loc='upper left',fontsize='x-small')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.set_xlabel('time in '+date)
ax.set_xlim(x)
