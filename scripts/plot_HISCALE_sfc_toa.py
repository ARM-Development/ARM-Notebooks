"""
script to generate all plots for HISCALE surface data

"""
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import esmac_diags.plotting.plot_esmac_diags as plot
import esmac_diags.plotting.calc_statistics as calc
import matplotlib.dates as mdates

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# settings
# set site name and datapath

# set site name.
site = 'HISCALE'

# path of prepared files
prep_model_path = '../prep_data/'+site+'/model/'
prep_sfc_path = '../prep_data/'+site+'/surface/'
prep_sat_path = '../prep_data/'+site+'/satellite/'
# path of output figures
figpath= '../figures/'+site+'/sfc_toa/'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in data
filename = prep_sfc_path + 'sfc_ACSM_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_acsm = obsdata['time'].load()
org = obsdata['org'].load()
so4 = obsdata['so4'].load()
nh4 = obsdata['nh4'].load()
no3 = obsdata['no3'].load()
chl = obsdata['chl'].load()
obsdata.close()

lst = sorted(glob.glob(prep_sfc_path + 'sfc_CCN_'+site+'_*.nc'))
obsdata = xr.open_mfdataset(lst)
time_ccn = obsdata['time'].load()
ccn1 = obsdata['CCN1'].load()
ccn2 = obsdata['CCN2'].load()
ccn5 = obsdata['CCN5'].load()
obsdata.close()

filename = prep_sfc_path + 'sfc_CPC_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cpc = obsdata['time'].load()
cpc3 = obsdata['cpc3'].load()
cpc10 = obsdata['cpc10'].load()
obsdata.close()

filename = prep_sfc_path + 'cod_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_cod = obsdata['time'].load()
cod = obsdata['cod'].load()
obsdata.close()

filename = prep_sfc_path + 'LWP_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_lwp = obsdata['time'].load()
lwp_armbe = obsdata['lwp_armbe'].load()
lwp_mfrsr = obsdata['lwp_mfrsr'].load()
obsdata.close()

filename = prep_sfc_path + 'LTS_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_lts = obsdata['time'].load()
LTS700 = obsdata['LTS700'].load()
LTS850 = obsdata['LTS850'].load()
obsdata.close()

# filename = prep_sfc_path + 'Nd_ARMretrieval_'+site+'.nc'
# obsdata = xr.open_dataset(filename)
# time_nd_arm = obsdata['time'].load()
# Nd_arm = obsdata['cdnc'].load()
# obsdata.close()

filename = prep_sfc_path + 'Ndrop_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_ndrop = obsdata['time'].load()
ndrop = obsdata['cdnc'].load()
obsdata.close()

filename = prep_sfc_path + 'reff_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_reff = obsdata['time'].load()
reff = obsdata['reff'].load()
obsdata.close()

filename = prep_sfc_path + 'precip_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_pr = obsdata['time'].load()
precip = obsdata['precip'].load()
obsdata.close()

filename = prep_sfc_path + 'sfc_radiation_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_rad = obsdata['time'].load()
lwdnsfc = obsdata['lwdn'].load()
swdnsfc = obsdata['swdn'].load()
lwupsfc = obsdata['lwup'].load()
swupsfc = obsdata['swup'].load()
obsdata.close()

filename = prep_sfc_path + 'sfc_SMPS_'+site+'_IOP1.nc'
obsdata = xr.open_dataset(filename)
time1 = obsdata['time'].load()
smpsall_1 = obsdata['dN_dlogDp'].load()
smps100_1 = obsdata['smps100_dlogDp'].load()
size1 = obsdata['size'].load()
obsdata.close()
filename = prep_sfc_path + 'sfc_SMPS_'+site+'_IOP2.nc'
obsdata = xr.open_dataset(filename)
time2 = obsdata['time'].load()
smpsall_2 = obsdata['dN_dlogDp'].load()
smps100_2 = obsdata['smps100_dlogDp'].load()
size2 = obsdata['size'].load()
obsdata.close()
time_smps = xr.concat((time1,time2),dim='time')
smpsall_2_resize = np.full((smpsall_2.shape[0], smpsall_1.shape[1]), np.nan)
smpsall_2_resize[:, 73:183] = smpsall_2
smps_all = np.vstack((smpsall_1.data, smpsall_2_resize))
size_smps = size1
smps_all = xr.DataArray(data=smps_all, coords=dict(time=time_smps,size=size_smps), attrs=smpsall_1.attrs)
smps100 = xr.concat((smps100_1,smps100_2),dim='time')

filename = prep_sfc_path + 'totcld_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_totcld = obsdata['time'].load()
cld_arscl = obsdata['tot_cld_arscl'].load()
cld_tsi = obsdata['tot_cld_tsi'].load()
cld_visst = obsdata['tot_cld_visst'].load()
obsdata.close()

filename = prep_sfc_path + 'sfc_UHSAS_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_uhsas = obsdata['time'].load()
uhsas_all = obsdata['uhsas_all'].load()
uhsas100 = obsdata['uhsas100'].load()
size_uhsas = obsdata['size'].load()
sizel_uhsas = obsdata['size_low'].load()
sizeh_uhsas = obsdata['size_high'].load()
obsdata.close()

# satellite data
filename = prep_sat_path + 'cod_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
cod_sat = obsdata['cod'].load()
obsdata.close()

filename = prep_sat_path + 'LWP_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwp_sat = obsdata['lwp'].load()
obsdata.close()

filename = prep_sat_path + 'Nd_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
nd_sat = obsdata['Nd'].load()
obsdata.close()


filename = prep_sat_path + 'Reff_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
reff_sat = obsdata['reff'].load()
obsdata.close()

filename = prep_sat_path + 'lwflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
lwnettoa = obsdata['lwnettoa'].load()
obsdata.close()

filename = prep_sat_path + 'swflx_VISSTgrid_'+site+'.nc'
obsdata = xr.open_dataset(filename)
time_visst = obsdata['time'].load()
swnettoa = obsdata['swnettoa'].load()
obsdata.close()

# E3SM data
filename = prep_model_path + 'E3SMv1_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m = modeldata['time'].load()
bc_m = modeldata['bc'].load()
dst_m = modeldata['dst'].load()
ncl_m = modeldata['ncl'].load()
pom_m = modeldata['pom'].load()
mom_m = modeldata['mom'].load()
so4_m = modeldata['so4'].load()
soa_m = modeldata['soa'].load()
ccn1_m = modeldata['CCN3'].load()
ccn2_m = modeldata['CCN4'].load()
ccn5_m = modeldata['CCN5'].load()
ncn3_m = modeldata['NCN3'].load()
ncn10_m = modeldata['NCN10'].load()
ncn100_m = modeldata['NCN100'].load()
ncn_m = modeldata['NCNall'].load()
cod_m = modeldata['cod'].load()
reff_m = modeldata['reff'].load()
lwp_m = modeldata['TGCLDLWP'].load()
nd_m = modeldata['Nd_mean'].load()
precip_m = modeldata['PRECT'].load()
cld_m = modeldata['CLDTOT'].load()
lwdnsfc_m = modeldata['FLDS'].load()
lwnetsfc_m = modeldata['FLNS'].load()
lwnettoa_m = modeldata['FLNT'].load()
lwuptoa_m = modeldata['FLUT'].load()
swdnsfc_m = modeldata['FSDS'].load()
swnetsfc_m = modeldata['FSNS'].load()
swdntoa_m = modeldata['SOLIN'].load()
swnettoa_m = modeldata['FSNT'].load()
swuptoa_m = modeldata['FSUTOA'].load()
modeldata.close()
# lwupsfc_m = lwnetsfc_m + lwdnsfc_m
# swupsfc_m = swdnsfc_m - swnetsfc_m
org_m = pom_m + mom_m + soa_m

filename = prep_model_path + 'E3SMv2_'+site+'_sfc.nc'
modeldata = xr.open_dataset(filename)
time_m2 = modeldata['time'].load()
bc_m2 = modeldata['bc'].load()
dst_m2 = modeldata['dst'].load()
ncl_m2 = modeldata['ncl'].load()
pom_m2 = modeldata['pom'].load()
mom_m2 = modeldata['mom'].load()
so4_m2 = modeldata['so4'].load()
soa_m2 = modeldata['soa'].load()
ccn1_m2 = modeldata['CCN3'].load()
ccn2_m2 = modeldata['CCN4'].load()
ccn5_m2 = modeldata['CCN5'].load()
ncn3_m2 = modeldata['NCN3'].load()
ncn10_m2 = modeldata['NCN10'].load()
ncn100_m2 = modeldata['NCN100'].load()
ncn_m2 = modeldata['NCNall'].load()
cod_m2 = modeldata['cod'].load()
reff_m2 = modeldata['reff'].load()
lwp_m2 = modeldata['TGCLDLWP'].load()
nd_m2 = modeldata['Nd_mean'].load()
precip_m2 = modeldata['PRECT'].load()
cld_m2 = modeldata['CLDTOT'].load()
lwdnsfc_m2 = modeldata['FLDS'].load()
lwnetsfc_m2 = modeldata['FLNS'].load()
lwnettoa_m2 = modeldata['FLNT'].load()
lwuptoa_m2 = modeldata['FLUT'].load()
swdnsfc_m2 = modeldata['FSDS'].load()
swnetsfc_m2 = modeldata['FSNS'].load()
swdntoa_m2 = modeldata['SOLIN'].load()
swnettoa_m2 = modeldata['FSNT'].load()
swuptoa_m2 = modeldata['FSUTOA'].load()
modeldata.close()
# lwupsfc_m2 = lwnetsfc_m2 + lwdnsfc_m2
# swupsfc_m2 = swdnsfc_m2 - swnetsfc_m2
org_m2 = pom_m2 + mom_m2 + soa_m2


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# specific data treatments

# divided by dlogDp in size distribution
dlogDp_uhsas = np.log10(sizeh_uhsas/sizel_uhsas)
uhsas_all = uhsas_all/dlogDp_uhsas

dlogDp_e3sm = np.log10(np.arange(2,3002)/np.arange(1,3001))
ncn_m = ncn_m.T/dlogDp_e3sm
ncn_m2 = ncn_m2.T/dlogDp_e3sm

# SMPS data is already dN/dlogDp, total number concentration must multiply by dlogDp
dlogDp_smps = np.log10(size_smps[1:].data/size_smps[0:-1].data)
smps100 = smps100 * np.mean(dlogDp_smps)

# remove LWP less than 20 g/m2
lwp_armbe[lwp_armbe<20]=np.nan
# lwp_armbe[lwp_armbe>200]=np.nan
lwp_mfrsr[lwp_mfrsr<20]=np.nan
# lwp_mfrsr[lwp_mfrsr>200]=np.nan
lwp_sat[lwp_sat<20]=np.nan
lwp_m[lwp_m<20]=np.nan
lwp_m2[lwp_m2<20]=np.nan

# remove Nd less than 1 cm-3
ndrop[ndrop<11] = np.nan
nd_sat[nd_sat<11] = np.nan
nd_m[nd_m<11] = np.nan
nd_m2[nd_m2<11] = np.nan

# remove Nd greater than 2000
ndrop[ndrop>2000] = np.nan
nd_m[nd_m>2000] = np.nan
nd_m2[nd_m2>2000] = np.nan

# trim for the same time period
# IOP = 'IOP1'
# time1 = np.datetime64('2016-04-25')
# time2 = np.datetime64('2016-05-22')
# time = pd.date_range(start='2016-04-25', end='2016-05-22', freq="H")
IOP = 'IOP2'
time1 = np.datetime64('2016-08-28')
time2 = np.datetime64('2016-09-23')
time = pd.date_range(start='2016-08-28', end='2016-09-23', freq="H")

# surface
org = org[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
so4 = so4[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
nh4 = nh4[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
no3 = no3[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
chl = chl[np.logical_and(time_acsm>=time1, time_acsm<=time2)]
ccn1 = ccn1[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
ccn2 = ccn2[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
ccn5 = ccn5[np.logical_and(time_ccn>=time1, time_ccn<=time2)]
cpc3 = cpc3[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
cpc10 = cpc10[np.logical_and(time_cpc>=time1, time_cpc<=time2)]
cod = cod[np.logical_and(time_cod>=time1, time_cod<=time2)]
lwp_armbe = lwp_armbe[np.logical_and(time_lwp>=time1, time_lwp<=time2)]
lwp_mfrsr = lwp_mfrsr[np.logical_and(time_lwp>=time1, time_lwp<=time2)]
LTS700 = LTS700[np.logical_and(time_lts>=time1, time_lts<=time2)]
LTS850 = LTS850[np.logical_and(time_lts>=time1, time_lts<=time2)]
ndrop = ndrop[np.logical_and(time_ndrop>=time1, time_ndrop<=time2)]
reff = reff[np.logical_and(time_reff>=time1, time_reff<=time2)]
precip = precip[np.logical_and(time_pr>=time1, time_pr<=time2)]
lwdnsfc = lwdnsfc[np.logical_and(time_rad>=time1, time_rad<=time2)]
swdnsfc = swdnsfc[np.logical_and(time_rad>=time1, time_rad<=time2)]
lwupsfc = lwupsfc[np.logical_and(time_rad>=time1, time_rad<=time2)]
swupsfc = swupsfc[np.logical_and(time_rad>=time1, time_rad<=time2)]
smps100 = smps100[np.logical_and(time_smps>=time1, time_smps<=time2)]
uhsas100 = uhsas100[np.logical_and(time_uhsas>=time1, time_uhsas<=time2)]
cld_arscl = cld_arscl[np.logical_and(time_totcld>=time1, time_totcld<=time2)]
cld_tsi = cld_tsi[np.logical_and(time_totcld>=time1, time_totcld<=time2)]
cld_visst = cld_visst[np.logical_and(time_totcld>=time1, time_totcld<=time2)]
lwnetsfc = lwupsfc - lwdnsfc
swnetsfc = swdnsfc - swupsfc

uhsas_all = uhsas_all[np.logical_and(time_uhsas>=time1, time_uhsas<=time2), :]
smps_all = smps_all[np.logical_and(time_smps>=time1, time_smps<=time2), :]

# satellite
cod_sat = cod_sat[np.logical_and(time_visst>=time1, time_visst<=time2)]
lwp_sat = lwp_sat[np.logical_and(time_visst>=time1, time_visst<=time2)]
nd_sat = nd_sat[np.logical_and(time_visst>=time1, time_visst<=time2)]
reff_sat = reff_sat[np.logical_and(time_visst>=time1, time_visst<=time2)]
lwnettoa = lwnettoa[np.logical_and(time_visst>=time1, time_visst<=time2)]
swnettoa = swnettoa[np.logical_and(time_visst>=time1, time_visst<=time2)]

# E3SMv1 data
bc_m = bc_m[np.logical_and(time_m>=time1, time_m<=time2)]
dst_m = dst_m[np.logical_and(time_m>=time1, time_m<=time2)]
ncl_m = ncl_m[np.logical_and(time_m>=time1, time_m<=time2)]
so4_m = so4_m[np.logical_and(time_m>=time1, time_m<=time2)]
org_m = org_m[np.logical_and(time_m>=time1, time_m<=time2)]
ccn1_m = ccn1_m[np.logical_and(time_m>=time1, time_m<=time2)]
ccn2_m = ccn2_m[np.logical_and(time_m>=time1, time_m<=time2)]
ccn5_m = ccn5_m[np.logical_and(time_m>=time1, time_m<=time2)]
ncn3_m = ncn3_m[np.logical_and(time_m>=time1, time_m<=time2)]
ncn10_m = ncn10_m[np.logical_and(time_m>=time1, time_m<=time2)]
ncn100_m = ncn100_m[np.logical_and(time_m>=time1, time_m<=time2)]
cod_m = cod_m[np.logical_and(time_m>=time1, time_m<=time2)]
reff_m = reff_m[np.logical_and(time_m>=time1, time_m<=time2)]
lwp_m = lwp_m[np.logical_and(time_m>=time1, time_m<=time2)]
nd_m = nd_m[np.logical_and(time_m>=time1, time_m<=time2)]
precip_m = precip_m[np.logical_and(time_m>=time1, time_m<=time2)]
cld_m = cld_m[np.logical_and(time_m>=time1, time_m<=time2)]
lwdnsfc_m = lwdnsfc_m[np.logical_and(time_m>=time1, time_m<=time2)]
swdnsfc_m = swdnsfc_m[np.logical_and(time_m>=time1, time_m<=time2)]
swnetsfc_m = swnetsfc_m[np.logical_and(time_m>=time1, time_m<=time2)]
lwnetsfc_m = lwnetsfc_m[np.logical_and(time_m>=time1, time_m<=time2)]
swnettoa_m = swnettoa_m[np.logical_and(time_m>=time1, time_m<=time2)]
lwnettoa_m = lwnettoa_m[np.logical_and(time_m>=time1, time_m<=time2)]

ncn_m = ncn_m[np.logical_and(time_m>=time1, time_m<=time2)]

# E3SMv2 data
bc_m2 = bc_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
dst_m2 = dst_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ncl_m2 = ncl_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
so4_m2 = so4_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
org_m2 = org_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ccn1_m2 = ccn1_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ccn2_m2 = ccn2_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ccn5_m2 = ccn5_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ncn3_m2 = ncn3_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ncn10_m2 = ncn10_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
ncn100_m2 = ncn100_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
cod_m2 = cod_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
reff_m2 = reff_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
lwp_m2 = lwp_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
nd_m2 = nd_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
precip_m2 = precip_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
cld_m2 = cld_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
lwdnsfc_m2 = lwdnsfc_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
swdnsfc_m2 = swdnsfc_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
swnetsfc_m2 = swnetsfc_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
lwnetsfc_m2 = lwnetsfc_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
swnettoa_m2 = swnettoa_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]
lwnettoa_m2 = lwnettoa_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]

ncn_m2 = ncn_m2[np.logical_and(time_m2>=time1, time_m2<=time2)]


# unit change:
precip_m = precip_m*3600*1000   # m/s to mm/hr
precip_m2 = precip_m2*3600*1000   # m/s to mm/hr

# set a small threshold of E3SM precipitation
precip_m[precip_m<0.02] = 0
precip_m2[precip_m2<0.02] = 0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not os.path.exists(figpath):
    os.makedirs(figpath)

# #%% bar plot
# datagroup0 = [org,so4,nh4,no3,chl, [], []]
# datagroup1 = [org_m, so4_m, [], [], [], bc_m, dst_m]
# datagroup2 = [org_m2, so4_m2, [], [], [], bc_m2, dst_m2]
# dataall=[datagroup0,datagroup1, datagroup2,]
# labelall = ['Organic', 'SO$_4$', 'NH$_4$', 'NO$_3$', 'Chl', 'BC', 'Dust']
# colorall = ['limegreen', 'red', 'lightblue', 'orange', 'cyan', 'k', 'silver']
# fig,ax = plot.bar(dataall, datalabel=['Obs','E3SMv1','E3SMv2',], xlabel=None, ylabel='unit: $\mu$g/m$^3$', 
#                   title='Aerosol Composition  '+site+' '+IOP, varlabel= labelall, colorall=colorall)
# fig.savefig(figpath+'bar_composition_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% timeseries
# fig,ax = plot.timeseries([time,time,time], [org,org_m,org_m2], 
#                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                           title='Total Organic '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [so4,so4_m,so4_m2], 
#                           legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],  
#                           title='Sulfate '+site+' '+IOP, xlabel=None, ylabel='${\mu}$g/m$^{3}$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time], [ccn1,ccn1_m,ccn1_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='0.1%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_CCN1_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [ccn2,ccn2_m,ccn2_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='0.2%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [ccn5,ccn5_m,ccn5_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                          title='0.5%CCN '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_CCN5_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# if IOP=='IOP1':
#     fig,ax = plot.timeseries([time,time,time], [cpc3,ncn3_m,ncn3_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                              title='CN(>3nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
#     ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
#     ax.set_xlim(time1,time2)
#     fig.savefig(figpath+'timeseries_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     fig,ax = plot.timeseries([time,time,time], [cpc10,ncn10_m,ncn10_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                              title='CN(>10nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
#     ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
#     ax.set_xlim(time1,time2)
#     fig.savefig(figpath+'timeseries_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     fig,ax = plot.timeseries([time,time,time,time], [smps100,uhsas100,ncn100_m,ncn100_m2], 
#                             legend = ['SMPS','UHSAS','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                             title='CN(>100nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
#     ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
#     ax.set_xlim(time1,time2)
#     fig.savefig(figpath+'timeseries_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# elif IOP=='IOP2':
#     fig,ax = plot.timeseries([time,time,time], [smps100,ncn100_m,ncn100_m2], 
#                              legend = ['SMPS','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                             title='CN(>100nm) '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
#     ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
#     ax.set_xlim(time1,time2)
#     fig.savefig(figpath+'timeseries_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time,time], [cod,cod_sat,cod_m,cod_m2], marker='.',
#                          legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                         title='cloud optical depth '+site+' '+IOP, xlabel=None, ylabel=None)
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time,time,time,time], [lwp_mfrsr, lwp_sat, lwp_m, lwp_m2], 
                        legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                        title='LWP '+site+' '+IOP,xlabel=None, ylabel="g/m$^2$")
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
ax.set_xlim(time1,time2)
fig.savefig(figpath+'timeseries_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time], [LTS700,LTS850], legend = ['Theta(700hPa-sfc)','Theta(850hPa-sfc)'], 
#                         color=['b','r'], marker='.', figsize=(10,3), xlabel=None, ylabel="K")
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_LTS_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.timeseries([time,time,time,time], [ndrop,nd_sat, nd_m, nd_m2], marker='.',
                          legend = ['Ndrop (ARM)','Nd (satellite)','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                          title='Nd '+site+' '+IOP, xlabel=None, ylabel='cm$^{-3}$')
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
ax.set_xlim(time1,time2)
fig.savefig(figpath+'timeseries_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time,time], [reff,reff_sat,reff_m,reff_m2],  marker='.',
#                         legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                         title='Reff '+site+' '+IOP,xlabel=None, ylabel='$\mu$m')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time], [precip,precip_m,precip_m2],  
#                          legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='Precip '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='mm/hr')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time], [lwnetsfc,lwnetsfc_m,lwnetsfc_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='Sfc. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [swnetsfc,swnetsfc_m,swnetsfc_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='Sfc. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [lwnettoa,lwnettoa_m,lwnettoa_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='TOA. net LW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.timeseries([time,time,time], [swnettoa,swnettoa_m,swnettoa_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='TOA. net SW Flux '+site+' '+IOP, figsize=(10,3), xlabel=None, ylabel='W/m$^2$')
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.timeseries([time,time,time,time], [cld_arscl,cld_visst,cld_m,cld_m2], 
#                         legend = ['ARSCL','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                         title='Cloud fraction '+site+' '+IOP,xlabel=None, ylabel="%")
# ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# ax.set_xlim(time1,time2)
# fig.savefig(figpath+'timeseries_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# if IOP=='IOP1':
#     fig,ax = plot.timeseries_size([time,time,time,time], [size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
#                                   [uhsas_all.T.data, smps_all.T.data, ncn_m.T.data, ncn_m2.T.data], 
#                                   legend = ['UHSAS','SMPS','E3SMv1','E3SMv2'],
#                               ylabel='Diameter (nm)', ylimit=(3,1000),
#                               title = 'Aerosol Size Distribution (dN/dlogDp, cm$^{-3}$)')
# elif IOP=='IOP2':
#     fig,ax = plot.timeseries_size([time,time,time], [size_smps, np.arange(1,3001), np.arange(1,3001)], 
#                                   [smps_all.T.data, ncn_m.T.data, ncn_m2.T.data], legend = ['SMPS','E3SMv1','E3SMv2'],
#                               ylabel='Diameter (nm)', ylimit=(3,1000),
#                               title = 'Aerosol Size Distribution (dN/dlogDp, cm$^{-3}$)')    
# for ax_i in ax:
#     ax_i.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#     ax_i.xaxis.set_major_locator(mdates.DayLocator(interval=5))
# fig.savefig(figpath+'aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% diurnal cycle
# fig,ax = plot.diurnalcycle([org,org_m,org_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                           title='Organic '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
# fig.savefig(figpath+'diurnalcycle_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([so4,so4_m,so4_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                           title='Sulfate '+site+' '+IOP, xlabel='Time (UTC)', ylabel='${\mu}$g/m$^{3}$')
# fig.savefig(figpath+'diurnalcycle_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([ccn1,ccn1_m,ccn1_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='0.1%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
# fig.savefig(figpath+'diurnalcycle_CCN1_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([ccn2,ccn2_m,ccn2_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='0.2%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
# fig.savefig(figpath+'diurnalcycle_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([ccn5,ccn5_m,ccn5_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='0.5%CCN '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
# fig.savefig(figpath+'diurnalcycle_CCN5_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# if IOP=='IOP1':
#     fig,ax = plot.diurnalcycle([cpc3,ncn3_m,ncn3_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                             title='CN(>3nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
#     fig.savefig(figpath+'diurnalcycle_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     fig,ax = plot.diurnalcycle([cpc10,ncn10_m,ncn10_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                             title='CN(>10nm) '+site+' '+IOP, xlabel='Time (UTC)', ylabel='cm$^{-3}$')
#     fig.savefig(figpath+'diurnalcycle_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     fig,ax = plot.diurnalcycle([smps100,uhsas100,ncn100_m,ncn100_m2], legend = ['SMPS100','UHSAS100','E3SMv1','E3SMv2'], 
#                             title='CN(>100nm) '+site+' '+IOP, color=['k','gray','r','b'], xlabel='Time (UTC)',ylabel='cm$^{-3}$')
#     fig.savefig(figpath+'diurnalcycle_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# elif IOP=='IOP2':
#     fig,ax = plot.diurnalcycle([smps100,ncn100_m,ncn100_m2], legend = ['SMPS','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                             title='CN(>100nm) '+site+' '+IOP, xlabel='Time (UTC)',ylabel='cm$^{-3}$')
#     fig.savefig(figpath+'diurnalcycle_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle( [cod, cod_sat, cod_m, cod_m2], legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'], 
#                         title='Cloud optical depth '+site+' '+IOP, xlabel='Time (UTC)', ylabel=None)
# fig.savefig(figpath+'diurnalcycle_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([lwp_mfrsr,lwp_sat, lwp_m, lwp_m2], legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], 
#                         title='LWP '+site+' '+IOP, color=['k','gray','r','b'], xlabel='Time (UTC)',ylabel="g/m$^2$")
# fig.savefig(figpath+'diurnalcycle_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# # fig,ax = plot.diurnalcycle([LTS700,LTS850], legend = ['Theta(700hPa-sfc)','Theta(850hPa-sfc)'], 
# #                         color=['b','r'], xlabel='Time (UTC)', ylabel="K")
# # fig.savefig(figpath+'diurnalcycle_LTS_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([ndrop, nd_sat, nd_m,nd_m2], legend = ['Ndrop (ARM)', 'Nd (satellite)', 'E3SMv1', 'E3SMv2'], 
#                           title='Nd '+site+' '+IOP, color=['k','gray','r','b'], xlabel='Time (UTC)', ylabel='cm$^{-3}$')
# fig.savefig(figpath+'diurnalcycle_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([reff, reff_sat, reff_m, reff_m2], legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                         title='droplet effective radius '+site+' '+IOP, xlabel='Time (UTC)', ylabel='$\mu$m')
# fig.savefig(figpath+'diurnalcycle_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle( [precip,precip_m,precip_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         nozero_percentile=True, title='Precipitation '+site+' '+IOP, xlabel='Time (UTC)',ylabel='mm/hr')
# fig.savefig(figpath+'diurnalcycle_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([lwnetsfc,lwnetsfc_m,lwnetsfc_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='Sfc. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)',ylabel='W/m$^2$')
# fig.savefig(figpath+'diurnalcycle_LWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([swnetsfc,swnetsfc_m,swnetsfc_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='Sfc. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
# fig.savefig(figpath+'diurnalcycle_SWsfc_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([lwnettoa,lwnettoa_m,lwnettoa_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='TOA. net LW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
# fig.savefig(figpath+'diurnalcycle_LWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# fig,ax = plot.diurnalcycle([swnettoa,swnettoa_m,swnettoa_m2], legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                         title='TOA. net SW Flux '+site+' '+IOP, xlabel='Time (UTC)', ylabel='W/m$^2$')
# fig.savefig(figpath+'diurnalcycle_SWtoa_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# fig,ax = plot.diurnalcycle([cld_arscl,cld_visst,cld_m,cld_m2], legend = ['ARSCL','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                            title='Total cloud fraction '+site+' '+IOP, xlabel='Time (UTC)', ylabel="%")
# fig.savefig(figpath+'diurnalcycle_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# if IOP=='IOP1':
#     fig,ax = plot.diurnalcycle_2d([uhsas_all.T, smps_all.T, ncn_m.T, ncn_m2.T], 
#                                   y=[size_uhsas,size_smps, np.arange(1,3001), np.arange(1,3001)], 
#                                   title= ['UHSAS','SMPS','E3SMv1','E3SMv2'],
#                                   levellist=np.arange(0,11500,200), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
#                                   ylimit=(3,1000),cmap='jet')
# elif IOP=='IOP2':
#     fig,ax = plot.diurnalcycle_2d([smps_all.T, ncn_m.T, ncn_m2.T], y=[size_smps, np.arange(1,3001), np.arange(1,3001)], 
#                                   title= ['SMPS','E3SMv1','E3SMv2'],
#                                   levellist=np.arange(0,10000,200), xlabel='Time (UTC)', ylabel='Diameter (nm)', 
#                                   ylimit=(3,1000),cmap='jet')
# for ax_i in ax:
#     ax_i.set_yscale('log')
# fig.savefig(figpath+'diurnalcycle_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# #%% 1d histogram
# w0 = np.ones_like(org)/sum(~np.isnan(org.data))
# w1 = np.ones_like(org_m)/sum(~np.isnan(org_m.data))
# w2 = np.ones_like(org_m2)/sum(~np.isnan(org_m2.data))
# fig,ax = plot.hist([org,org_m,org_m2], weights=[w0,w1,w2],
#                         legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,10.2,0.5), 
#                           title='Total Organic '+site+' '+IOP,ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
# fig.savefig(figpath+'hist_org_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# w0 = np.ones_like(so4)/sum(~np.isnan(so4.data))
# w1 = np.ones_like(so4_m)/sum(~np.isnan(so4_m.data))
# w2 = np.ones_like(so4_m2)/sum(~np.isnan(so4_m2.data))
# fig,ax = plot.hist([so4,so4_m,so4_m2], weights=[w0,w1,w2],
#                         legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,8.2,0.5), 
#                           title='Sulfate '+site+' '+IOP, ylabel='Fraction', xlabel='${\mu}$g/m$^{3}$')
# fig.savefig(figpath+'hist_so4_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(ccn1)/sum(~np.isnan(ccn1.data))
# w1 = np.ones_like(ccn1_m)/sum(~np.isnan(ccn1_m.data))
# w2 = np.ones_like(ccn1_m2)/sum(~np.isnan(ccn1_m2.data))
# fig,ax = plot.hist([ccn1,ccn1_m,ccn1_m2], weights=[w0,w1,w2],
#                     legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,1100,50), 
#                     title='CCN (SS=0.1%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath+'hist_CCN1_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# w0 = np.ones_like(ccn2)/sum(~np.isnan(ccn2.data))
# w1 = np.ones_like(ccn2_m)/sum(~np.isnan(ccn2_m.data))
# w2 = np.ones_like(ccn2_m2)/sum(~np.isnan(ccn2_m2.data))
# fig,ax = plot.hist([ccn2,ccn2_m,ccn2_m2], weights=[w0,w1,w2],
#                     legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,1800,50), 
#                     title='CCN (SS=0.2%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath+'hist_CCN2_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# w0 = np.ones_like(ccn5)/sum(~np.isnan(ccn5.data))
# w1 = np.ones_like(ccn5_m)/sum(~np.isnan(ccn5_m.data))
# w2 = np.ones_like(ccn5_m2)/sum(~np.isnan(ccn5_m2.data))
# fig,ax = plot.hist([ccn5,ccn5_m,ccn5_m2], weights=[w0,w1,w2],
#                     legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], bins=np.arange(0,2500,100), 
#                     title='CCN (SS=0.5%) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
# fig.savefig(figpath+'hist_CCN5_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# if IOP=='IOP1':
#     w0 = np.ones_like(cpc3)/sum(~np.isnan(cpc3.data))
#     w1 = np.ones_like(ncn3_m)/sum(~np.isnan(ncn3_m.data))
#     w2 = np.ones_like(ncn3_m2)/sum(~np.isnan(ncn3_m2.data))
#     fig,ax = plot.hist([cpc3,ncn3_m,ncn3_m2], weights=[w0,w1,w2], bins=np.arange(0,29000,1000),
#                        legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'], 
#                             title='Aerosol number (>3nm) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
#     fig.savefig(figpath+'hist_CPC3_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     w0 = np.ones_like(cpc10)/sum(~np.isnan(cpc10.data))
#     w1 = np.ones_like(ncn10_m)/sum(~np.isnan(ncn10_m.data))
#     w2 = np.ones_like(ncn10_m2)/sum(~np.isnan(ncn10_m2.data))
#     fig,ax = plot.hist([cpc10,ncn10_m,ncn10_m2], weights=[w0,w1,w2], bins=np.arange(0,22000,1000),
#                        legend = ['Obs','E3SMv1','E3SMv2'], color=['k','r','b'],
#                         title='Aerosol number (>10nm) '+site+' '+IOP,ylabel='Fraction', xlabel='cm$^{-3}$')
#     fig.savefig(figpath+'hist_CPC10_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
#     w0 = np.ones_like(smps100)/sum(~np.isnan(smps100.data))
#     w00 = np.ones_like(uhsas100)/sum(~np.isnan(uhsas100.data))
#     w1 = np.ones_like(ncn100_m)/sum(~np.isnan(ncn100_m.data))
#     w2 = np.ones_like(ncn100_m2)/sum(~np.isnan(ncn100_m2.data))
#     fig,ax = plot.hist([smps100,uhsas100,ncn100_m,ncn100_m2], weights=[w0,w00,w1,w2],bins=np.arange(0,2100,100), 
#                        legend = ['SMPS100','UHSAS100','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                         title='Aerosol number (>100nm) '+site+' '+IOP, ylabel='Fraction', xlabel='cm$^{-3}$')
#     fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# elif IOP=='IOP2':
#     w0 = np.ones_like(smps100)/sum(~np.isnan(smps100.data))
#     w1 = np.ones_like(ncn100_m)/sum(~np.isnan(ncn100_m.data))
#     w2 = np.ones_like(ncn100_m2)/sum(~np.isnan(ncn100_m2.data))
#     fig,ax = plot.hist([smps100,ncn100_m,ncn100_m2], weights=[w0,w1,w2], bins=np.arange(0,2100,100),
#                        legend = ['SMPS100','E3SMv1','E3SMv2'], color=['k','r','b'],
#                         title='Aerosol number (>100nm) '+site+' '+IOP,  ylabel='Fraction', xlabel='cm$^{-3}$')
#     fig.savefig(figpath+'hist_CN100_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(cod)/sum(~np.isnan(cod.data))
# w00 = np.ones_like(cod_sat)/sum(~np.isnan(cod_sat.data))
# w1 = np.ones_like(cod_m)/sum(~np.isnan(cod_m.data))
# w2 = np.ones_like(cod_m2)/sum(~np.isnan(cod_m2.data))
# fig,ax = plot.hist( [cod, cod_sat, cod_m, cod_m2], weights=[w0,w00,w1,w2], 
#                     legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                     title='Cloud Optical Depth '+site+' '+IOP, bins=np.arange(0,90,5), ylabel='Fraction', xlabel=None)
# fig.savefig(figpath+'hist_cod_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w1 = np.ones_like(lwp_armbe)/sum(~np.isnan(lwp_armbe.data))
w0 = np.ones_like(lwp_mfrsr)/sum(~np.isnan(lwp_mfrsr.data))
w00 = np.ones_like(lwp_sat)/sum(~np.isnan(lwp_sat.data))
w1 = np.ones_like(lwp_m)/sum(~np.isnan(lwp_m.data))
w2 = np.ones_like(lwp_m2)/sum(~np.isnan(lwp_m2.data))
fig,ax = plot.hist([lwp_mfrsr, lwp_sat, lwp_m, lwp_m2], weights=[w0,w00,w1,w2], 
                    legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title='LWP '+site+' '+IOP, bins=np.arange(10,440,20), ylabel='Fraction', xlabel="g/m$^2$")
fig.savefig(figpath+'hist_LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w1 = np.ones_like(LTS700)/sum(~np.isnan(LTS700.data))
# w2 = np.ones_like(LTS850)/sum(~np.isnan(LTS850.data))
# fig,ax = plot.hist([LTS700,LTS850], weights=[w1,w2], legend = ['Theta(700hPa-sfc)','Theta(850hPa-sfc)'], 
#                     color=['b','r'], ylabel='Fraction', xlabel="K")
# fig.savefig(figpath+'hist_LTS_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

w0 = np.ones_like(ndrop)/sum(~np.isnan(ndrop.data))
w00 = np.ones_like(nd_sat)/sum(~np.isnan(nd_sat.data))
w1 = np.ones_like(nd_m)/sum(~np.isnan(nd_m.data))
w2 = np.ones_like(nd_m2)/sum(~np.isnan(nd_m2.data))
fig,ax = plot.hist([ndrop,nd_sat,nd_m,nd_m2],  weights=[w0,w00,w1,w2], 
                    legend = ['Ndrop','Satellite','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
                    title = 'Nd '+site+' '+IOP, bins=np.arange(0,510,20), ylabel='Fraction', xlabel='cm$^{-3}$')
fig.savefig(figpath+'hist_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(reff)/sum(~np.isnan(reff.data))
# w00 = np.ones_like(reff_sat)/sum(~np.isnan(reff_sat.data))
# w1 = np.ones_like(reff_m)/sum(~np.isnan(reff_m.data))
# w2 = np.ones_like(reff_m2)/sum(~np.isnan(reff_m2.data))
# fig,ax = plot.hist([reff,reff_sat,reff_m,reff_m2], weights=[w0,w00,w1,w2], 
#                    legend = ['MFRSR','VISST','E3SMv1','E3SMv2'], color=['k','gray','r','b'],
#                     title = 'Cloud Effective Radius '+site+' '+IOP, bins=np.arange(0,39,2), ylabel='Fraction', xlabel='$\mu$m')
# fig.savefig(figpath+'hist_reff_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# pr0 = precip[precip!=0]
# prm = precip_m[precip_m!=0]
# prm2 = precip_m[precip_m2!=0]
# w0 = np.ones_like(pr0)/sum(~np.isnan(pr0.data))
# w1 = np.ones_like(prm)/sum(~np.isnan(prm.data))
# w2 = np.ones_like(prm2)/sum(~np.isnan(prm2.data))
# fig,ax = plot.hist( [pr0,prm,prm2], weights=[w0,w1,w2], legend = ['Obs','E3SMv1','E3SMv2'], 
#                    color=['k','r','b'],  bins=np.arange(0,8,.2), 
#                     title = 'Precipitation '+site+' '+IOP, ylabel='Fraction', xlabel='mm/hr')
# fig.savefig(figpath+'hist_precip_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

# w0 = np.ones_like(cld_arscl)/sum(~np.isnan(cld_arscl.data))
# w00 = np.ones_like(cld_visst)/sum(~np.isnan(cld_visst.data))
# w1 = np.ones_like(cld_m)/sum(~np.isnan(cld_m.data))
# w2 = np.ones_like(cld_m2)/sum(~np.isnan(cld_m2.data))
# fig,ax = plot.hist([cld_arscl,cld_visst,cld_m,cld_m2], weights=[w0,w00,w1,w2], legend = ['ARSCL','VISST','E3SMv1','E3SMv2'], 
#                         color=['k','gray','r','b'], title = 'Cloud Fraction '+site+' '+IOP, ylabel='Fraction', xlabel="%")
# fig.savefig(figpath+'hist_totcld_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)


# #%% mean size distribution
# if IOP=='IOP1':
#     fig,ax = plot.mean_size([size_uhsas,size_smps,np.arange(1,3001),np.arange(1,3001)], 
#                 [np.nanmean(uhsas_all,axis=0), np.nanmean(smps_all,axis=0), np.nanmean(ncn_m,axis=0), np.nanmean(ncn_m2,axis=0)], 
#                 legend = ['UHSAS','SMPS','E3SMv1','E3SMv2'],color=['k','gray','r','b'], 
#                 marker=['o','+',None,None], linestyles=['none','none','-','-'],
#                 xlimit=(2, 2e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
#                 title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
#     fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)
# elif IOP=='IOP2':
#     fig,ax = plot.mean_size([size_smps,np.arange(1,3001),np.arange(1,3001)], 
#                             [np.nanmean(smps_all,axis=0), np.nanmean(ncn_m,axis=0), np.nanmean(ncn_m2,axis=0)], 
#                       legend = ['SMPS','E3SMv1','E3SMv2'],color=['k','r','b'],
#                       marker=['.',None,None], linestyles=['none','-','-'],
#                       xlimit=(1, 3e3), ylimit=(1e-2,1e4), xlabel='Diameter (nm)', ylabel='dN/dlogDp (cm$^{-3}$)', 
#                         title = 'Mean Aerosol Size Distribution '+site+' '+IOP)
#     fig.savefig(figpath+'mean_aerosol_size_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% scatter

fig,ax = plot.scatter([ndrop.data, nd_sat.data, nd_m.data, nd_m2.data], 
                      [lwp_mfrsr.data, lwp_sat.data, lwp_m.data, lwp_m2.data], 
                      title=['MFRSR','satellite','E3SMv1','E3SMv2'], xlimit=(0,400), ylimit=(0,500),
                    xlabel='Nd (cm$^{-3}$)', ylabel='LWP (g/m$^{2}$)', 
                    linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_LWP_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

fig,ax = plot.scatter([ndrop.data, nd_sat.data, nd_m.data, nd_m2.data], 
                      [reff.data*10, reff_sat.data*10, reff_m.data*10, reff_m2.data*10], 
                      title=['MFRSR','satellite','E3SMv1','E3SMv2'], xlimit=(0,400), ylimit=(0,300),
                    xlabel='Nd (cm$^{-3}$)', ylabel='Reff ($\mu$m *10)', 
                    linear_fit=True, intercept=True)
fig.savefig(figpath+'scatter_Reff_Nd_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% joint histogram
fig,ax = plot.jointhist([ccn2.data, ccn2.data, ccn2_m.data, ccn2_m2.data], [ndrop.data, nd_sat.data, nd_m.data, nd_m2.data], 
                    xedges=np.arange(0,1200,100), yedges=np.arange(0,500,50), 
                    normalize_x=False, vmin=0, vmax=0.2,
                    xlabel=['SFC 0.2%CCN (cm$^{-3}$)','SFC 0.2%CCN (cm$^{-3}$)','SFC 0.2%CCN (cm$^{-3}$)','SFC 0.2%CCN (cm$^{-3}$)'], 
                    ylabel='Nd (cm$^{-3}$)', title=['ARM','Satellite','E3SMv1','E3SMv2'], )
fig.savefig(figpath+'jointhist_Nd_CCN_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% heatmaps
fig,ax = plot.heatmap([ccn2.data, ccn2.data, ccn2_m.data, ccn2_m2.data], 
                      [lwp_mfrsr.data, lwp_sat.data, lwp_m.data, lwp_m2.data], 
                      [ndrop.data, nd_sat.data, nd_m.data, nd_m2.data], 
                    xedges=np.arange(50,1500,150), yedges=np.arange(0,500,50), 
                    vmin=0, vmax=200,
                    xlabel='0.2%CCN (cm$^{-3}$)', ylabel='LWP (g/m$^2$)', 
                    # title=['Ndrop_ARM (cm$^{-3}$)','Nd_Satellite (cm$^{-3}$)','Nd_E3SMv1 (cm$^{-3}$)','Nd_E3SMv2 (cm$^{-3}$)'])
                    title=['Ndrop_ARM','Nd_Satellite','Nd_E3SMv1','Nd_E3SMv2'])
fig.text(.94, .85,'cm$^{-3}$')      # unit of color variable
fig.savefig(figpath+'heatmap_Nd_CCN&LWP_'+site+'_'+IOP+'.png',dpi=fig.dpi,bbox_inches='tight', pad_inches=1)

#%% calculate statistics
# calc.mean_std_percentiles([org,org_m,org_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_ORG_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([so4, so4_m, so4_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_SO4_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn1,ccn1_m,ccn1_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN1_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn2,ccn2_m,ccn2_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN2_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ccn5,ccn5_m,ccn5_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                           outfile=figpath+'statistics_1var_CCN5_'+site+'_'+IOP+'.txt')
# if IOP=='IOP1':
#     calc.mean_std_percentiles([cpc3,ncn3_m,ncn3_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                               outfile=figpath+'statistics_1var_CPC3_'+site+'_'+IOP+'.txt')
#     calc.mean_std_percentiles([cpc10,ncn10_m,ncn10_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                               outfile=figpath+'statistics_1var_CPC10_'+site+'_'+IOP+'.txt')
#     calc.mean_std_percentiles([uhsas100, ncn100_m, ncn100_m2],legend=['Obs','E3SMv1','E3SMv2'], 
#                               outfile=figpath+'statistics_1var_UHSAS100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([smps100, ncn100_m, ncn100_m2],legend=['Obs','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_SMPS100_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cod,cod_sat, cod_m, cod_m2],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_COD_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([reff,reff_sat,reff_m,reff_m2],legend=['MFRSR','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Reff_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([lwp_mfrsr,lwp_armbe,lwp_sat,lwp_m,lwp_m2],legend=['MFRSR','ARMBE','Satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_LWP_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([ndrop,nd_sat,nd_m,nd_m2],legend=['Ndrop','Nd_satellite','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Nd_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([precip,precip_m,precip_m2],legend=['Obs','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_Precip_'+site+'_'+IOP+'.txt')
# calc.mean_std_percentiles([cld_arscl,cld_visst,cld_tsi,cld_m,cld_m2],legend=['ARSCL','VISST','TSI','E3SMv1','E3SMv2'],
#                           outfile=figpath+'statistics_1var_totcld_'+site+'_'+IOP+'.txt')


# calc.bias_corrcoef_RMSE(org,org_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_ORG_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(org,org_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_ORG_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(so4, so4_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_SO4_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(so4, so4_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_SO4_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn1,ccn1_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN1_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn1,ccn1_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN1_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn2,ccn2_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN2_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn2,ccn2_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN2_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ccn5,ccn5_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CCN5_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ccn5,ccn5_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CCN5_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')
# if IOP=='IOP1':
#     calc.bias_corrcoef_RMSE(cpc3,ncn3_m,label1='Obs',label2='E3SMv1', 
#                             outfile=figpath+'statistics_CN3nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
#     calc.bias_corrcoef_RMSE(cpc3,ncn3_m2,label1='Obs',label2='E3SMv2', 
#                             outfile=figpath+'statistics_CN3nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')
#     calc.bias_corrcoef_RMSE(cpc10,ncn10_m,label1='Obs',label2='E3SMv1', 
#                             outfile=figpath+'statistics_CN10nm_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
#     calc.bias_corrcoef_RMSE(cpc10,ncn10_m2,label1='Obs',label2='E3SMv2', 
#                             outfile=figpath+'statistics_CN10nm_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(smps100, ncn100_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_CN100_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(smps100, ncn100_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_CN100_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(lwp_mfrsr, lwp_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_lwp_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(lwp_mfrsr, lwp_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_lwp_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(ndrop, nd_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Nd_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(ndrop, nd_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Nd_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

# calc.bias_corrcoef_RMSE(reff, reff_m,label1='Obs',label2='E3SMv1', 
#                         outfile=figpath+'statistics_Reff_E3SMv1vsOBS_'+site+'_'+IOP+'.txt')
# calc.bias_corrcoef_RMSE(reff, reff_m2,label1='Obs',label2='E3SMv2', 
#                         outfile=figpath+'statistics_Reff_E3SMv2vsOBS_'+site+'_'+IOP+'.txt')

