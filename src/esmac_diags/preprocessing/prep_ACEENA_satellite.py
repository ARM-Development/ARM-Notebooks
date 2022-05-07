"""
prepare satellite data from ACEENA
options of output data into coarser resolution
"""

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
import time as ttt
import esmac_diags
from esmac_diags.subroutines.time_resolution_change import avg_time_1d, avg_time_2d
from esmac_diags.subroutines.time_format_change import datetime2cday
from esmac_diags.subroutines.specific_data_treatment import calc_cdnc_VISST, calc_clouddepth_VISST, insolation

# visstgridpath = '../../../data/ACEENA/obs/satellite/visst/grid/'
# visstpixpath = '../../../data/ACEENA/obs/satellite/visst/pix_3x3/'
# predatapath = 'C:/Users/tang357/Downloads/ACEENA/'
# dt=3600


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_VISST_grid(visstgridpath, predatapath, dt=3600):
    """
    prepare VISST-satellite data in grid level (0.5x0.5 degrees)

    Parameters
    ----------
    visstgridpath : str
        input datapath
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                               
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% settings
    
    # index of lon(x), lat(y) at the site
    # site=='ENA':
    x_idx = 9
    y_idx = 7
    
    #%% read in data
    lst = glob.glob(os.path.join(visstgridpath, '*visstgrid*.cdf'))
    filetime = [a.split('.c1.')[1] for a in lst]
    sortidx = np.argsort(filetime)
    # first data
    visstdata = xr.open_dataset(lst[sortidx[0]])
    vissttime = visstdata['time']
    lat = visstdata['latitude'][y_idx]
    lon = visstdata['longitude'][x_idx]
    # check in case the index is incorrect
    if np.abs(lat-39.09527)>0.5 or np.abs(lon+28.0339)>0.5:
        print(lat, lon)
        raise ValueError('index at ENA may not right, check x_idx and y_idx')
            
    solar_zenith = visstdata['solar_zenith_angle'][:,y_idx,x_idx]
    clearsky_vis_reflectance = visstdata['clearsky_vis_reflectance'][:,y_idx,x_idx]
    vis_reflectance_all = visstdata['visible_reflectance'][:,y_idx,x_idx,0]
    vis_reflectance_clr = visstdata['visible_reflectance'][:,y_idx,x_idx,1]
    lwp = visstdata['water_path'][:,y_idx,x_idx,1]
    iwp = visstdata['water_path'][:,y_idx,x_idx,0]
    sfc_net_sw = visstdata['surface_net_shortwave_flux'][:,y_idx,x_idx]
    sfc_net_lw = visstdata['surface_net_longwave_flux'][:,y_idx,x_idx]
    sfc_down_sw = visstdata['surface_down_shortwave_flux'][:,y_idx,x_idx]
    sfc_down_lw = visstdata['surface_down_longwave_flux'][:,y_idx,x_idx]
    reff_liq = visstdata['particle_size'][:,y_idx,x_idx,1]
    cod_liq_linavg = visstdata['optical_depth_linear'][:,y_idx,x_idx,2]
    cod_liq_logavg = visstdata['optical_depth_log'][:,y_idx,x_idx,2]
    ctt_liq = visstdata['cloud_temperature'][:,y_idx,x_idx,2]
    ctp_liq = visstdata['cloud_pressure_top'][:,y_idx,x_idx,2]
    cth_liq = visstdata['cloud_height_top'][:,y_idx,x_idx,2]
    cf_all = visstdata['cloud_percentage'][:,y_idx,x_idx,0]
    cf_liq = visstdata['cloud_percentage'][:,y_idx,x_idx,2]
    cf_allz = visstdata['cloud_percentage_level'][:,y_idx,x_idx,0]
    cf_low = visstdata['cloud_percentage_level'][:,y_idx,x_idx,1]
    cf_mid = visstdata['cloud_percentage_level'][:,y_idx,x_idx,2]
    cf_high = visstdata['cloud_percentage_level'][:,y_idx,x_idx,3]
    bb_lw_all = visstdata['broadband_longwave_flux'][:,y_idx,x_idx,0]
    bb_sw_albedo_all = visstdata['broadband_shortwave_albedo'][:,y_idx,x_idx,0]
    bb_lw_clr = visstdata['broadband_longwave_flux'][:,y_idx,x_idx,1]
    bb_sw_albedo_clr = visstdata['broadband_shortwave_albedo'][:,y_idx,x_idx,1]
    visstdata.close()
    for ii in range(1,len(lst)):
        file = lst[sortidx[ii]]
        print(file)
        visstdata = xr.open_dataset(file)
        vissttime = xr.concat([vissttime, visstdata['time']], dim="time")
        solar_zenith = xr.concat([solar_zenith, visstdata['solar_zenith_angle'][:,y_idx,x_idx]], dim="time")
        clearsky_vis_reflectance = xr.concat([clearsky_vis_reflectance, visstdata['clearsky_vis_reflectance'][:,y_idx,x_idx]], dim="time")
        vis_reflectance_all = xr.concat([vis_reflectance_all, visstdata['visible_reflectance'][:,y_idx,x_idx,0]], dim="time")
        vis_reflectance_clr = xr.concat([vis_reflectance_clr, visstdata['visible_reflectance'][:,y_idx,x_idx,1]], dim="time")
        lwp = xr.concat([lwp, visstdata['water_path'][:,y_idx,x_idx,1]], dim="time")
        iwp = xr.concat([iwp, visstdata['water_path'][:,y_idx,x_idx,0]], dim="time")
        sfc_net_sw = xr.concat([sfc_net_sw, visstdata['surface_net_shortwave_flux'][:,y_idx,x_idx]], dim="time")
        sfc_net_lw = xr.concat([sfc_net_lw, visstdata['surface_net_longwave_flux'][:,y_idx,x_idx]], dim="time")
        sfc_down_sw = xr.concat([sfc_down_sw, visstdata['surface_down_shortwave_flux'][:,y_idx,x_idx]], dim="time")
        sfc_down_lw = xr.concat([sfc_down_lw, visstdata['surface_down_longwave_flux'][:,y_idx,x_idx]], dim="time")
        reff_liq = xr.concat([reff_liq, visstdata['particle_size'][:,y_idx,x_idx,1]], dim="time")
        cod_liq_linavg = xr.concat([cod_liq_linavg, visstdata['optical_depth_linear'][:,y_idx,x_idx,2]], dim="time")
        cod_liq_logavg = xr.concat([cod_liq_logavg, visstdata['optical_depth_log'][:,y_idx,x_idx,2]], dim="time")
        ctt_liq = xr.concat([ctt_liq, visstdata['cloud_temperature'][:,y_idx,x_idx,2]], dim="time")
        ctp_liq = xr.concat([ctp_liq, visstdata['cloud_pressure_top'][:,y_idx,x_idx,2]], dim="time")
        cth_liq = xr.concat([cth_liq, visstdata['cloud_height_top'][:,y_idx,x_idx,2]], dim="time")
        cf_all = xr.concat([cf_all, visstdata['cloud_percentage'][:,y_idx,x_idx,0]], dim="time")
        cf_liq = xr.concat([cf_liq, visstdata['cloud_percentage'][:,y_idx,x_idx,2]], dim="time")
        cf_allz = xr.concat([cf_allz, visstdata['cloud_percentage_level'][:,y_idx,x_idx,0]], dim="time")
        cf_low = xr.concat([cf_low, visstdata['cloud_percentage_level'][:,y_idx,x_idx,1]], dim="time")
        cf_mid = xr.concat([cf_mid, visstdata['cloud_percentage_level'][:,y_idx,x_idx,2]], dim="time")
        cf_high = xr.concat([cf_high, visstdata['cloud_percentage_level'][:,y_idx,x_idx,3]], dim="time")
        bb_lw_all = xr.concat([bb_lw_all, visstdata['broadband_longwave_flux'][:,y_idx,x_idx,0]], dim="time")
        bb_lw_clr = xr.concat([bb_lw_clr, visstdata['broadband_longwave_flux'][:,y_idx,x_idx,1]], dim="time")
        bb_sw_albedo_all = xr.concat([bb_sw_albedo_all, visstdata['broadband_shortwave_albedo'][:,y_idx,x_idx,0]], dim="time")
        bb_sw_albedo_clr = xr.concat([bb_sw_albedo_clr, visstdata['broadband_shortwave_albedo'][:,y_idx,x_idx,1]], dim="time")
        visstdata.close()
        
    #%% calculate TOA SW flux from albedo
    
    # change time to calendar day
    calday = datetime2cday(vissttime.data)
    # calculate insolation
    ins = insolation(calday, lon.data, lat.data, leap_year='leap')
    
    # calculate net SW flux
    bb_sw_all = ins * (1 - bb_sw_albedo_all*0.01)
    
    #%% retrieve CDNC
    lwp = lwp.data
    ctt = ctt_liq.data
    cod = cod_liq_linavg.data
    H = calc_clouddepth_VISST(lwp, ctt, adiabaticity=0.8)
    H_ad = calc_clouddepth_VISST(lwp, ctt, adiabaticity=1.0)
    Nd = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=0.8)
    Nd_ad = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=1.0)
    
    #filter out columns with ice and bad retrievals
    H_array = np.array(H)
    H_ad_array = np.array(H_ad)
    Nd_array = np.array(Nd)
    Nd_ad_array = np.array(Nd_ad)
    
    ind = np.array(iwp > 0)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    ind = np.isinf(Nd_array)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(vissttime[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(vissttime[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    Nd_new = avg_time_1d(vissttime, Nd_array, time_new)
    H_new = avg_time_1d(vissttime, H, time_new)
    lwp_new = avg_time_1d(vissttime, lwp, time_new)
    iwp_new = avg_time_1d(vissttime, iwp, time_new)
    swnetsfc_new = avg_time_1d(vissttime, sfc_net_sw, time_new)
    lwnetsfc_new = avg_time_1d(vissttime, sfc_net_lw, time_new)
    swdnsfc_new = avg_time_1d(vissttime, sfc_down_sw, time_new)
    lwdnsfc_new = avg_time_1d(vissttime, sfc_down_lw, time_new)
    reff_new = avg_time_1d(vissttime, reff_liq, time_new)
    cod_new = avg_time_1d(vissttime, cod_liq_linavg, time_new)
    codlog_new = avg_time_1d(vissttime, cod_liq_logavg, time_new)
    cf_all_new = avg_time_1d(vissttime, cf_allz, time_new)
    cf_low_new = avg_time_1d(vissttime, cf_low, time_new)
    cf_mid_new = avg_time_1d(vissttime, cf_mid, time_new)
    cf_high_new = avg_time_1d(vissttime, cf_high, time_new)
    ctt_new = avg_time_1d(vissttime, ctt_liq, time_new)
    ctp_new = avg_time_1d(vissttime, ctp_liq, time_new)
    cth_new = avg_time_1d(vissttime, cth_liq, time_new)
    lw_new = avg_time_1d(vissttime, bb_lw_all, time_new)
    sw_new = avg_time_1d(vissttime, bb_sw_all, time_new)
    
    
    #%% output file
    outfile = predatapath + 'Nd_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Nd': (['time'], np.float32(Nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Nd'].attrs["long_name"] = 'cloud droplet number concentration'
    ds['Nd'].attrs["units"] = '#/cm3'
    
    ds.attrs["title"] = 'cloud droplet number concentration retrieved from VISST 0.5x0.5 data'
    ds.attrs["description"] = 'retrieved following Bennartz 2007, assuming adiabaticity = 0.8'
    ds.attrs["reference"] = 'https://doi.org/10.1029/2006JD007547'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Hcld_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Hcld': (['time'], np.float32(H_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Hcld'].attrs["long_name"] = 'cloud depth for liquid cloud only'
    ds['Hcld'].attrs["units"] = 'm'
    
    ds.attrs["title"] = 'liquid cloud depth retrieved from VISST 0.5x0.5 data'
    ds.attrs["description"] = 'retrieved from LWP and cloud top temperature assuming adiabaticity = 0.8'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'LWP_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwp': (['time'], np.float32(lwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwp'].attrs["long_name"] = 'liquid water path'
    ds['lwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'liquid water path from VISST 0.5x0.5 data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'IWP_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'iwp': (['time'], np.float32(iwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['iwp'].attrs["long_name"] = 'ice water path'
    ds['iwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'ice water path from VISST 0.5x0.5 data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Reff_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'reff': (['time'], np.float32(reff_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['reff'].attrs["long_name"] = 'effective radius for liquid clouds'
    ds['reff'].attrs["units"] = 'um'
    
    ds.attrs["title"] = 'liquid clouds effective radius from VISST 0.5x0.5 data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cod_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cod': (['time'], np.float32(cod_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cod'].attrs["long_name"] = 'cloud optical depth for liquid clouds'
    ds['cod'].attrs["units"] = 'N/A'
    
    ds.attrs["title"] = 'liquid clouds opical depth from VISST 0.5x0.5 data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cloudfraction_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cldtot': (['time'], np.float32(cf_all_new)),
                    'cldhigh': (['time'], np.float32(cf_high_new)),
                    'cldmid': (['time'], np.float32(cf_mid_new)),
                    'cldlow': (['time'], np.float32(cf_low_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cldtot'].attrs["long_name"] = 'cloud fraction for all heights'
    ds['cldtot'].attrs["units"] = '%'
    ds['cldhigh'].attrs["long_name"] = 'cloud fraction for high clouds (>6km)'
    ds['cldhigh'].attrs["units"] = '%'
    ds['cldmid'].attrs["long_name"] = 'cloud fraction for middle clouds (2-6km)'
    ds['cldmid'].attrs["units"] = '%'
    ds['cldlow'].attrs["long_name"] = 'cloud fraction for low clouds (0-2km)'
    ds['cldlow'].attrs["units"] = '%'
    
    ds.attrs["title"] = 'cloud fraction from VISST 0.5x0.5 data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cloudtop_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'ctt': (['time'], np.float32(ctt_new)),
                    'cth': (['time'], np.float32(cth_new)),
                    'ctp': (['time'], np.float32(ctp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['ctt'].attrs["long_name"] = 'cloud top temperature for liquid clouds'
    ds['ctt'].attrs["units"] = 'K'
    ds['ctp'].attrs["long_name"] = 'cloud top pressure for liquid clouds'
    ds['ctp'].attrs["units"] = 'hPa'
    ds['cth'].attrs["long_name"] = 'cloud top height for liquid clouds'
    ds['cth'].attrs["units"] = 'km'
    
    ds.attrs["title"] = 'cloud top temperature, pressure and height from VISST 0.5x0.5 data'
    ds.attrs["description"] = 'for liquid clouds only'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'lwflx_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwnettoa': (['time'], np.float32(lw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwnettoa'].attrs["long_name"] = 'net LW flux at TOA'
    ds['lwnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net longwave flux at TOA from VISST 0.5x0.5 data'
    ds.attrs["description"] = 'upward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'swflx_VISSTgrid_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'swnettoa': (['time'], np.float32(sw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['swnettoa'].attrs["long_name"] = 'net SW flux at TOA'
    ds['swnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net shortwave flux at TOA from VISST 0.5x0.5 data'
    ds.attrs["description"] = 'calculate from insolation and SW albedo, downward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def prep_VISST_pixel(visstpixpath, predatapath, dt=3600):
    """
    prepare VISST-satellite data in pixel level (4km)

    Parameters
    ----------
    visstpixpath : str
        input datapath
    predatapath : str
        output datapath
    dt : float
        time resolution (unit: sec) of output

    Returns
    -------
    None.

    """
                               
    if not os.path.exists(predatapath):
        os.makedirs(predatapath)
    
    #%% settings
    
    # read data from extracted pixel-level data in 3x3 grid centered at SGP
    x_idx = 1
    y_idx = 1
    
    #%% read in data
    lst = glob.glob(os.path.join(visstpixpath, '*visstpx2d*.cdf'))
    # lst = glob.glob(os.path.join(visstpixpath, 'enavisstpx2*.c1.201802*.cdf'))
    filetime = [a.split('.c1.')[1] for a in lst]
    sortidx = np.argsort(filetime)
    # first data
    visstdata = xr.open_dataset(lst[sortidx[0]])
    vissttime = visstdata['time_offset']
    lat = visstdata['latitude'][x_idx, y_idx]
    lon = visstdata['longitude'][x_idx, y_idx]
    # check in case the index is incorrect
    if np.abs(lat-39.09527)>0.5 or np.abs(lon+28.0339)>0.5:
        print(lat, lon)
        raise ValueError('index at ENA may not right, check x_idx and y_idx')
    vis_reflectance = visstdata['reflectance_vis'][x_idx, y_idx]
    wp = visstdata['cloud_lwp_iwp'][x_idx, y_idx]
    phase = visstdata['cloud_phase'][x_idx, y_idx] #0:snow, 1:water, 2:ice, 3:no retrieval, 4:clear, 5:bad data, 6:suspected water, 7:suspected ice, 13:cleaned data
    particle_size = visstdata['cloud_particle_size'][x_idx, y_idx]
    cod = visstdata['cloud_visible_optical_depth'][x_idx, y_idx]
    ctt = visstdata['cloud_top_temperature'][x_idx, y_idx]
    ctp = visstdata['cloud_top_pressure'][x_idx, y_idx]
    cth = visstdata['cloud_top_height'][x_idx, y_idx]
    bb_lw = visstdata['broadband_longwave_flux'][x_idx, y_idx]
    bb_sw_albedo = visstdata['broadband_shortwave_albedo'][x_idx, y_idx]
    visstdata.close()
    for file in lst[1:]:
        print(file)
        visstdata = xr.open_dataset(file)
        lat = visstdata['latitude'][x_idx, y_idx]
        lon = visstdata['longitude'][x_idx, y_idx]
        # check each file since grid information in pixel-level data may change
        if np.abs(lat-39.09527)>0.5 or np.abs(lon+28.0339)>0.5:
            print(lat, lon)
            raise ValueError('index at ENA may not right, check x_idx and y_idx')
        vissttime = xr.concat([vissttime, visstdata['time_offset']], dim="time")
        vis_reflectance = xr.concat([vis_reflectance, visstdata['reflectance_vis'][x_idx, y_idx]], dim="time")
        wp = xr.concat([wp, visstdata['cloud_lwp_iwp'][x_idx, y_idx]], dim="time")
        phase = xr.concat([phase, visstdata['cloud_phase'][x_idx, y_idx]], dim="time")
        particle_size = xr.concat([particle_size, visstdata['cloud_particle_size'][x_idx, y_idx]], dim="time")
        cod = xr.concat([cod, visstdata['cloud_visible_optical_depth'][x_idx, y_idx]], dim="time")
        ctt = xr.concat([ctt, visstdata['cloud_top_temperature'][x_idx, y_idx]], dim="time")
        ctp = xr.concat([ctp, visstdata['cloud_top_pressure'][x_idx, y_idx]], dim="time")
        cth = xr.concat([cth, visstdata['cloud_top_height'][x_idx, y_idx]], dim="time")
        bb_lw = xr.concat([bb_lw, visstdata['broadband_longwave_flux'][x_idx, y_idx]], dim="time")
        bb_sw_albedo = xr.concat([bb_sw_albedo, visstdata['broadband_shortwave_albedo'][x_idx, y_idx]], dim="time")
        visstdata.close()
        
    #%% calculate TOA SW flux from albedo
    
    # change time to calendar day
    calday = datetime2cday(vissttime.data)
    # calculate insolation
    ins = insolation(calday, lon.data, lat.data, leap_year='leap')
    
    # calculate net SW flux
    bb_sw = ins * (1 - bb_sw_albedo*0.01)
    
    #%% retrieve CDNC
    lwp = wp.data
    ctt = ctt.data
    cod = cod.data
    H = calc_clouddepth_VISST(lwp, ctt, adiabaticity=0.8)
    H_ad = calc_clouddepth_VISST(lwp, ctt, adiabaticity=1.0)
    Nd = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=0.8)
    Nd_ad = calc_cdnc_VISST(lwp, ctt, cod, adiabaticity=1.0)
    
    #filter out columns with ice and bad retrievals
    H_array = np.array(H)
    H_ad_array = np.array(H_ad)
    Nd_array = np.array(Nd)
    Nd_ad_array = np.array(Nd_ad)
    
    ind = np.isinf(Nd_array)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    #%% only choose liquid water clouds. cloud phase value:
    #       cloud_phase:value_0 = "snow" ;
    # 		cloud_phase:value_1 = "water" ;
    # 		cloud_phase:value_2 = "ice" ;
    # 		cloud_phase:value_3 = "no retrieval" ;
    # 		cloud_phase:value_4 = "clear" ;
    # 		cloud_phase:value_5 = "bad data" ;
    # 		cloud_phase:value_6 = "suspected water" ;
    # 		cloud_phase:value_7 = "suspected ice" ;
    # 		cloud_phase:value_13 = "cleaned data" ;
    
    ind = np.array(phase != 1)
    H_array[ind] = np.nan
    H_ad_array[ind] = np.nan
    Nd_array[ind] = np.nan
    Nd_ad_array[ind] = np.nan
    
    # effective radius
    reff = particle_size.data
    reff[ind] = np.nan
    
    # lwp and iwp
    lwp = np.array(wp.data)
    iwp = np.array(wp.data)
    lwp[phase!=1] = np.nan
    iwp[phase!=2] = np.nan
    
    #%% re-shape the data into coarser resolution
    startdate = np.datetime_as_string(np.datetime64(vissttime[0].data))[:10]
    enddate = np.datetime_as_string(np.datetime64(vissttime[-1].data))[:10]
    
    time_new = pd.date_range(start=startdate, end=enddate, freq=str(int(dt))+"s")
    
    Nd_new = avg_time_1d(vissttime, Nd_array, time_new)
    H_new = avg_time_1d(vissttime, H, time_new)
    lwp_new = avg_time_1d(vissttime, lwp, time_new)
    iwp_new = avg_time_1d(vissttime, iwp, time_new)
    reff_new = avg_time_1d(vissttime, reff, time_new)
    cod_new = avg_time_1d(vissttime, cod, time_new)
    ctt_new = avg_time_1d(vissttime, ctt, time_new)
    ctp_new = avg_time_1d(vissttime, ctp, time_new)
    cth_new = avg_time_1d(vissttime, cth, time_new)
    lw_new = avg_time_1d(vissttime, bb_lw, time_new)
    sw_new = avg_time_1d(vissttime, bb_sw, time_new)
    
    #%% output file
    outfile = predatapath + 'Nd_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Nd': (['time'], np.float32(Nd_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Nd'].attrs["long_name"] = 'cloud droplet number concentration'
    ds['Nd'].attrs["units"] = '#/cm3'
    
    ds.attrs["title"] = 'cloud droplet number concentration retrieved from VISST 4x4km data'
    ds.attrs["description"] = 'retrieved following Bennartz 2007, assuming adiabaticity = 0.8'
    ds.attrs["reference"] = 'https://doi.org/10.1029/2006JD007547'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Hcld_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'Hcld': (['time'], np.float32(H_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['Hcld'].attrs["long_name"] = 'cloud depth for liquid cloud only'
    ds['Hcld'].attrs["units"] = 'm'
    
    ds.attrs["title"] = 'liquid cloud depth retrieved from VISST 4x4km data'
    ds.attrs["description"] = 'retrieved from LWP and cloud top temperature assuming adiabaticity = 0.8'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'LWP_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwp': (['time'], np.float32(lwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwp'].attrs["long_name"] = 'liquid water path'
    ds['lwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'liquid water path from VISST 4x4km data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'IWP_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'iwp': (['time'], np.float32(iwp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['iwp'].attrs["long_name"] = 'ice water path'
    ds['iwp'].attrs["units"] = 'g/m2'
    
    ds.attrs["title"] = 'ice water path from VISST 4x4km data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'Reff_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'reff': (['time'], np.float32(reff_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['reff'].attrs["long_name"] = 'effective radius for liquid clouds'
    ds['reff'].attrs["units"] = 'um'
    
    ds.attrs["title"] = 'liquid clouds effective radius from VISST 4x4km data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cod_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'cod': (['time'], np.float32(cod_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['cod'].attrs["long_name"] = 'cloud optical depth for liquid clouds'
    ds['cod'].attrs["units"] = 'N/A'
    
    ds.attrs["title"] = 'liquid clouds opical depth from VISST 4x4km data'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'cloudtop_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'ctt': (['time'], np.float32(ctt_new)),
                    'cth': (['time'], np.float32(cth_new)),
                    'ctp': (['time'], np.float32(ctp_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['ctt'].attrs["long_name"] = 'cloud top temperature for liquid clouds'
    ds['ctt'].attrs["units"] = 'K'
    ds['ctp'].attrs["long_name"] = 'cloud top pressure for liquid clouds'
    ds['ctp'].attrs["units"] = 'hPa'
    ds['cth'].attrs["long_name"] = 'cloud top height for liquid clouds'
    ds['cth'].attrs["units"] = 'km'
    
    ds.attrs["title"] = 'cloud top temperature, pressure and height from VISST 4x4km data'
    ds.attrs["description"] = 'for any cloud type'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'lwflx_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'lwnettoa': (['time'], np.float32(lw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['lwnettoa'].attrs["long_name"] = 'net LW flux at TOA'
    ds['lwnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net longwave flux at TOA from VISST 4x4km data'
    ds.attrs["description"] = 'upward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
    #
    outfile = predatapath + 'swflx_VISSTpix_ACEENA.nc'
    print('output file '+outfile)
    ds = xr.Dataset({
                    'swnettoa': (['time'], np.float32(sw_new)),
                    },
                     coords={'time': ('time', time_new)})
    #assign attributes
    ds['time'].attrs["long_name"] = "Time"
    ds['time'].attrs["standard_name"] = "time"
    ds['swnettoa'].attrs["long_name"] = 'net SW flux at TOA'
    ds['swnettoa'].attrs["units"] = 'W/m2'
    
    ds.attrs["title"] = 'net shortwave flux at TOA from VISST 4x4km data'
    ds.attrs["description"] = 'calculate from insolation and SW albedo, downward positive'
    ds.attrs["date"] = ttt.ctime(ttt.time())
    
    ds.to_netcdf(outfile, mode='w')
    
