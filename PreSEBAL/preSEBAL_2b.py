# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 19:15:02 2018

@author: tih
"""
def main(VegetationExcel, rasters_path_METEO, name_format_METEO, cores = 1):

    # import general modules
    from openpyxl import load_workbook
    import os
    import glob
    import gdal
    import pandas as pd
    
    # import modules SEBAL
    from SEBAL.hants.wa_gdal import main as HANTS  # from hants.wa_arcpy import *
    
    # Open Excelfile
    wb_veg = load_workbook(VegetationExcel, data_only=True)
    
    # Open worksheets
    ws_veg = wb_veg['General_Input']
    ws_inputs = wb_veg['HANTS_Input']
    
    # input parameters
    input_folder = r"%s" %str(ws_veg['B7'].value)
    start_date = "%s" %str(ws_veg['B2'].value)
    end_date = "%s" %str(ws_veg['B3'].value)
    Dates = pd.date_range(start_date, end_date, freq = "D")

    # Create input datasets for HANTS
    rasters_path_LST = os.path.join(input_folder, "PreSEBAL_SEBAL_out","Surface_Temperature")
    name_format_LST = 'Surface_Temperature_VIIRS_{yyyy}{mm}{dd}.tif'
    
    # Find extend LST
    os.chdir(rasters_path_LST)
    file_example_LST = glob.glob("*.tif")[0]
    filename_in_ex_LST = os.path.join(rasters_path_LST, file_example_LST)
    
    dest_LST = gdal.Open(filename_in_ex_LST)
    xsize_LST = dest_LST.RasterXSize
    ysize_LST = dest_LST.RasterYSize
    geo_LST = dest_LST.GetGeoTransform()
    
    latlim_LST = [geo_LST[3] + geo_LST[5] * ysize_LST, geo_LST[3]]
    lonlim_LST = [geo_LST[0], geo_LST[0] + geo_LST[1] * xsize_LST]
    
    cellsize_LST = geo_LST[1]

    # Find extend Meteo
    os.chdir(rasters_path_METEO)
    file_example_METEO = glob.glob("*.tif")[0]
    filename_in_ex_METEO = os.path.join(rasters_path_METEO, file_example_METEO)
    
    dest_METEO = gdal.Open(filename_in_ex_METEO)
    xsize_METEO = dest_METEO.RasterXSize
    ysize_METEO = dest_METEO.RasterYSize
    geo_METEO = dest_METEO.GetGeoTransform()
    
    latlim_METEO = [geo_METEO[3] + geo_METEO[5] * ysize_METEO, geo_METEO[3]]
    lonlim_METEO= [geo_METEO[0], geo_METEO[0] + geo_METEO[1] * xsize_METEO]
    
    cellsize_METEO = geo_METEO[1]
    
    # HANTS general paramters
    nb = len(Dates)
    
    # HANTS parameters LST
    nf_lst = int(ws_inputs['D2'].value)   # number of frequencies to be considered above the zero frequency
    low_lst = float(ws_inputs['D3'].value) # valid range minimum
    high_lst = float(ws_inputs['D4'].value) # valid range maximum
    HiLo_lst = str(ws_inputs['D5'].value) # 2-character string indicating rejection of high or low outliers
    fet_lst = float(ws_inputs['D6'].value)  # fit error tolerance (point eviating more than fet from curve fit are rejected)
    delta_lst = float(ws_inputs['D7'].value) # small positive number e.g. 0.1 to supress high amplitudes
    dod_lst = float(ws_inputs['D8'].value) # degree of overdeterminedness (iteration stops if number of points reaches the minimum required for curve fitting, plus dod). This is a safety measure

    # HANTS parameters Meteo 
    nf_meteo = int(ws_inputs['E2'].value)   # number of frequencies to be considered above the zero frequency
    low_meteo = float(ws_inputs['E3'].value) # valid range minimum
    high_meteo = float(ws_inputs['E4'].value) # valid range maximum
    HiLo_meteo = str(ws_inputs['E5'].value) # 2-character string indicating rejection of high or low outliers
    fet_meteo = float(ws_inputs['E6'].value)  # fit error tolerance (point eviating more than fet from curve fit are rejected)
    delta_meteo = float(ws_inputs['E7'].value) # small positive number e.g. 0.1 to supress high amplitudes
    dod_meteo = float(ws_inputs['E8'].value) # degree of overdeterminedness (iteration stops if number of points reaches the minimum required for curve fitting, plus dod). This is a safety measure

    # Define outputs HANTS    
    nc_path_lst = os.path.join(input_folder, "PreSEBAL_SEBAL_out", "LST_HANTS_nf%s_l%s_h%s_%s_fet%s_delta%s_dod%s.nc" %(nf_lst, low_lst, high_lst, HiLo_lst, fet_lst, delta_lst, dod_lst))  
    nc_path_meteo = os.path.join(input_folder, "PreSEBAL_SEBAL_out", "METEO_HANTS_nf%s_l%s_h%s_%s_fet%s_delta%s_dod%s.nc" %(nf_meteo, low_meteo, high_meteo, HiLo_meteo, fet_meteo, delta_meteo, dod_meteo))    

    import time
    
    time_now = time.time()
    # Run
    HANTS.run_HANTS(rasters_path_LST, name_format_LST,
              start_date, end_date, latlim_LST, lonlim_LST, cellsize_LST, nc_path_lst,
              nb, nf_lst, HiLo_lst, low_lst, high_lst, fet_lst, dod_lst, delta_lst, 0.001,
              4326, cores)
    
    time_end = time.time()
    print('HANTS LST executed in ', time_end - time_now, 'Seconds')

    time_now = time.time()
    # Run
    HANTS.run_HANTS(rasters_path_METEO, name_format_METEO,
              start_date, end_date, latlim_METEO, lonlim_METEO, cellsize_METEO, nc_path_meteo,
              nb, nf_meteo, HiLo_meteo, low_meteo, high_meteo, fet_meteo, dod_meteo, delta_meteo, 0.001,
              4326, cores)
    
    time_end = time.time()
    print('HANTS METEO executed in ', time_end - time_now, 'Seconds')
 
    return(nc_path_lst, nc_path_meteo)
    
    
'''
# Data parameters
rasters_path = r'G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_SEBAL_out\Surface_Temperature'
name_format = 'Surface_Temperature_VIIRS_{yyyy}{mm}{dd}.tif'
start_date = '2016-10-01'
end_date = '2017-09-30'
latlim = [31.9471388888888903, 32.7931388888888904]
lonlim = [-7.2748388888888886, -5.7718388888888885]
cellsize = 0.00375
nc_path = r'G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_SEBAL_out\Thermal_HANTS.nc'
#rasters_path_out = r'G:\SEBAL_Tadla\Output_rasters'

# HANTS parameters
nb = 365
nf = 1
low = 270
high = 350
HiLo = 'Lo'
fet = 10
delta = 1
dod = 1
cores = 3

import time

time_now = time.time()
# Run
run_HANTS(rasters_path, name_format,
          start_date, end_date, latlim, lonlim, cellsize, nc_path,
          nb, nf, HiLo, low, high, fet, dod, delta, 0.001,
          4326, cores)

time_end = time.time()
print('HANTS executed in ', time_end - time_now, 'Seconds')

import watools
watools.Collect.GLDAS.three_hourly("G:\SEBAL_Tadla\Daily_SEBAL\Dataset_in",["avgsurft_inst"], "2016-10-01","2017-09-30", [31.75, 33.25],[-7.50,-5.25], Periods=[4,5])


from SEBAL.hants.wa_gdal import *  # from hants.wa_arcpy import *

# Data parameters
rasters_path = r'G:\SEBAL_Tadla\Daily_SEBAL\Dataset_in\Weather_Data\Model\GLDAS\three_hourly\avgsurft_inst'
name_format = 'T_GLDAS-NOAH_C_3hour_{yyyy}.{mm}.{dd}_5.tif'
start_date = '2016-10-01'
end_date = '2017-09-30'
latlim = [31.75, 33.25]
lonlim = [-7.50,-5.25]
cellsize = 0.25
nc_path = r'G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_SEBAL_out\GLDAS_HANTS5.nc'
#rasters_path_out = r'G:\SEBAL_Tadla\Output_rasters'

# HANTS parameters
nb = 365
nf = 5
low = -3
high = 70
HiLo = 'Lo'
fet = 10
delta = 1
dod = 1
cores = 3

import time
'''


