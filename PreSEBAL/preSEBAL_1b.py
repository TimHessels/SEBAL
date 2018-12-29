# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:48:36 2018

@author: tih
"""
def main(VegetationExcel, cores = 1):

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
    rasters_path_albedo = os.path.join(input_folder, "PreSEBAL_SEBAL_out","Albedo")
    name_format_albedo = 'Albedo_PROBAV_{yyyy}{mm}{dd}.tif'
    rasters_path_ndvi = os.path.join(input_folder, "PreSEBAL_SEBAL_out","NDVI")
    name_format_ndvi = 'NDVI_PROBAV_{yyyy}{mm}{dd}.tif'    
    
    # Find extend
    os.chdir(rasters_path_albedo)
    file_example_PV = glob.glob("*.tif")[0]
    filename_in_ex_PV = os.path.join(rasters_path_albedo, file_example_PV)
    
    dest = gdal.Open(filename_in_ex_PV)
    xsize = dest.RasterXSize
    ysize = dest.RasterYSize
    geo = dest.GetGeoTransform()
    
    latlim = [geo[3] + geo[5] * ysize, geo[3]]
    lonlim = [geo[0], geo[0] + geo[1] * xsize]
    
    cellsize = geo[1]
    
    # HANTS general paramters
    nb = len(Dates)
    
    # HANTS parameters Albedo
    nf_alb = int(ws_inputs['B2'].value)   # number of frequencies to be considered above the zero frequency
    low_alb = float(ws_inputs['B3'].value) # valid range minimum
    high_alb = float(ws_inputs['B4'].value) # valid range maximum
    HiLo_alb = str(ws_inputs['B5'].value) # 2-character string indicating rejection of high or low outliers
    fet_alb = float(ws_inputs['B6'].value)  # fit error tolerance (point eviating more than fet from curve fit are rejected)
    delta_alb = float(ws_inputs['B7'].value) # small positive number e.g. 0.1 to supress high amplitudes
    dod_alb = float(ws_inputs['B8'].value) # degree of overdeterminedness (iteration stops if number of points reaches the minimum required for curve fitting, plus dod). This is a safety measure

    # HANTS parameters NDVI 
    nf_ndvi = int(ws_inputs['C2'].value)   # number of frequencies to be considered above the zero frequency
    low_ndvi = float(ws_inputs['C3'].value) # valid range minimum
    high_ndvi = float(ws_inputs['C4'].value) # valid range maximum
    HiLo_ndvi = str(ws_inputs['C5'].value) # 2-character string indicating rejection of high or low outliers
    fet_ndvi = float(ws_inputs['C6'].value)  # fit error tolerance (point eviating more than fet from curve fit are rejected)
    delta_ndvi = float(ws_inputs['C7'].value) # small positive number e.g. 0.1 to supress high amplitudes
    dod_ndvi = float(ws_inputs['C8'].value) # degree of overdeterminedness (iteration stops if number of points reaches the minimum required for curve fitting, plus dod). This is a safety measure

    # Define outputs HANTS    
    nc_path_alb = os.path.join(input_folder, "PreSEBAL_SEBAL_out", "Albedo_HANTS_nf%s_l%s_h%s_%s_fet%s_delta%s_dod%s.nc" %(nf_alb, low_alb, high_alb, HiLo_alb, fet_alb, delta_alb, dod_alb))  
    nc_path_ndvi = os.path.join(input_folder, "PreSEBAL_SEBAL_out", "NDVI_HANTS_nf%s_l%s_h%s_%s_fet%s_delta%s_dod%s.nc" %(nf_ndvi, low_ndvi, high_ndvi, HiLo_ndvi, fet_ndvi, delta_ndvi, dod_ndvi))    

    import time
    
    time_now = time.time()
    # Run
    HANTS.run_HANTS(rasters_path_albedo, name_format_albedo,
              start_date, end_date, latlim, lonlim, cellsize, nc_path_alb,
              nb, nf_alb, HiLo_alb, low_alb, high_alb, fet_alb, dod_alb, delta_alb, 0.001,
              4326, cores)
    
    time_end = time.time()
    print('HANTS ALBEDO executed in ', time_end - time_now, 'Seconds')

    time_now = time.time()    
       
    # Run
    HANTS.run_HANTS(rasters_path_ndvi, name_format_ndvi,
              start_date, end_date, latlim, lonlim, cellsize, nc_path_ndvi,
              nb, nf_ndvi, HiLo_ndvi, low_ndvi, high_ndvi, fet_alb, dod_ndvi, delta_ndvi, 0.001,
              4326, cores)
    
    time_end = time.time()
    print('HANTS NDVI executed in ', time_end - time_now, 'Seconds')

    return(nc_path_alb, nc_path_ndvi)




