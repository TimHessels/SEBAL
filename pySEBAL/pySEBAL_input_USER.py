# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Thu Apr  4 12:27:04 2019
"""
from osgeo import gdal
import os
import numpy as np
import datetime

def Get_Time_Info(workbook, number):
    
    ws = workbook['USER_Input']
    
    UTM_Zone = int(ws['G%d' %number].value) 
        
    Time_GMT_String = "%s" %str(ws['H%d' %number].value)      
    hour_GTM = int(Time_GMT_String.split(":")[0])
    minutes_GTM = int(Time_GMT_String.split(":")[1])
    
    Date_string = "%s" %str(ws['I%d' %number].value)       
    try:
       Date_stamp = datetime.datetime.strptime(Date_string, "%d/%m/%Y")
    except:
        Date_stamp = datetime.datetime.strptime(Date_string, "%Y-%m-%d 00:00:00")       
        
    year = Date_stamp.year
    DOY = int(Date_stamp.strftime("%j"))
    
    return(year, DOY, hour_GTM, minutes_GTM, UTM_Zone)

def Get_USER_Para_Veg(workbook, number, Example_fileName, year, month, day):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    output_folder = r"%s" %str(ws['C%d' %number].value)

    ws = workbook['Additional_Input']
    
    # Define the bands that will be used
    res2 = '10m'
    res3 = '10m'    
 
    # Get General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    ######################### Calculate Vegetation Parameters Based on VIS data #####################################

    # Extract the name of the NDVI USER image from the excel file
    ws = workbook['USER_Input']
    src_FileName_NDVI = r"%s" %str(ws['B%d' %number].value)                #NDVI

    # Open the Additional input excel sheet
    ws = workbook['Additional_Input']

    # Check NDVI and Calculate NDVI
    try:
        if (ws['B%d' % number].value) is None:    
            
            # Open User Albedo
            dest_NDVI = SEBAL.reproject_dataset_example(src_FileName_NDVI, Example_fileName)[0]
            NDVI = dest_NDVI.GetRasterBand(1).ReadAsArray()     

            # Create Water mask based on PROBA-V
            water_mask_temp = np.where(NDVI < 0.0, 1.0, 0.0) 
            
        else:
            # Output folder NDVI
            if not os.path.exists(os.path.join(output_folder, 'Output_vegetation')):
                os.makedirs(os.path.join(output_folder, 'Output_vegetation'))
            ndvi_fileName_user = os.path.join(output_folder, 'Output_vegetation', 'User_NDVI_%s_%s%02d%02d.tif' %(res3, year, month, day))
            NDVI = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['B%d' % number].value),ndvi_fileName_user,Example_fileName)
            water_mask_temp = np.where(NDVI < 0.0, 1.0, 0.0)           
            
    except:
        print("Please check the NDVI input path")

    # Check Water Mask and replace if it is filled in the additianal data sheet
    try:
        if (ws['E%d' % number].value) is not None:
            
            if not os.path.exists(os.path.join(output_folder, 'Output_soil_moisture')):
                os.makedirs(os.path.join(output_folder, 'Output_soil_moisture')) 
                
            # Overwrite the Water mask and change the output name
            water_mask_temp_fileName = os.path.join(output_folder, 'Output_soil_moisture', 'User_Water_mask_temporary_%s_%s%02d%02d.tif' %(res2, year, month, day))
            water_mask_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['E%d' % number].value), water_mask_temp_fileName, Example_fileName)
            SEBAL.save_GeoTiff_proy(lsc, water_mask_temp, water_mask_temp_fileName, shape_lsc, nband=1)

    except:
        print("Please check the Water Mask input path")

    # Extract the name of the surface albedo USER image from the excel file
    ws = workbook['USER_Input']
    src_FileName_Surf_albedo = r"%s" %str(ws['C%d' %number].value)                #surface albedo
    
    ws = workbook['Additional_Input']
    # Check Surface albedo
    try:
        if (ws['C%d' % number].value) is not None:

            # Output folder surface albedo
            if not os.path.exists(os.path.join(output_folder, 'Output_vegetation')):
                os.makedirs(os.path.join(output_folder, 'Output_vegetation'))
            surface_albedo_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_albedo_%s_%s%02d%02d.tif' %(res2, year, month, day))
            Surf_albedo=SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['C%d' % number].value),surface_albedo_fileName,Example_fileName)
            SEBAL.save_GeoTiff_proy(lsc, Surf_albedo, surface_albedo_fileName, shape_lsc, nband=1)

        else:

            # Open User Albedo
            dest_Surf_albedo = SEBAL.reproject_dataset_example(src_FileName_Surf_albedo, Example_fileName)[0]
            Surf_albedo = dest_Surf_albedo.GetRasterBand(1).ReadAsArray()     
            
            # Set limit surface albedo
            Surf_albedo = np.minimum(Surf_albedo, 0.6)

    except:
          print("Please check the Albedo input path")

    # calculate vegetation properties
    FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity=SEBAL.Calc_vegt_para(NDVI, water_mask_temp)

    # create quality map
    QC_Map = np.zeros(NDVI.shape)
    QC_Map[np.isnan(NDVI)] = 1

    print('Average NDVI = %s' %np.nanmean(NDVI))
    print('Average Surface Albedo = %s' %np.nanmean(Surf_albedo))
    print('Average LAI = %s' %np.nanmean(LAI))
    print('Average Vegetation Cover = %s' %np.nanmean(vegt_cover))
    print('Average FPAR = %s' %np.nanmean(FPAR))

    return(Surf_albedo, NDVI, LAI, vegt_cover, FPAR, Nitrogen, tir_emis, water_mask_temp, QC_Map)


def Get_USER_Para_Thermal(workbook, number, Example_fileName, year, month, day, water_mask_temp, surf_temp_offset, Thermal_Sharpening_not_needed):
    
    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    output_folder = r"%s" %str(ws['C%d' %number].value)

    # Open General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    # General sensor resolution and names'
    res2 = '250m'

    # Extract the name of the thermal and quality MODIS image from the excel file
    ws = workbook['USER_Input']
    src_FileName_LST = r"%s" %str(ws['D%d' %number].value)                #land surface temperature

    ws = workbook['Additional_Input']

    try:
        # If all additional fields are filled in than do not open the datasets
        if ws['D%d' % number].value is None:

            print('...................... Open USER Thermal ........................')

            # Open User LST
            dest_Surface_temp = SEBAL.reproject_dataset_example(src_FileName_LST, Example_fileName)[0]
            Surface_temp = dest_Surface_temp.GetRasterBand(1).ReadAsArray()
            Thermal_Sharpening_not_needed = 1
            
        else:

            # Output folder surface temperature
            surf_temp_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_temp_%s_%s%02d%02d.tif' %(res2, year, month, day))
            Surface_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['D%d' % number].value), surf_temp_fileName, Example_fileName)
            Thermal_Sharpening_not_needed = 1

    except:
        assert "Please check the surface temperature input path"


    # Cloud mask:
    temp_water = np.zeros((shape_lsc[1], shape_lsc[0]))
    temp_water = np.copy(Surface_temp)
    temp_water[water_mask_temp == 0.0] = np.nan
    temp_water_sd = np.nanstd(temp_water)     # Standard deviation
    temp_water_mean = np.nanmean(temp_water)  # Mean
    print('Mean water temperature = ', '%0.3f (Kelvin)' % temp_water_mean)
    print('SD water temperature = ', '%0.3f (Kelvin)' % temp_water_sd)
    cloud_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
    cloud_mask_temp[Surface_temp < np.minimum((temp_water_mean - 1.0 * temp_water_sd -
               surf_temp_offset),290)] = 1.0
    cloud_mask_temp[np.isnan(Surface_temp)] = 1

    # remove wrong values VIIRS defined by user
    Surface_temp[cloud_mask_temp == 1] = np.nan
    print('Mean Surface Temperature = %s Kelvin' %np.nanmean(Surface_temp))
    
    return(Surface_temp, cloud_mask_temp, Thermal_Sharpening_not_needed)