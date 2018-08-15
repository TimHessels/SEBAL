# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:27:54 2018

@author: tih
"""
import os
import gdal
import numpy as np

def Get_Time_Info(workbook, number):

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)

    # Open the VIIRS_PROBAV_Input sheet
    ws = workbook['MODIS_Input']

    # Extract the name of the thermal and quality VIIRS image from the excel file
    Name_MODIS_Image_Ref = str(ws['B%d' %number].value)                #reflectance
    Name_MODIS_Image_NDVI = str(ws['D%d' %number].value)               #ndvi
    Name_MODIS_Image_LST = str(ws['C%d' %number].value)                #land surface temperature

    # Create complete path to data
    src_FileName_LST = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_LST)
    src_FileName_NDVI = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_NDVI)
    src_FileName_Ref = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_Ref)

    # UTM Zone of the end results
    UTM_Zone = float(ws['G%d' %number].value)

	 #Get time from the MODIS dataset name (IMPORTANT TO KEEP THE TEMPLATE OF THE MODIS NAME CORRECT example: MOD13Q1.A2008129.h18v05.006.2015175090913.hdf)
    Total_Day_MODIS = Name_MODIS_Image_LST.split('.')[-4][1:]

    # Get the information out of the VIIRS name
    year = int(Total_Day_MODIS[0:4])
    DOY =  int(Total_Day_MODIS[4:7])

    # Print data used from sheet General_Input
    print('MODIS Input:')
    print('Path to MODIS LST image = %s' %str(src_FileName_LST))
    print('Path to MODIS NDVI image = %s' %str(src_FileName_NDVI))
    print('Path to MODIS Reflectance image = %s' %str(src_FileName_Ref))

    return(year, DOY, UTM_Zone)

def Get_MODIS_Para_Veg(workbook, number, Example_fileName, year, DOY, path_radiance, Apparent_atmosf_transm, cos_zn, dr, DEM_resh, epsg_to):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)

    # Open General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    # General sensor resolution and names
    res2 = '250m'

    # Extract the name of the thermal and quality MODIS image from the excel file
    ws = workbook['MODIS_Input']
    Name_MODIS_Image_Ref = str(ws['B%d' %number].value)                #reflectance
    Name_MODIS_Image_NDVI = str(ws['D%d' %number].value)               #ndvi

    # Create complete path to data
    src_FileName_NDVI = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_NDVI)
    src_FileName_Ref = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_Ref)

    ws = workbook['Additional_Input']

    # Calculate NDVI
    try:
        if (ws['B%d' % number].value) is not None:

            # Output folder NDVI	defined by the user
            ndvi_fileName = os.path.join(output_folder, 'Output_vegetation', 'User_NDVI_%s_%s_%s.tif' %(res2, year, DOY))

            # Reproject and reshape users NDVI
            NDVI = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['B%d' % number].value),ndvi_fileName, Example_fileName)
            NDVI_MAX = np.nanmax(NDVI)
            NDVI_SD =  np.nanstd(NDVI)
            print('NDVI User max ' , NDVI_MAX)
            print('NDVI User sd' , NDVI_SD)

            # Create Water mask based on PROBA-V
            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
            water_mask_temp[NDVI<0.0]=1

            # Create Quality map
            QC_Map = np.zeros((shape_lsc[1], shape_lsc[0]))

        else:
            # Calculate the NDVI based on MODIS
            NDVI = Open_reprojected_hdf(src_FileName_NDVI, 0, epsg_to, 0.0001, Example_fileName)

            # NDVI stats
            NDVI_MODIS_MAX = np.nanmax(NDVI)
            NDVI_MODIS_SD =  np.nanstd(NDVI)
            print('NDVI MODIS max ' , NDVI_MODIS_MAX)
            print('NDVI MODIS sd' , NDVI_MODIS_SD)

            # Create Water mask based on MODIS
            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
            water_mask_temp[NDVI < 0]=1

            # Create Quality map
            QC_Map = np.zeros((shape_lsc[1], shape_lsc[0]))
            QC_Map[np.logical_and(NDVI<-0.2, NDVI>1.0)] = 1
    except:
        assert "Please check the MODIS path, was not able to create NDVI"

    # Check Water Mask and replace the temporary water mask if needed
    try:
        if (ws['E%d' % number].value) is not None:
            # Overwrite the Water mask
             water_mask_temp_fileName = os.path.join(output_folder, 'Output_soil_moisture', 'User_Water_mask_temporary_%s_%s_%s.tif' %(res2, year, DOY))
             water_mask_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['E%d' % number].value), water_mask_temp_fileName, Example_fileName)
             SEBAL.save_GeoTiff_proy(lsc, water_mask_temp, water_mask_temp_fileName, shape_lsc, nband=1)
    except:
        assert "Please check the Water Mask input path"


    # Check surface albedo
    try:
        if (ws['C%d' % number].value) is not None:
            # Output folder surface albedo
            surface_albedo_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_albedo_%s_%s_%s.tif' %(res2, year, DOY))

            # Reproject and reshape users surface albedo
            Surf_albedo=SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['C%d' % number].value), surface_albedo_fileName, Example_fileName)

        # if the users surface albedo data cannot be reprojected than use the original MODIS data as imported into SEBAL
        else:

            # Calculate the MOD9 based on MODIS by opening and reproject bands
            B1_modis = Open_reprojected_hdf(src_FileName_Ref, 11, epsg_to, 0.0001, Example_fileName)
            B2_modis = Open_reprojected_hdf(src_FileName_Ref, 12, epsg_to, 0.0001, Example_fileName)
            B3_modis = Open_reprojected_hdf(src_FileName_Ref, 13, epsg_to, 0.0001, Example_fileName)
            B4_modis = Open_reprojected_hdf(src_FileName_Ref, 14, epsg_to, 0.0001, Example_fileName)
            B5_modis = Open_reprojected_hdf(src_FileName_Ref, 15, epsg_to, 0.0001, Example_fileName)
            B6_modis = Open_reprojected_hdf(src_FileName_Ref, 16, epsg_to, 0.0001, Example_fileName)
            B7_modis = Open_reprojected_hdf(src_FileName_Ref, 17, epsg_to, 0.0001, Example_fileName)

            # Calc surface albedo within shortwave domain using a weighting function (Tasumi et al 2008)
            Surf_albedo = 0.215 * B1_modis + 0.215 * B2_modis + 0.242 * B3_modis + 0.129 * B4_modis + 0.101 * B5_modis + 0.062 * B6_modis + 0.036 * B7_modis

            # Define bad pixels
            QC_Map[np.logical_and(Surf_albedo < 0.0, Surf_albedo > 1.0)] = 1

    except:
         assert "Please check the PROBA-V path, was not able to create Albedo"

    # Calculate the Fpar, TIR, Nitrogen, Vegetation Cover, LAI and b10_emissivity based on PROBA-V
    FPAR, tir_emis, Nitrogen, vegt_cover, LAI, b10_emissivity = SEBAL.Calc_vegt_para(NDVI, water_mask_temp, shape_lsc)

    print('Average NDVI = %s' %np.nanmean(NDVI))
    print('Average Surface Albedo = %s' %np.nanmean(Surf_albedo))
    print('Average LAI = %s' %np.nanmean(LAI))
    print('Average Vegetation Cover = %s' %np.nanmean(vegt_cover))
    print('Average FPAR = %s' %np.nanmean(FPAR))

    return(Surf_albedo, NDVI, LAI, vegt_cover, FPAR, Nitrogen, tir_emis, b10_emissivity, water_mask_temp, QC_Map)


def Get_MODIS_Para_Thermal(workbook, number, Example_fileName, year, DOY, water_mask_temp, b10_emissivity, Temp_inst,  Rp, tau_sky, surf_temp_offset, Thermal_Sharpening_not_needed, epsg_to):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    output_folder = r"%s" %str(ws['C%d' %number].value)
    input_folder = r"%s" %str(ws['B%d' %number].value)

    # Open General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    # General sensor resolution and names'
    res2 = '250m'

    # Extract the name of the thermal and quality MODIS image from the excel file
    ws = workbook['MODIS_Input']
    Name_MODIS_Image_LST = str(ws['C%d' %number].value)                #land surface temperature
    src_FileName_LST = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_LST)

    ws = workbook['Additional_Input']

    try:
        # If all additional fields are filled in than do not open the datasets
        if ws['D%d' % number].value is None:

            print('...................... Open MODIS Thermal ........................')

            # Calculate the MOD9 based on MODIS
            Surface_temp = Open_reprojected_hdf(src_FileName_LST, 0, epsg_to, 0.02, Example_fileName)

        else:

            # Output folder surface temperature
            surf_temp_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_temp_%s_%s_%s.tif' %(res2, year, DOY))
            Surface_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['D%d' % number].value), surf_temp_fileName, Example_fileName)
            Thermal_Sharpening_not_needed = 0

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

    # remove wrong values VIIRS defined by user
    Surface_temp[cloud_mask_temp == 1] = np.nan
    print('Mean Surface Temperature = %s Kelvin' %np.nanmean(Surface_temp))

    return(Surface_temp, cloud_mask_temp, Thermal_Sharpening_not_needed)


#------------------------------------------------------------------------------
def reproject_MODIS(input_name, output_name, epsg_to):

    '''
    Reproject the merged data file

    Keywords arguments:
    output_folder -- 'C:/file/to/path/'
    '''
    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Get environmental variable
    SEBAL_env_paths = os.environ["SEBAL"].split(';')
    GDAL_env_path = SEBAL_env_paths[0]
    GDALWARP_PATH = os.path.join(GDAL_env_path, 'gdalwarp.exe')

    split_input = input_name.split('hdf":')
    inputname = '%shdf":"%s"' %(split_input[0],split_input[1])

    # find path to the executable
    fullCmd = ' '.join(["%s" %(GDALWARP_PATH), '-overwrite -s_srs "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"', '-t_srs EPSG:%s -of GTiff' %(epsg_to), inputname, output_name])
    SEBAL.Run_command_window(fullCmd)

    return()

#------------------------------------------------------------------------------
def Open_reprojected_hdf(input_name, Band, epsg_to, scale_factor, Example_fileName):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    g=gdal.Open(input_name, gdal.GA_ReadOnly)

    folder_out = os.path.dirname(input_name)
    output_name_temp = os.path.join(folder_out, "temporary.tif")

    # Open and reproject
    name_in = g.GetSubDatasets()[Band][0]

    reproject_MODIS(name_in, output_name_temp, epsg_to)

    dest, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset_example(output_name_temp, Example_fileName)
    Array = dest.GetRasterBand(1).ReadAsArray() * scale_factor

    os.remove(output_name_temp)

    return(Array)

#------------------------------------------------------------------------------

def Modis_Time(workbook, epsg_to, number, Example_fileName):

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)

    # Open the VIIRS_PROBAV_Input sheet
    ws = workbook['MODIS_Input']

    # Extract the name of the thermal and quality VIIRS image from the excel file
    Name_MODIS_Image_LST = str(ws['C%d' %number].value)                #land surface temperature

    # Create complete path to data
    src_FileName_LST = os.path.join(input_folder, '%s.hdf' %Name_MODIS_Image_LST)

    Time = Open_reprojected_hdf(src_FileName_LST, 2, epsg_to, 0.1, Example_fileName)

    hour = np.floor(Time)
    hour[hour == 0] = np.nan
    minutes = (Time - hour) * 60

    return(hour, minutes)