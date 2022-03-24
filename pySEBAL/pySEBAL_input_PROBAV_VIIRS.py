# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:26:35 2018

@author: tih
"""

import time
import os
import re
from osgeo import gdal, osr
import numpy as np
import datetime


def Get_Time_Info(workbook, number):

    # Open the VIIRS_PROBAV_Input sheet
    ws = workbook['VIIRS_PROBAV_Input']

    # Extract the name of the thermal and quality VIIRS image from the excel file
    Name_VIIRS_Image_TB = '%s' %str(ws['B%d' %number].value)
    Name_VIIRS_Image_QC = '%s' %str(ws['C%d' %number].value)

    # Extract the name to the PROBA-V image from the excel file
    Name_PROBAV_Image = '%s' %str(ws['D%d' %number].value)    # Must be a tiff file

    # UTM Zone of the end results
    UTM_Zone = float(ws['G%d' %number].value)

	# Get time from the VIIRS dataset name (IMPORTANT TO KEEP THE TEMPLATE OF THE VIIRS NAME CORRECT example: VIIRS_SVIO5_npp_20160601_1103128_e1108532_b23808_c20160601170854581426_noaa_ops.tif npp_viirs_i05_20150701_124752_wgs84_fit.tif)
    #Total_Day_VIIRS = Name_VIIRS_Image_TB.split('_')[3]
    #Total_Time_VIIRS = Name_VIIRS_Image_TB.split('_')[4]
    # "NPP_VIAES_L1.A2018001.0854.001.2018009150807.pssgrpgs_000501315598.BrightnessTemperature_I5-BrightnessTemperature_I5.tif"
    Total_Day_VIIRS = datetime.datetime.strptime(Name_VIIRS_Image_TB.split('.')[1][1:], "%Y%j").strftime("%Y%m%d")
    Total_Time_VIIRS = Name_VIIRS_Image_TB.split('.')[2]
    
    # Get the information out of the VIIRS name
    year = int(Total_Day_VIIRS[:4])
    month = int(Total_Day_VIIRS[4:6])
    day = int(Total_Day_VIIRS[6:8])
    Startdate = '%d-%02d-%02d' % (year,month,day)
    DOY = datetime.datetime.strptime(Startdate,'%Y-%m-%d').timetuple().tm_yday
    hour = int(Total_Time_VIIRS[0:2])
    minutes = int(Total_Time_VIIRS[2:4])

    # Print data used from sheet General_Input
    print('VIIRS PROBA-V Input:')
    print('Path to Thermal VIIRS image = %s' %str(Name_VIIRS_Image_TB))
    print('Path to Quality VIIRS image = %s' %str(Name_VIIRS_Image_QC))
    print('Name of PROBA-V image = %s' %str(Name_PROBAV_Image))

    return(year, DOY, hour, minutes, UTM_Zone)

def Get_PROBAV_Para_Veg(workbook, number, Example_fileName, year, month, day, path_radiance, Apparent_atmosf_transm, cos_zn, dr, DEM_resh):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)

    ws = workbook['Additional_Input']
    
    # Define the bands that will be used
    sensor1 = 'PROBAV'
    sensor2 = 'VIIRS'
    res1 = '375m'
    res2 = '100m'
    res3 = '30m'    

    # If all additional fields are filled in than do not open the datasets
    if ws['B%d' % number].value is None or ws['C%d' % number].value is None:

        print('--------------------- Open PROBA-V VIS ------------------------')

        # Open the Landsat_Input sheet
        ws = workbook['VIIRS_PROBAV_Input']

        # read PROBAV name input
        Name_PROBAV_Image = '%s' %str(ws['D%d' %number].value)    # Must be a tiff or hdf file

##############################################################################################################################################
############################################# maak hier een functie van zodat deze ook in preSEBAL_1.py gebruikt kan worden ##################
##############################################################################################################################################

        spectral_reflectance_PROBAV, cloud_mask_temp = Open_PROBAV_Reflectance(Name_PROBAV_Image, input_folder, output_folder, Example_fileName)
  
##############################################################################################################################################
############################    RETURN spectral_reflectance_PROBAV, cloud_mask_temp      #####################################################
##############################################################################################################################################
    
    # Get General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    ######################### Calculate Vegetation Parameters Based on VIS data #####################################

    # Open the Additional input excel sheet
    ws = workbook['Additional_Input']

    # Check NDVI and Calculate NDVI
    try:
        if (ws['B%d' % number].value) is None:          
            n218_memory = spectral_reflectance_PROBAV[:, :, 2] + spectral_reflectance_PROBAV[:, :, 3]
            NDVI = np.zeros((shape_lsc[1], shape_lsc[0]))
            NDVI[n218_memory != 0] =  ( spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] - spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] )/ ( spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] + spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] )
 
            # Create Water mask based on PROBA-V
            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
            water_mask_temp[np.logical_and(spectral_reflectance_PROBAV[:, :, 2] >= spectral_reflectance_PROBAV[:, :, 3],DEM_resh>0)]=1
            
        else:
            # Output folder NDVI
            if not os.path.exists(os.path.join(output_folder, 'Output_vegetation')):
                os.makedirs(os.path.join(output_folder, 'Output_vegetation'))
            ndvi_fileName_user = os.path.join(output_folder, 'Output_vegetation', 'User_NDVI_%s_%s%02d%02d.tif' %(res3, year, month, day))
            NDVI = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['B%d' % number].value),ndvi_fileName_user,Example_fileName)
            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))            
            water_mask_temp[NDVI < 0.0] = 1.0               
            SEBAL.save_GeoTiff_proy(lsc, NDVI, ndvi_fileName_user, shape_lsc, nband=1)
         
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
            # Calculate surface albedo based on PROBA-V
            # Surf_albedo = 0.429 * spectral_reflectance_PROBAV[:, :, 1] + 0.333 * spectral_reflectance_PROBAV[:, :, 2] + 0.133 * spectral_reflectance_PROBAV[:, :, 3] + 0.105 * spectral_reflectance_PROBAV[:, :, 4]  # WAPOR          
            Surf_albedo = 0.219 * spectral_reflectance_PROBAV[:, :, 1] + 0.361 * spectral_reflectance_PROBAV[:, :, 2] + 0.379 * spectral_reflectance_PROBAV[:, :, 3] + 0.041 * spectral_reflectance_PROBAV[:, :, 4]
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

    return(Surf_albedo, NDVI, LAI, vegt_cover, FPAR, Nitrogen, tir_emis, b10_emissivity, water_mask_temp, QC_Map)

def Get_VIIRS_Para_Thermal(workbook, number, Example_fileName, year, month, day, water_mask_temp, b10_emissivity, Temp_inst,  Rp, tau_sky, surf_temp_offset, Thermal_Sharpening_not_needed):

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

    # Open the VIIRS_PROBAV_Input sheet
    ws = workbook['VIIRS_PROBAV_Input']

    # Extract VIIRS name from input excel
    Name_VIIRS_Image = str(ws['B%d' %number].value)
    Name_VIIRS_Image_QC = str(ws['C%d' %number].value)
    k1=606.399172
    k2=1258.78
    sensor1 = 'PROBAV'
    sensor2 = 'VIIRS'
    res1 = '375m'
    res2 = '100m'
    res3 = '30m'

    # Open the Landsat_Input sheet
    ws = workbook['Additional_Input']

    try:
        # If all additional fields are filled in than do not open the datasets
        if ws['D%d' % number].value is None:

            print('...................... Open VIIRS Thermal ........................')

             # Define the VIIRS thermal data name
            VIIRS_data_name=os.path.join(input_folder, '%s' % (Name_VIIRS_Image))

            # Reproject VIIRS thermal data
            VIIRS, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset_example(
       	                        VIIRS_data_name, Example_fileName)

            # Open VIIRS thermal data
            data_VIIRS = VIIRS.GetRasterBand(1).ReadAsArray()

            # Set the conditions for the brightness temperature
            brightness_temp=np.where(data_VIIRS>=250, data_VIIRS, np.nan)

            # Constants
            L_lambda_b10 = ((2*6.63e-34*(3.0e8)**2)/((11.45e-6)**5*(np.exp((6.63e-34*3e8)/(1.38e-23*(11.45e-6)*brightness_temp))-1)))*1e-6

            # Get Surface Temperature
            Surface_temp = SEBAL.Get_Thermal(L_lambda_b10, Rp, Temp_inst, tau_sky, b10_emissivity, k1, k2)
            Surface_temp = Surface_temp.clip(230.0, 360.0)

        else:

            # Output folder surface temperature
            surf_temp_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_temp_%s_%s%02d%02d.tif' %(res2, year, month, day))
            Surface_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['D%d' % number].value),surf_temp_fileName, Example_fileName)
            Thermal_Sharpening_not_needed = 1

    except:
        assert "Please check the surface temperature input path"

    try:
        # Check Quality
        if Name_VIIRS_Image_QC != 'None':

            # Define the VIIRS Quality data name
            VIIRS_data_name=os.path.join(input_folder, '%s' % (Name_VIIRS_Image_QC))

            # Reproject VIIRS Quality data
            dest_VIIRS_QC, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset_example(
                                           VIIRS_data_name, Example_fileName)

            # Open VIIRS Quality data
            cloud_mask_temp = dest_VIIRS_QC.GetRasterBand(1).ReadAsArray()

        else:

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
    except:
        assert "Please check the QC temperature input path"

    # remove wrong values VIIRS defined by user
    Surface_temp[cloud_mask_temp == 1] = np.nan
    print('Mean Surface Temperature = %s Kelvin' %np.nanmean(Surface_temp))

    return(Surface_temp, cloud_mask_temp, Thermal_Sharpening_not_needed)
    
def Open_PROBAV_Reflectance(Name_PROBAV_Image, input_folder, output_folder, Example_fileName):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL
    
    # Define bands probav
    bands=['SM', 'B1', 'B2', 'B3', 'B4']  #'SM', 'BLUE', 'RED', 'NIR', 'SWIR'
    
    # Get extension of PROBAV dataset
    Name_PROBAV, Name_PROBAV_exe = os.path.splitext(Name_PROBAV_Image)

    # Set the index number at 0
    index=0

    # constants
    n188_float = 248       

    if (Name_PROBAV_exe == '.hdf5' or Name_PROBAV_exe == '.HDF5'):

        # write the data one by one to the spectral_reflectance_PROBAV
        for bandnmr in bands:

            # Translate the PROBA-V names to the Landsat band names
            Band_number = {'SM':7,'B1':8,'B2':10,'B3':9,'B4':11}

            # Open the hdf file
            Band_PROBAVhdf_fileName = os.path.join(input_folder, Name_PROBAV_Image)
            g = gdal.Open(Band_PROBAVhdf_fileName, gdal.GA_ReadOnly)

            #  Define temporary file out and band name in
            name_out = os.path.join(input_folder, '%s_test.tif' % (Name_PROBAV_Image))
            name_in = g.GetSubDatasets()[Band_number[bandnmr]][0]

            # Get environmental variable
            SEBAL_env_paths = os.environ["SEBAL"].split(';')
            GDAL_env_path = SEBAL_env_paths[0]
            GDAL_TRANSLATE = os.path.join(GDAL_env_path, 'gdal_translate.exe')

            # run gdal translate command
            FullCmd = '%s -of GTiff %s %s' %(GDAL_TRANSLATE, name_in, name_out)
            SEBAL.Run_command_window(FullCmd)

            # Get the data array
            dest_PV = gdal.Open(name_out)
            Data = dest_PV.GetRasterBand(1).ReadAsArray()
            dest_PV = None
 
            # Define the x and y spacing
            Meta_data = g.GetMetadata()
            Lat_Top = float(Meta_data['LEVEL3_GEOMETRY_TOP_RIGHT_LATITUDE'])
            Lon_Left = float(Meta_data['LEVEL3_GEOMETRY_BOTTOM_LEFT_LONGITUDE'])
            Pixel_size = float((Meta_data['LEVEL3_GEOMETRY_VNIR_VAA_MAPPING']).split(' ')[-3])
            x_size_pv = int(Data.shape[1])
            y_size_pv = int(Data.shape[0])
            ulx = Lon_Left - 0.5*Pixel_size
            uly = Lat_Top + 0.5*Pixel_size
            lrx = ulx + x_size_pv * Pixel_size
            lry = uly - y_size_pv * Pixel_size      
            
            # Define the georeference of the PROBA-V data
            geo_PROBAV=[Lon_Left-0.5*Pixel_size, Pixel_size, 0, Lat_Top+0.5*Pixel_size, 0, -Pixel_size] #0.000992063492063

            ################################# Create a MEMORY file ##############################
            # create memory output with the PROBA-V band
            fmt = 'MEM'
            driver = gdal.GetDriverByName(fmt)
            dst_dataset = driver.Create('', int(Data.shape[1]), int(Data.shape[0]), 1,gdal.GDT_Float32)
            dst_dataset.SetGeoTransform(geo_PROBAV)

            # set the reference info
            srs = osr.SpatialReference()
            srs.SetWellKnownGeogCS("WGS84")
            dst_dataset.SetProjection(srs.ExportToWkt())

            # write the array in the geotiff band
            dst_dataset.GetRasterBand(1).WriteArray(Data)

            ################################# reproject PROBAV MEMORY file ##############################

            # Reproject the PROBA-V band  to match DEM's resolution
            PROBAV, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset_example(
                          dst_dataset, Example_fileName)

            dst_dataset = None

            #################################### Get example information ################################

            if not "shape_lsc" in locals():
                nrow = PROBAV.RasterYSize
                ncol = PROBAV.RasterXSize
                shape_lsc = [ncol, nrow]

            # Open the reprojected PROBA-V band data
            data_PROBAV_DN = PROBAV.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)

            # Define the filename to store the cropped Landsat image
            dst_FileName = os.path.join(output_folder, 'Output_radiation_balance','proy_PROBAV_%s.tif' % bandnmr)

            # close the PROBA-V
            g=None

            if not "spectral_reflectance_PROBAV" in locals():
                spectral_reflectance_PROBAV=np.zeros([shape_lsc[1], shape_lsc[0], 5])

            # If the band data is the SM band than write the data into the spectral_reflectance_PROBAV and create cloud mask
            if bandnmr is 'SM':
                cloud_mask_temp = np.zeros(data_PROBAV_DN.shape)
                cloud_mask_temp[data_PROBAV_DN[:,:]!=n188_float]=1
                spectral_reflectance_PROBAV[:, :, index] = cloud_mask_temp

            # If the band data is not SM change the DN values into PROBA-V values and write into the spectral_reflectance_PROBAV
            else:
                data_PROBAV=data_PROBAV_DN/2000
                spectral_reflectance_PROBAV[:, :, index]=data_PROBAV[:, :]

            # Change the spectral reflectance to meet certain limits
            spectral_reflectance_PROBAV[:, :, index]=np.where(spectral_reflectance_PROBAV[:, :, index]<=0,np.nan,spectral_reflectance_PROBAV[:, :, index])
            spectral_reflectance_PROBAV[:, :, index]=np.where(spectral_reflectance_PROBAV[:, :, index]>=150,np.nan,spectral_reflectance_PROBAV[:, :, index])

            # Save the PROBA-V as a tif file
            SEBAL.save_GeoTiff_proy(PROBAV, spectral_reflectance_PROBAV[:, :, index], dst_FileName, shape_lsc, nband=1)

            # Go to the next index
            index=index+1
            
            # Remove temporary file
            os.remove(name_out)
            
    else:
        for bandnmr in bands:

            # Open information
            Band_PROBAVhdf_fileName = os.path.join(input_folder, "%s_%s.tif" %(Name_PROBAV_Image, bandnmr))
            g = gdal.Open(Band_PROBAVhdf_fileName, gdal.GA_ReadOnly) 
            
            # Open original raster band
            x_size_pv = g.RasterXSize
            y_size_pv = g.RasterYSize
            geo_out = g.GetGeoTransform()
            ulx = geo_out[0]
            uly = geo_out[3]
            xDist = geo_out[1]
            yDist = geo_out[5]   
            lrx = ulx + x_size_pv * xDist
            lry = uly - y_size_pv * yDist                  
            
            # Reproject the PROBA-V band  to match DEM's resolution
            PROBAV, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset_example(
                          g, Example_fileName)
        
            data_PROBAV = PROBAV.GetRasterBand(1).ReadAsArray()
            
            # Define example array properties
            if not "shape_lsc" in locals():
                nrow = PROBAV.RasterYSize
                ncol = PROBAV.RasterXSize
                shape_lsc = [ncol, nrow]
            
            # Define spectal reflactance empty array
            if not "spectral_reflectance_PROBAV" in locals():
                spectral_reflectance_PROBAV=np.zeros([shape_lsc[1], shape_lsc[0], 5])
   

            # If the band data is the SM band than write the data into the spectral_reflectance_PROBAV and create cloud mask
            if bandnmr is 'SM':
                cloud_mask_temp = np.zeros(data_PROBAV.shape)
                cloud_mask_temp[data_PROBAV[:,:]!=n188_float]=1
                spectral_reflectance_PROBAV[:, :, index] = cloud_mask_temp
         
            # If the band data is not SM change the DN values into PROBA-V values and write into the spectral_reflectance_PROBAV
            else:
                spectral_reflectance_PROBAV[:, :, index]=data_PROBAV[:, :]/2000

            # Define the filename to store the cropped Landsat image
            dst_FileName = os.path.join(output_folder, 'Output_radiation_balance','proy_PROBAV_%s.tif' % bandnmr)

            # Define shape
            if not "shape_lsc" in locals():
                nrow = PROBAV.RasterYSize
                ncol = PROBAV.RasterXSize
                shape_lsc = [ncol, nrow]   
          
            SEBAL.save_GeoTiff_proy(PROBAV, spectral_reflectance_PROBAV[:, :, index], dst_FileName, shape_lsc, nband=1)
        
            # Go to the next index
            index=index+1    

    # Original size PROBAV dataset
    print('Original PROBA-V Image - ')
    print('  Size :', x_size_pv, y_size_pv)
    print('  Upper Left corner x, y: ', ulx, ', ', uly)
    print('  Lower right corner x, y: ', lrx, ', ', lry)

    print('Reprojected PROBA-V Image - ')
    print('  Size :', shape_lsc[1], shape_lsc[0])
    print('  Upper Left corner x, y: ', ulx_dem, ', ', uly_dem)
    print('  Lower right corner x, y: ', lrx_dem, ', ', lry_dem)
    
    return(spectral_reflectance_PROBAV, cloud_mask_temp)    
    
    
    