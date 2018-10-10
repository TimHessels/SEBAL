# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:24:48 2018
test
@author: tih
"""
import time
import os
import re
import gdal
import numpy as np

def Get_Time_Info(workbook, number):

   # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)

    # Open the Landsat_Input sheet
    ws = workbook['Landsat_Input']

    # Extract Landsat name, number and amount of thermal bands from excel file
    Name_Landsat_Image = str(ws['B%d' %number].value)
    Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)

    # the path to the MTL file of landsat
    Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' % Name_Landsat_Image)

    # read out the general info out of the MTL file
    year, DOY, hour, minutes, UTM_Zone, Sun_elevation = info_general_metadata(Landsat_meta_fileName) # call definition info_general_metadata

    return(year, DOY, hour, minutes, UTM_Zone, Sun_elevation, Landsat_nr)


def Get_LS_Para_Veg(workbook, number, Example_fileName, year, DOY, path_radiance, Apparent_atmosf_transm, cos_zn, dr):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)

    ws = workbook['Additional_Input']

    # If all additional fields are filled in than do not open the datasets
    if ws['B%d' % number].value is None or ws['C%d' % number].value is None:

        print('-------------------- Open Landsat VIS -----------------------')

        # Open the Landsat_Input sheet
        ws = workbook['Landsat_Input']

        # Extract Landsat name, number and amount of thermal bands from excel file
        Name_Landsat_Image = str(ws['B%d' %number].value)
        Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)
                                                             # temperature: 1 = Band 6 for LS_5 & 7, Band 10 for LS_8 (optional)
        # Define bands used for each Landsat number
        if Landsat_nr == 5 or Landsat_nr == 7:
            Bands = np.array([1, 2, 3, 4, 5, 7, 6])
        elif Landsat_nr == 8:
           Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        # Open MTL landsat and get the correction parameters
        Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' %Name_Landsat_Image)
        Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)
        print('Lmin= ', Lmin)
        print('Lmax= ', Lmax)
        print('k1= ', k1_c)
        print('k2= ', k2_c)

        sensor1 = 'LS%d' %Landsat_nr
        res1 = '30m'
        res2 = '30m'
        res3 = '30m'

        # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
        # for the different Landsat images (L5, L7, or L8)
        ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
        ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
        ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])

        # Open one band - To get the metadata of the landsat images only once (to get the extend)
        src_FileName = os.path.join(input_folder, '%s_B2.TIF' %Name_Landsat_Image)  # before 10!
        ls, band_data, ulx, uly, lrx, lry, x_size_ls, y_size_ls = Get_Extend_Landsat(src_FileName)
        print('Original LANDSAT Image - ')
        print('  Size :', x_size_ls, y_size_ls)
        print('  Upper Left corner x, y: ', ulx, ', ', uly)
        print('  Lower right corner x, y: ', lrx, ', ', lry)

        lsc, ulx, uly, lrx, lry, epsg_to = SEBAL.reproject_dataset_example(src_FileName, Example_fileName)

        #	Get the extend of the remaining landsat file	after clipping based on the DEM file
        y_size_lsc = lsc.RasterYSize
        x_size_lsc = lsc.RasterXSize
        shape_lsc = [x_size_lsc, y_size_lsc]

        print('--- ')
        print('Cropped LANDSAT Image - ')
        print('  Size :', x_size_lsc, y_size_lsc)
        print('  Upper Left corner x, y: ', ulx, ', ',  uly)
        print('  Lower right corner x, y: ', lrx, ', ', lry)

        # if landsat 5 or 7 is used then first create a mask for removing the no data stripes
        if Landsat_nr == 5 or Landsat_nr == 7:
            src_FileName = os.path.join(input_folder, '%s_B6.TIF' % (Name_Landsat_Image)) #open smallest band
            if not os.path.exists(src_FileName):
                src_FileName = os.path.join(input_folder, '%s_B6_VCID_2.TIF' % (Name_Landsat_Image))
            src_FileName_2 = os.path.join(input_folder, '%s_B1.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_3 = os.path.join(input_folder, '%s_B3.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_4 = os.path.join(input_folder, '%s_B4.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_5 = os.path.join(input_folder, '%s_B7.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_6 = os.path.join(input_folder, '%s_B2.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_7 = os.path.join(input_folder, '%s_B5.TIF' % (Name_Landsat_Image)) #open smallest band
            ls_data=Open_landsat(src_FileName,Example_fileName)
            ls_data_2=Open_landsat(src_FileName_2,Example_fileName)
            ls_data_3=Open_landsat(src_FileName_3,Example_fileName)
            ls_data_4=Open_landsat(src_FileName_4,Example_fileName)
            ls_data_5=Open_landsat(src_FileName_5,Example_fileName)
            ls_data_6=Open_landsat(src_FileName_6,Example_fileName)
            ls_data_7=Open_landsat(src_FileName_7,Example_fileName)

            # create and save the landsat mask for all images based on band 11 (smallest map)
            QC_Map=np.zeros((shape_lsc[1], shape_lsc[0]))
            QC_Map=np.where(np.logical_or.reduce((ls_data==0,ls_data_2==0,ls_data_3==0,ls_data_4==0,ls_data_5==0,ls_data_6==0,ls_data_7==0)),1,0)

        # If landsat 8 then use landsat band 10 and 11
        elif Landsat_nr == 8:
             src_FileName_11 = os.path.join(input_folder, '%s_B11.TIF' % (Name_Landsat_Image)) #open smallest band
             ls_data_11=Open_landsat(src_FileName_11,Example_fileName)

             src_FileName_10 = os.path.join(input_folder, '%s_B10.TIF' % (Name_Landsat_Image)) #open smallest band
             ls_data_10=Open_landsat(src_FileName_10, Example_fileName)

             # create and save the landsat mask for all images based on band 10 and 11
             QC_Map=np.zeros((shape_lsc[1], shape_lsc[0]))
             QC_Map=np.where(np.logical_or(ls_data_11==0, ls_data_10==0),1,0)

        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        # Open data of the landsat mask
        ls_data=Open_landsat(src_FileName, Example_fileName)

        # Create 3D array to store Spectral radiance and Reflectivity for each band
        Reflect, Spec_Rad = Landsat_Reflect(Bands, input_folder, Name_Landsat_Image, output_folder, shape_lsc, QC_Map, Lmax, Lmin, ESUN_L5, ESUN_L7, ESUN_L8, cos_zn, dr, Landsat_nr, Example_fileName)

        # save spectral data
        for i in range(0,6):
            spec_ref_fileName = os.path.join(output_folder, 'Output_radiation_balance','%s_spectral_reflectance_B%s_%s_%s_%s.tif' %(sensor1, Bands[i], res3, year, DOY))
            SEBAL.save_GeoTiff_proy(lsc, Reflect[:, :, i], spec_ref_fileName, shape_lsc, nband=1)

    else:
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
        if (ws['B%d' % number].value) is not None:

            # Output folder NDVI
            ndvi_fileName_user = os.path.join(output_folder, 'Output_vegetation', 'User_NDVI_%s_%s_%s.tif' %(res3, year, DOY))
            NDVI=SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['B%d' % number].value),ndvi_fileName_user,Example_fileName)

            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
            water_mask_temp[NDVI < 0.0] = 1.0
            SEBAL.save_GeoTiff_proy(lsc, NDVI, ndvi_fileName_user, shape_lsc, nband=1)

        else:
            # use the Landsat reflectance to calculate the surface albede, NDVI
            NDVI = SEBAL.Calc_NDVI(Reflect)

            # Calculate temporal water mask
            water_mask_temp=SEBAL.Water_Mask(shape_lsc,Reflect)

    except:
        assert "Please check the NDVI input path"

    # Check Water Mask and replace if it is filled in the additianal data sheet
    try:
        if (ws['E%d' % number].value) is not None:

            # Overwrite the Water mask and change the output name
            water_mask_temp_fileName = os.path.join(output_folder, 'Output_soil_moisture', 'User_Water_mask_temporary_%s_%s_%s.tif' %(res2, year, DOY))
            water_mask_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['E%d' % number].value), water_mask_temp_fileName, Example_fileName)
            SEBAL.save_GeoTiff_proy(lsc, water_mask_temp, water_mask_temp_fileName, shape_lsc, nband=1)

    except:
        assert "Please check the Water Mask input path"

    # Check Surface albedo
    try:
        if (ws['C%d' % number].value) is not None:

            # Output folder surface albedo
            surface_albedo_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_albedo_%s_%s_%s.tif' %(res2, year, DOY))
            Surf_albedo=SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['C%d' % number].value),surface_albedo_fileName,Example_fileName)
            SEBAL.save_GeoTiff_proy(lsc, Surf_albedo, surface_albedo_fileName, shape_lsc, nband=1)

        else:

            # use the Landsat reflectance to calculate the surface albede, NDVI
            Surf_albedo = SEBAL.Calc_albedo(Reflect, path_radiance, Apparent_atmosf_transm)

    except:
          assert "Please check the Albedo input path"

    # calculate vegetation properties
    FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity=SEBAL.Calc_vegt_para(NDVI, water_mask_temp,shape_lsc)

    print('Average NDVI = %s' %np.nanmean(NDVI))
    print('Average Surface Albedo = %s' %np.nanmean(Surf_albedo))
    print('Average LAI = %s' %np.nanmean(LAI))
    print('Average Vegetation Cover = %s' %np.nanmean(vegt_cover))
    print('Average FPAR = %s' %np.nanmean(FPAR))

    return(Surf_albedo, NDVI, LAI, vegt_cover, FPAR, Nitrogen, tir_emis, b10_emissivity, water_mask_temp, QC_Map)

def Get_LS_Para_Thermal(workbook, number, Example_fileName, year, DOY,  water_mask_temp, b10_emissivity, Temp_inst, Rp, tau_sky, surf_temp_offset, Thermal_Sharpening_not_needed, DEM_fileName, UTM_Zone, eact_inst, QC_Map):

    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)
    Image_Type = 1

    # Open General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    # Open the Landsat_Input sheet
    ws = workbook['Landsat_Input']

    # Extract Landsat name, number and amount of thermal bands from excel file
    Name_Landsat_Image = str(ws['B%d' %number].value)
    Bands_thermal = int(ws['D%d' %number].value)
    Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)
                                                         # temperature: 1 = Band 6 for LS_5 & 7, Band 10 for LS_8 (optional)
    # Define bands used for each Landsat number
    if Landsat_nr == 5 or Landsat_nr == 7:
        Bands = np.array([1, 2, 3, 4, 5, 7, 6])
    elif Landsat_nr == 8:
        Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
    else:
        print('Landsat image not supported, use Landsat 7 or 8')

    # Open MTL landsat and get the correction parameters
    Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' %Name_Landsat_Image)
    Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)

    sensor1 = 'LS%d' %Landsat_nr
    sensor2 = 'LS%d' %Landsat_nr
    res1 = '30m'
    res2 = '30m'
    res3 = '30m'

    # Open the Landsat_Input sheet
    ws = workbook['Additional_Input']

    # If all additional fields are filled in than do not open the datasets
    if ws['D%d' % number].value is None:

        # Define bands used for each Landsat number
        if Landsat_nr == 5 or Landsat_nr == 7:
            Bands = np.array([1, 2, 3, 4, 5, 7, 6])
        elif Landsat_nr == 8:
           Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        print('...................... Open Landsat Thermal ........................')

        # Check if a surface temperature dataset is defined. If so use this one instead of the Landsat, otherwise Landsat
        therm_data = Landsat_therm_data(Bands, input_folder, Name_Landsat_Image, output_folder, shape_lsc, QC_Map, Example_fileName)

        # Create Cloud mask if BQA map is available (newer version Landsat images)
        BQA_LS_Available = 0
        if os.path.exists(os.path.join(input_folder, '%s_BQA.TIF' %Name_Landsat_Image)):
            src_FileName_BQA = os.path.join(input_folder, '%s_BQA.TIF' %Name_Landsat_Image)
            ls_data_BQA = Open_landsat(src_FileName_BQA, Example_fileName)
            if Landsat_nr == 8:
                Cloud_Treshold = 3000  #2720
            if Landsat_nr == 5 or Landsat_nr == 7:
                Cloud_Treshold = 700
            QC_mask_Cloud = np.copy(ls_data_BQA)
            QC_mask_Cloud[ls_data_BQA<Cloud_Treshold] = 0
            QC_mask_Cloud[ls_data_BQA>=Cloud_Treshold] = 1
            BQA_LS_Available = 1

        # Calculate surface temperature and create a cloud mask
        Surface_temp, cloud_mask_temp = SEBAL.Calc_surface_water_temp(Temp_inst, Landsat_nr, Lmax, Lmin, therm_data, b10_emissivity, k1_c, k2_c, eact_inst, shape_lsc, water_mask_temp, Bands_thermal, Rp, tau_sky, surf_temp_offset, Image_Type)

        # Replace clouds mask calculated by SEBAL by the official BQA file if this exists
        if BQA_LS_Available == 1:
            cloud_mask_temp = QC_mask_Cloud

        Surface_temp[cloud_mask_temp == 1] = np.nan
        print('Mean Surface Temperature = %s Kelvin' %np.nanmean(Surface_temp))

    else:
        try:
            # Output folder surface temperature
            surf_temp_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_temp_%s_%s_%s.tif' %(res2, year, DOY))
            Surface_temp = SEBAL.Reshape_Reproject_Input_data(r'%s' %str(ws['D%d' % number].value),surf_temp_fileName,Example_fileName)
            cloud_mask_temp = np.zeros([int(np.shape(Surface_temp)[1]),int(np.shape(Surface_temp)[0])])
            Thermal_Sharpening_not_needed = 0

        except:
            assert "Please check the surface temperature input path"

    return(Surface_temp, cloud_mask_temp, Thermal_Sharpening_not_needed)

#------------------------------------------------------------------------------
def info_band_metadata(filename, Bands):
    """
    This function retrieves Landsat band information (minimum and maximum
    radiance) from the metadata file.

    """
    Lmin = np.zeros(len(Bands))  # Minimum band radiance, for each band
    Lmax = np.zeros(len(Bands))  # Maximum band radiance, for each band
    k1_const = np.zeros(len(Bands)-6)  # TIRS_Thermal constant k1 ######
    k2_const = np.zeros(len(Bands)-6)  # TIRS_Thermal constant k2 ######
    for band in Bands:
        Landsat_meta = open(filename, "r")  # Open metadata file
        for line in Landsat_meta:
            if re.match("(.*)RADIANCE_MINIMUM_BAND_%1d(.*)" % band, line):
                words = line.split()
                value = float(words[2])
                Lmin[np.where(Bands == band)[0][0]] = value
            if re.match("(.*)RADIANCE_MAXIMUM_BAND_%1d(.*)" % band, line):
                words = line.split()
                value = float(words[2])
                Lmax[np.where(Bands == band)[0][0]] = value
            if re.match("(.*)K1_CONSTANT_BAND_%1d(.*)" % band, line):  # #####
                words = line.split()
                value = float(words[2])
                k1_const[np.where(Bands == band)[0][0]-6] = value
            if re.match("(.*)K2_CONSTANT_BAND_%1d(.*)" % band, line):  # #####
                words = line.split()
                value = float(words[2])
                k2_const[np.where(Bands == band)[0][0]-6] = value
    return Lmin, Lmax, k1_const, k2_const

#------------------------------------------------------------------------------
def Get_Extend_Landsat(src_FileName):
    """
    This function gets the extend of the landsat image
    """
    ls = gdal.Open(src_FileName)       # Open Landsat image
    geo_t_ls = ls.GetGeoTransform()    # Get the Geotransform vector
    x_size_ls = ls.RasterXSize         # Raster xsize - Columns
    y_size_ls = ls.RasterYSize         # Raster ysize - Rows
    (ulx, uly) = geo_t_ls[0], geo_t_ls[3]
    (lrx, lry) = (geo_t_ls[0] + geo_t_ls[1] * x_size_ls,
                  geo_t_ls[3] + geo_t_ls[5] * y_size_ls)
    band_data = ls.GetRasterBand(1)

    return(ls,band_data,ulx,uly,lrx,lry,x_size_ls,y_size_ls)

#------------------------------------------------------------------------------
def Landsat_Reflect(Bands,input_folder,Name_Landsat_Image,output_folder,shape_lsc,QC_Map,Lmax,Lmin,ESUN_L5,ESUN_L7,ESUN_L8,cos_zn,dr,Landsat_nr, proyDEM_fileName):
    """
    This function calculates and returns the reflectance and spectral radiation from the landsat image.
    """

    Spec_Rad = np.zeros((shape_lsc[1], shape_lsc[0], 7))
    Reflect = np.zeros((shape_lsc[1], shape_lsc[0], 7))
    for band in Bands[:-(len(Bands)-6)]:
        # Open original Landsat image for the band number
        src_FileName = os.path.join(input_folder, '%s_B%1d.TIF'
                                    % (Name_Landsat_Image, band))

        ls_data=Open_landsat(src_FileName, proyDEM_fileName)
        ls_data = ls_data * (-1 * QC_Map + 1)
        # stats = band_data.GetStatistics(0, 1)

        index = np.where(Bands[:-(len(Bands)-6)] == band)[0][0]
        if Landsat_nr == 8:
            # Spectral radiance for each band:
            L_lambda = Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L8, index, cos_zn, dr)
        elif Landsat_nr == 7:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L7, index, cos_zn, dr)
        elif Landsat_nr == 5:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda =Landsat_rho_lambda(L_lambda, ESUN_L5, index, cos_zn, dr)
        else:
            print('Landsat image not supported, use Landsat 5, 7 or 8')

        Spec_Rad[:, :, index] = L_lambda
        Reflect[:, :, index] = rho_lambda
    Reflect = Reflect.clip(0.0, 1.0)
    return(Reflect,Spec_Rad)


#------------------------------------------------------------------------------
def Landsat_L_lambda(Lmin,Lmax,ls_data,index,Landsat_nr):
    """
    Calculates the lambda from landsat
    """
    if Landsat_nr==8:
        L_lambda = ((Lmax[index] - Lmin[index]) / (65535 - 1) * ls_data + Lmin[index])
    elif Landsat_nr == 5 or Landsat_nr ==7:
        L_lambda = (Lmax[index] - Lmin[index]) / 255 * ls_data + Lmin[index]
    return(L_lambda)


#------------------------------------------------------------------------------
def Landsat_rho_lambda(L_lambda,ESUN,index,cos_zn,dr):
    """
    Calculates the rho from landsat
    """
    rho_lambda = np.pi * L_lambda / (ESUN[index] * cos_zn * dr)
    return(rho_lambda)


#------------------------------------------------------------------------------
def Landsat_therm_data(Bands,input_folder,Name_Landsat_Image,output_folder, shape_lsc, QC_Map, proyDEM_fileName):
    """
    This function calculates and returns the thermal data from the landsat image.
    """

    ClipLandsat = (-1 * QC_Map) + 1

    therm_data = np.zeros((shape_lsc[1], shape_lsc[0], len(Bands)-6))
    for band in Bands[-(len(Bands)-6):]:
        # Open original Landsat image for the band number
        src_FileName = os.path.join(input_folder, '%s_B%1d.TIF'
                                    % (Name_Landsat_Image, band))
        if not os.path.exists(src_FileName):
             src_FileName = os.path.join(input_folder, '%s_B%1d_VCID_2.TIF'
                                    % (Name_Landsat_Image, band))

        ls_data=Open_landsat(src_FileName, proyDEM_fileName)
        ls_data = ls_data*ClipLandsat
        index = np.where(Bands[:] == band)[0][0] - 6
        therm_data[:, :, index] = ls_data

    return(therm_data)

#------------------------------------------------------------------------------
def Open_landsat(src_FileName, proyDEM_fileName):
    """
    This function opens a landsat image and returns the data array of a specific landsat band.
    """
    import SEBAL.pySEBAL.pySEBAL_code as SEBAL

    # crop band to the DEM extent
    ls, ulx, uly, lrx, lry, epsg_to = SEBAL.reproject_dataset_example(src_FileName, proyDEM_fileName)

    # Open the cropped Landsat image for the band number
    ls_data = ls.GetRasterBand(1).ReadAsArray()
    return(ls_data)

#------------------------------------------------------------------------------
def info_general_metadata(filename):
    """
    This function retrieves general information of the Landsat image
    (date and time aquired, UTM zone, sun elevation) from the
    metadata file.

    """
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)SCENE_CENTER_TIME(.*)", line): # search in metadata for line SCENE_CENTER_TIME
            words = line.split()# make groups of words which are divided by an open space
            time_list = words[2].split(':', 2) # Take the second word of words and split the word which are divided by :
            if len(time_list[0])== 3:
                time_list[0]=time_list[0][1:3]
                time_list[2]=time_list[2][0:-1]
            hour = float(time_list[0]) # take the first word of time_list
            minutes = float(time_list[1]) + float(time_list[2][:-1]) / 60 # Take the second and third word of time_list and place :-1 to remove Z behined minutes
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)DATE_ACQUIRED(.*)", line):
            words = line.split()
            DOY = time.strptime(words[2], "%Y-%m-%d").tm_yday
            year = time.strptime(words[2], "%Y-%m-%d").tm_year
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)UTM_ZONE(.*)", line):
            words = line.split()
            UTM_Zone = int(words[2])
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)SUN_ELEVATION(.*)", line):
            words = line.split()
            Sun_elevation = float(words[2])

    return year, DOY, hour, minutes, UTM_Zone, Sun_elevation