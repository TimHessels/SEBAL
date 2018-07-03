# -*- coding: utf-8 -*-

"""
pySEBAL_3.3.8

@author: Tim Hessels, Jonna van Opstal, Patricia Trambauer, Wim Bastiaanssen,
         Mohamed Faouzi Smiej, Yasir Mohamed, and Ahmed Er-Raji
         UNESCO-IHE
         September 2017
"""
import os
import shutil
import numpy as np
from osgeo import osr
import gdal
from math import sin, cos, pi
import subprocess
from openpyxl import load_workbook
from pyproj import Proj, transform


def main(number, inputExcel, NDVI, Surf_albedo, Surface_temp, water_mask,
         proyDEM_fileName, QC_Map, year, DOY, hour, minutes, pixel_spacing,
         UTM_Zone, Cold_Pixel_Constant, Hot_Pixel_Constant):

    '''
    Additional input now from excel sheet:
       1. output_folder
       2. Image_Type (1,2,3: LS, VIIRS, MODIS)
       3. DEM_fileName (HydroSHED etc.)
       4. Temp_inst (Constant, GLDAS or CFSR)
       5. Temp_24 (Constant, GLDAS or CFSR)
       6. RH_inst (Constant, GLDAS or CFSR)
       7. RH_24 (Constant, GLDAS or CFSR)
       8. Wind_inst (Constant, GLDAS or CFSR)
       9. Wind_24 (Constant, GLDAS or CFSR)
       10. zx (Constant, GLDAS or CFSR)
       11. Rs_inst (Constant, GLDAS or CFSR)
       12. Rs_24 (Constant, GLDAS or CFSR)
       13. h_obst (constant or based on LU)
       14. Theta_sat_top (default value or soil map)
       15. Theta_sat_sub (default value or soil map)
       16. Theta_res_top (default value or soil map)
       17. Theta_res_sub (default value or soil map)
       18. Soil_moisture_wilting_point (default value or soil map)
       19. Depletion factor (constant or based on LU)
       20. Field_Capacity (default value or soil map)
       21. LUEmax (constant or based on LU)
    '''

    # Open Excel workbook to get the additional input mentioned above
    wb = load_workbook(inputExcel)

    # Open the General_Input sheet
    ws = wb['General_Input']

    # '---------------------------------------------------------'
    # '--------------------- Get Datasets ----------------------'
    # '---------------------------------------------------------'

    # ......................... General Data .......................... '

    # _________________________ 1 _______________________________
    # Extract the input and output folder, and Image type from the excel file
    output_folder = r"%s" % str(ws['C%d' % number].value)

    # Create or empty output folder if not exists
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    # _________________________ 2 _______________________________
    # Type of Image (1=Landsat & 2 = VIIRS & PROBA-V & 3 = MODIS)
    Image_Type = int(ws['D%d' % number].value)

    # _________________________ 3 _______________________________
    # Extract the Path to the DEM map from the excel file
    DEM_fileName = r"%s" % str(ws['E%d' % number].value)

    # ......................... Meteo Data ............................. '

    # Open the Meteo_Input sheet
    ws = wb['Meteo_Input']

    # _________________________ 4 _______________________________

    # ---------------------------- Instantaneous Air Temperature ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:  # Instantaneous Air Temperature (°C)
        Temp_inst = float(ws['B%d' % number].value)

        # If the data is a value than call this variable 0
        Temp_inst_kind_of_data = 0
        print 'Instantaneous Temperature constant value of \
        = %s (Celcius degrees)' % Temp_inst

    # if the data is not a value, than open as a string
    except:
        Temp_inst_name = '%s' % str(ws['B%d' % number].value)

        # If the data is a string than call this variable 1
        Temp_inst_kind_of_data = 1
        print 'Map to the Instantaneous Temperature = %s' % (Temp_inst_name)

    # _________________________ 5 _______________________________

    # ---------------------------- Daily Average Air Temperature ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:  # daily average Air Temperature (°C)
        Temp_24 = float(ws['C%d' % number].value)

        # If the data is a value than call this variable 0
        Temp_24_kind_of_data = 0
        print 'Daily average Temperature constant value of \
        = %s (Celcius degrees)' % Temp_24

    # if the data is not a value, than open as a string
    except:
        Temp_24_name = '%s' % str(ws['C%d' % number].value)

        # If the data is a string than call this variable 1
        Temp_24_kind_of_data = 1
        print 'Map to the Daily average Temperature = %s' % (Temp_24_name)

    # _________________________ 6 _______________________________

    # ---------------------------- Instantaneous Relative humidity ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:  # Instantaneous Relative humidity (%)
        RH_inst = float(ws['D%d' % number].value)

        # If the data is a value than call this variable 0
        RH_inst_kind_of_data = 0
        print 'Instantaneous Relative humidity constant value of \
        = %s (percentage)' % RH_inst

    # if the data is not a value, than open as a string
    except:
        RH_inst_name = '%s' % str(ws['D%d' % number].value)

        # If the data is a string than call this variable 1
        RH_inst_kind_of_data = 1
        print ('Map to the Instantaneous Relative humidity  = %s'
               % RH_inst_name)

    # _________________________ 7 _______________________________

    # ---------------------------- daily average Relative humidity ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:  # daily average Relative humidity (%)
        RH_24 = float(ws['E%d' % number].value)

        # If the data is a value than call this variable 0
        RH_24_kind_of_data = 0
        print ('Daily average Relative humidity constant value of \
        = %s (percentage)' % RH_24)

    # if the data is not a value, than open as a string
    except:
        RH_24_name = '%s' % str(ws['E%d' % number].value)

        # If the data is a string than call this variable 1
        RH_24_kind_of_data = 1
        print 'Map to the Daily average Relative humidity = %s' % RH_24_name

    # _________________________ 8 _______________________________

    # ---------------------------- instantaneous Wind Speed ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:  # instantaneous Wind Speed (m/s)
        Wind_inst = float(ws['G%d' % number].value)

        # If the data is a value than call this variable 0
        Wind_inst_kind_of_data = 0
        print ('Instantaneous Wind Speed constant value of = %s (m/s)'
               % Wind_inst)

    # if the data is not a value, than open as a string
    except:
        Wind_inst_name = '%s' % str(ws['G%d' % number].value)

        # If the data is a string than call this variable 1
        Wind_inst_kind_of_data = 1
        print 'Map to the Instantaneous Wind Speed = %s' % Wind_inst_name

    # _________________________ 9 _______________________________

    # ---------------------------- daily Wind Speed ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    try:
        Wind_24 = float(ws['H%d' % number].value)  # daily Wind Speed (m/s)

        # If the data is a value than call this variable 0
        Wind_24_kind_of_data = 0
        print 'Daily Wind Speed constant value of = %s (m/s)' % Wind_24

    # if the data is not a value, than open as a string
    except:
        Wind_24_name = '%s' % str(ws['H%d' % number].value)

        # If the data is a string than call this variable 1
        Wind_24_kind_of_data = 1
        print 'Map to the Daily Wind Speed = %s' % Wind_24_name

    # _________________________ 10 _______________________________

    # Height of the wind speed measurement
    zx = float(ws['F%d' % number].value)  # Height of wind speed measurement
    print 'Height at which wind speed is measured = %s (m)' % zx

    # _________________________ 11 _______________________________

    # Define the method of radiation (1 or 2)
    # 1=Transm_24 will be calculated Rs_24 must be given
    # 2=Rs_24 will be determined Transm_24 must be given
    Method_Radiation_24 = int(ws['I%d' % number].value)
    print ('Method for daily radiation (1=Rs_24, 2=Transm_24) = %s'
           % Method_Radiation_24)

    # if method radiation is 1
    # ---------------------------- daily Surface Solar Radiation ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    if Method_Radiation_24 == 1:
        try:
            # daily Surface Solar Radiation (W/m2)
            # only required when Method_Radiation_24 = 1
            Rs_24 = float(ws['J%d' % number].value)

            # If the data is a value than call this variable 0
            Rs_24_kind_of_data = 0
            print 'Daily Surface Solar Radiation constant value of = %s \
            (W/m2)' % Rs_24

        # If the data is not a value, than open as a string
        except:
            Rs_24_name = '%s' % str(ws['J%d' % number].value)

            #  If the data is a string than call this variable 1
            Rs_24_kind_of_data = 1
            print 'Map to the Daily Surface Solar Radiation = %s' % Rs_24_name

    # if method radiation is 2
    # ---------------------------- daily transmissivity ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    if Method_Radiation_24 == 2:
        try:
            # daily transmissivity, Typical values between 0.65 and 0.8
            # only required when Method_Radiation_24 = 2
            Transm_24 = float(ws['K%d' % number].value)

            # If the data is a value than call this variable 0
            Transm_24_kind_of_data = 0
            print 'Daily transmissivity constant value of = %s' % Transm_24

        # if the data is not a value, than open as a string
        except:
            Transm_24_name = '%s' % str(ws['K%d' % number].value)

            # If the data is a string than call this variable 1
            Transm_24_kind_of_data = 1
            print 'Map to the Daily transmissivity = %s' % Transm_24_name

    # _________________________ 12 _______________________________

    # Define the method of instataneous radiation (1 or 2)
    #  1=Transm_inst will be calculated Rs_inst must be given
    #  2=Rs_24 will be determined Transm_24 must be given
    Method_Radiation_inst = int(ws['L%d' % number].value)
    print 'Method for instantaneous radiation (1=Rs_inst, 2=Transm_inst) \
    = %s' % Method_Radiation_inst

    # if method instantaneous radiation is 1
    # --------------------- Instantaneous Surface Solar Radiation ------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    if Method_Radiation_inst == 1:
        try:
            # Instantaneous Surface Solar Radiation (W/m2)
            # only required when Method_Radiation_inst = 1
            Rs_in_inst = float(ws['M%d' % number].value)

            # If the data is a value than call this variable 0
            Rs_in_inst_kind_of_data = 0
            print 'Instantaneous Surface Solar Radiation constant value of\
            = %s (W/m2)' % Rs_in_inst

        # if the data is not a value, than open as a string
        except:
            Rs_in_inst_name = '%s' % str(ws['M%d' % number].value)

            # If the data is a string than call this variable 1
            Rs_in_inst_kind_of_data = 1
            print ('Map to the Instantaneous Surface Solar Radiation = %s'
                   % Rs_in_inst_name)

    # if method instantaneous radiation is 2
    # ---------------------------- Instantaneous transmissivity------------
    # Open meteo data, first try to open as value, otherwise as string (path)
    if Method_Radiation_inst == 2:
        try:
            # Instantaneous transmissivity,Typical values between 0.7 and 0.85
            # only required when Method_Radiation_inst = 2
            Transm_inst = float(ws['N%d' % number].value)

            # If the data is a value than call this variable 0
            Transm_inst_kind_of_data = 0

        # if the data is not a value, than open as a string
        except:
            Transm_inst_name = '%s' % str(ws['N%d' % number].value)

            # If the data is a string than call this variable 1
            Transm_inst_kind_of_data = 1

    radiation_inst_fileName = os.path.join(output_folder,
                                           'Output_radiation_balance',
                                           'Ra_inst_%s_%s.tif' % (year, DOY))

    # ......................... Obstacle height ..........................

    # _________________________ 13 _______________________________

    # Open obstacle height data, first try to open as value,
    # otherwise as string (path)
    try:
        h_obst = float(ws['O%d' % number].value)  # Obstacle height
        h_obst_kind_of_data = 0
        print 'Obstacle height constant value of = %s (Meter)' % h_obst
    except:
        h_obst_name = '%s' % str(ws['O%d' % number].value)
        h_obst_kind_of_data = 1  # Obstacle height (m)
        print 'Map to the Obstacle height = %s' % h_obst_name

    # ......................... Soil Properties  ..........................

    # Open soil input sheet
    ws = wb['Soil_Input']

    # _________________________ 14 _______________________________

    # ---------------------------- Saturated soil moisture content topsoil  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Theta_sat_top = float(ws['B%d' %number].value)                # Saturated soil moisture content topsoil

        # If the data is a value than call this variable 0
        Theta_sat_top_kind_of_data = 0
        print 'Saturated soil moisture content topsoil constant value of = %s (cm3/cm3)' %(Theta_sat_top)

    # if the data is not a value, than open as a string
    except:
        Theta_sat_top_name = '%s' %str(ws['B%d' %number].value)

        # If the data is a value than call this variable 1
        Theta_sat_top_kind_of_data = 1
        print 'Map to the Saturated soil moisture content topsoil = %s' %(Theta_sat_top_name)

    # _________________________ 15 _______________________________

    # ---------------------------- Saturated soil moisture content subsoil  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Theta_sat_sub = float(ws['C%d' %number].value)                # Saturated soil moisture content subsoil

        # If the data is a value than call this variable 0
        Theta_sat_sub_kind_of_data = 0
        print 'Saturated soil moisture content subsoil constant value of = %s (cm3/cm3)' %(Theta_sat_sub)

    # if the data is not a value, than open as a string
    except:
        Theta_sat_sub_name = '%s' %str(ws['C%d' %number].value)

        # If the data is a value than call this variable 1
        Theta_sat_sub_kind_of_data = 1
        print 'Map to the Saturated soil moisture content subsoil  = %s' %(Theta_sat_sub_name)

    # _________________________ 16 _______________________________

    # ---------------------------- Residual soil moisture content topsoil  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Theta_res_top = float(ws['D%d' %number].value)              # Residual soil moisture content

        # If the data is a value than call this variable 0
        Theta_res_top_kind_of_data = 0
        print 'Residual soil moisture content topsoil constant value of = %s (cm3/cm3)' %(Theta_res_top)

    # if the data is not a value, than open as a string
    except:
        Theta_res_top_name = '%s' %str(ws['D%d' %number].value)

        # If the data is a value than call this variable 1
        Theta_res_top_kind_of_data = 1
        print 'Map to the Residual soil moisture content topsoil = %s' %(Theta_res_top_name)

    # _________________________ 17 _______________________________

    # ---------------------------- Residual soil moisture content subsoil  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Theta_res_sub = float(ws['E%d' %number].value)              # Residual soil moisture content subsoil

        # If the data is a value than call this variable 0
        Theta_res_sub_kind_of_data = 0
        print 'Residual soil moisture content subsoil constant value of = %s (cm3/cm3)' %(Theta_res_sub)

    # if the data is not a value, than open as a string
    except:
        Theta_res_sub_name = '%s' %str(ws['E%d' %number].value)

        # If the data is a value than call this variable 1
        Theta_res_sub_kind_of_data = 1
        print 'Map to the residual soil moisture content subsoil = %s' %(Theta_res_sub_name)

    # _________________________ 18 _______________________________

    # ---------------------------- Wilting point  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Soil_moisture_wilting_point = float(ws['G%d' %number].value)              # Wilting point

        # If the data is a value than call this variable 0
        Soil_moisture_wilting_point_kind_of_data = 0
        print 'Soil moisture wilting point constant value of = %s (cm3/cm3)' %(Soil_moisture_wilting_point)

    # if the data is not a value, than open as a string
    except:
        Soil_moisture_wilting_point_name = '%s' %str(ws['G%d' %number].value)

        # If the data is a value than call this variable 1
        Soil_moisture_wilting_point_kind_of_data = 1
        print 'Map to the soil moisture wilting point = %s' %(Soil_moisture_wilting_point_name)

    # _________________________ 19 _______________________________

    # ---------------------------- Depletion factor  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        depl_factor = float(ws['H%d' %number].value)              # Depletion factor

        # If the data is a value than call this variable 0
        depl_factor_kind_of_data = 0
        print 'Depletion factor constant value of = %s' %(depl_factor)

    # if the data is not a value, than open as a string
    except:
        depl_factor_name = '%s' %str(ws['H%d' %number].value)

        # If the data is a value than call this variable 1
        depl_factor_kind_of_data = 1
        print 'Map to the Depletion factor = %s' %(depl_factor_name)

    # _________________________ 20 _______________________________

    # ---------------------------- Fraction field capacity  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        Field_Capacity = float(ws['F%d' %number].value)              # Field capacity divided by saturation capacity

        # If the data is a value than call this variable 0
        Field_Capacity_kind_of_data = 0
        print 'Fraction field capacity and saturation constant value= %s' %(Field_Capacity)

    # if the data is not a value, than open as a string
    except:
        Field_Capacity_name = '%s' %str(ws['F%d' %number].value)

        # If the data is a value than call this variable 1
        Field_Capacity_kind_of_data = 1
        print 'Map to the Fraction field capacity and saturation = %s' %(Field_Capacity_name)

    print 'General constants for Module 13:'

    # _________________________ 21 _______________________________

    # ---------------------------- Light Use Efficiency  ------------
    # Open soil data, first try to open as value, otherwise as string (path)
    try:
        LUEmax = float(ws['I%d' %number].value)              # Field capacity divided by saturation capacity

        # If the data is a value than call this variable 0
        LUEmax_kind_of_data = 0
        print 'Light use efficiency constant value = %s' %(LUEmax)

    # if the data is not a value, than open as a string
    except:
        LUEmax_name = '%s' %str(ws['I%d' %number].value)

        # If the data is a value than call this variable 1
        LUEmax_kind_of_data = 1
        print 'Map to the Light use efficiency = %s' %(LUEmax_name)

    # User parameters
    surf_roughness_equation_used = 1 # NDVI model = 1, Raupach model = 2




    print '---------------------------------------------------------'
    print '-------------------- Open DEM ---------------------------'
    print '---------------------------------------------------------'

    # Open DEM and create Latitude and longitude files
    lat, lon, lat_fileName, lon_fileName = DEM_lat_lon(DEM_fileName, output_folder)

    # Reproject from Geog Coord Syst to UTM -
    # 1) DEM - Original DEM coordinates is Geographic: lat, lon
    lsc, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(
                DEM_fileName, pixel_spacing, UTM_Zone = UTM_Zone)
    band = lsc.GetRasterBand(1)   # Get the reprojected dem band
    ncol = lsc.RasterXSize        # Get the reprojected dem column size
    nrow = lsc.RasterYSize        # Get the reprojected dem row size
    shape = [ncol, nrow]

    # Read out the DEM band and print the DEM properties
    DEM_resh = band.ReadAsArray(0, 0, ncol, nrow)
    #data_DEM[data_DEM<0] = 1
    print 'Projected DEM - '
    print '   Size: ', ncol, nrow
    print '   Upper Left corner x, y: ', ulx_dem, ',', uly_dem
    print '   Lower right corner x, y: ', lrx_dem, ',', lry_dem

    # 2) Latitude File - reprojection
    # Define output name of the latitude file

    # reproject latitude to the landsat projection	 and save as tiff file
    lat_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(
               lat_fileName, pixel_spacing, UTM_Zone=UTM_Zone)

    # Get the reprojected latitude data
    lat_proy = lat_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)

    # 3) Longitude file - reprojection

    # reproject longitude to the landsat projection	 and save as tiff file
    lon_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(lon_fileName, pixel_spacing, UTM_Zone)

    # Get the reprojected longitude data
    lon_proy = lon_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)

    # Calculate slope and aspect from the reprojected DEM
    deg2rad, rad2deg, slope, aspect = Calc_Gradient(DEM_resh, pixel_spacing)

    print '---------------------------------------------------------'
    print '-------------------- Radiation --------------------------'
    print '---------------------------------------------------------'

    # Calculation of extraterrestrial solar radiation for slope and aspect
    Ra_mountain_24, Ra_inst, cos_zn, dr, phi, delta = Calc_Ra_Mountain(lon, DOY, hour, minutes, lon_proy, lat_proy, slope, aspect)
    save_GeoTiff_proy(lsc, Ra_inst, radiation_inst_fileName, shape, nband = 1 )

    if Image_Type == 2 or Image_Type == 3:
        Sun_elevation = 90 - (np.nanmean(cos_zn) * 180/np.pi)

    # Resample DEM related maps to resolution of clipped Landsat images - 30 m

    # 2) DEM
    lsc = gdal.Open(proyDEM_fileName)

    #	Get the extend of the remaining landsat file	after clipping based on the DEM file
    y_size_lsc = lsc.RasterYSize
    x_size_lsc = lsc.RasterXSize
    shape_lsc = [x_size_lsc, y_size_lsc]

    # 4) Reshaped instantaneous radiation
    Ra_inst = Reshape_Reproject_Input_data(radiation_inst_fileName, proyDEM_fileName)

    # 6a) Instantaneous Temperature
    if Temp_inst_kind_of_data is 1:
        try:
            Temp_inst = Reshape_Reproject_Input_data(Temp_inst_name, proyDEM_fileName)
        except:
            print 'ERROR: Check the instantenious Temperature input path in the meteo excel tab'

    # 6b) Daily Temperature
    if Temp_24_kind_of_data is 1:
        try:
             Temp_24 = Reshape_Reproject_Input_data(Temp_24_name, proyDEM_fileName)
        except:
            print 'ERROR: Check the daily Temperature input path in the meteo excel tab'

    # 6c) Daily Relative Humidity
    if RH_24_kind_of_data is 1:
        try:
            RH_24 = Reshape_Reproject_Input_data(RH_24_name, proyDEM_fileName)
        except:
            print 'ERROR: Check the instantenious Relative Humidity input path in the meteo excel tab'

     # 6d) Instantaneous Relative Humidity
    if RH_inst_kind_of_data is 1:
        try:
            RH_inst = Reshape_Reproject_Input_data(RH_inst_name, proyDEM_fileName)
        except:
            print 'ERROR: Check the daily Relative Humidity input path in the meteo excel tab'

     # 6e) Daily wind speed
    if Wind_24_kind_of_data is 1:
        try:
            Wind_24 = Reshape_Reproject_Input_data(Wind_24_name, proyDEM_fileName)
            Wind_24[Wind_24 < 1.5] = 1.5
        except:
            print 'ERROR: Check the daily wind input path in the meteo excel tab'

     # 6f) Instantaneous wind speed
    if Wind_inst_kind_of_data is 1:
        try:
            Wind_inst = Reshape_Reproject_Input_data(Wind_inst_name, proyDEM_fileName)
            Wind_inst[Wind_inst < 1.5] = 1.5
        except:
            print 'ERROR: Check the instantenious wind input path in the meteo excel tab'

    # 6g) Daily incoming Radiation
    if Method_Radiation_24 == 1:
        if Rs_24_kind_of_data is 1:
            try:
                Rs_24 = Reshape_Reproject_Input_data(Rs_24_name, proyDEM_fileName)
            except:
                print 'ERROR: Check the daily net radiation input path in the meteo excel tab'

    # 6h) Instantaneous incoming Radiation
    if Method_Radiation_inst == 1:
        if Rs_in_inst_kind_of_data is 1:
            try:
                Rs_in_inst = Reshape_Reproject_Input_data(Rs_in_inst_name, proyDEM_fileName)
            except:
                print 'ERROR: Check the instanenious net radiation input path in the meteo excel tab'

    # 6i) Daily Transmissivity
    if Method_Radiation_24 == 2:
        if Transm_24_kind_of_data is 1:
            try:
                Transm_24 = Reshape_Reproject_Input_data(Transm_24_name, proyDEM_fileName)
            except:
                print 'ERROR: Check the daily transmissivity input path in the meteo excel tab'

    # 6j) Instantaneous Transmissivity
    if Method_Radiation_inst == 2:
        if Transm_inst_kind_of_data is 1:
            try:
                 Transm_inst = Reshape_Reproject_Input_data(Transm_inst_name, proyDEM_fileName)
            except:
                print 'ERROR: Check the instantenious transmissivity input path in the meteo excel tab'

    # 6k) Theta saturated topsoil
    if Theta_sat_top_kind_of_data is 1:
        try:
           Theta_sat_top = Reshape_Reproject_Input_data(Theta_sat_top_name, proyDEM_fileName)
        except:
           print 'ERROR: Check the saturated top soil input path in the soil excel tab'

    # 6l) Theta saturated subsoil
    if Theta_sat_sub_kind_of_data is 1:
        try:
            Theta_sat_sub  =Reshape_Reproject_Input_data(Theta_sat_sub_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the saturated sub soil input path in the soil excel tab'

    # 6m) Theta residual topsoil
    if Theta_res_top_kind_of_data is 1:
        try:
            Theta_res_top=Reshape_Reproject_Input_data(Theta_res_top_name, proyDEM_fileName)
        except:
            print 'ERROR: Check the residual top soil input path in the soil excel tab'

    # 6n) Theta residual subsoil
    if Theta_res_sub_kind_of_data is 1:
        try:
            Theta_res_sub=Reshape_Reproject_Input_data(Theta_res_sub_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the residual sub soil input path in the soil excel tab'

    # 6o) Wilting point
    if Soil_moisture_wilting_point_kind_of_data is 1:
        try:
            Soil_moisture_wilting_point=Reshape_Reproject_Input_data(Soil_moisture_wilting_point_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the wilting point input path in the soil excel tab'

    # 6p) Fraction field capacity
    if Field_Capacity_kind_of_data is 1:
        try:
            Field_Capacity=Reshape_Reproject_Input_data(Field_Capacity_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the field capacity input path in the soil excel tab'

    # 6q) Light Use Efficiency
    if LUEmax_kind_of_data is 1:
        try:
            LUEmax=Reshape_Reproject_Input_data(LUEmax_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the LUE input path in the soil excel tab'

    # 6r) Obstacle height
    if h_obst_kind_of_data is 1:
        try:
            h_obst=Reshape_Reproject_Input_data(h_obst_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the obstacle height input path in the soil excel tab'

    # 6s) deplection factor
    if depl_factor_kind_of_data is 1:
        try:
            depl_factor=Reshape_Reproject_Input_data(depl_factor_name,proyDEM_fileName)
        except:
            print 'ERROR: Check the depletion factor input path in the soil excel tab'

    print '---------------------------------------------------------'
    print '-------------------- Meteo part 1 -----------------------'
    print '---------------------------------------------------------'

    # Computation of some vegetation properties
    # 1)
    #constants:
    Temp_lapse_rate = 0.0065  # Temperature lapse rate (°K/m)
    SB_const = 5.6703E-8  # Stefan-Bolzmann constant (watt/m2/°K4)

    # Atmospheric pressure for altitude:
    Pair = 101.3 * np.power((293 - Temp_lapse_rate * DEM_resh) / 293, 5.26)

    # Psychrometric constant (kPa / °C), FAO 56, eq 8.:
    Psychro_c = 0.665E-3 * Pair

    # Saturation Vapor Pressure at the air temperature (kPa):
    esat_inst = 0.6108 * np.exp(17.27 * Temp_inst / (Temp_inst + 237.3))
    esat_24 = 0.6108 * np.exp(17.27 * Temp_24 / (Temp_24 + 237.3))

    # Actual vapour pressure (kPa), FAO 56, eq 19.:
    eact_inst = RH_inst * esat_inst / 100
    eact_24 = RH_24 * esat_24 / 100
    print 'Instantaneous Saturation Vapor Pressure = ', '%0.3f (kPa)' % np.nanmean(esat_inst)
    print 'Instantaneous Actual vapour pressure =  ', '%0.3f (kPa)' % np.nanmean(eact_inst)
    print 'Daily Saturation Vapor Pressure = ', '%0.3f (kPa)' % np.nanmean(esat_24)
    print 'Daily Actual vapour pressure =  ', '%0.3f (kPa)' % np.nanmean(eact_24)


    print '---------------------------------------------------------'
    print '----------------------- Meteo ---------------------------'
    print '---------------------------------------------------------'

    # calculate vegetation properties
    FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity=Calc_vegt_para(NDVI, water_mask, shape_lsc)

    # Slope of satur vapour pressure curve at air temp (kPa / °C)
    sl_es_24 = 4098 * esat_24 / np.power(Temp_24 + 237.3, 2)

    # calculate the daily radiation or daily transmissivity or daily surface radiation based on the method defined by the user
    if Method_Radiation_24==1:
       Transm_24 = Rs_24/Ra_mountain_24

    if Method_Radiation_24==2:
       Rs_24 = Ra_mountain_24 * Transm_24

    # If method of instantaneous radiation 1 is used than calculate the Transmissivity
    if Method_Radiation_inst==1:
        Transm_corr=Rs_in_inst/Ra_inst

    # If method of instantaneous radiation 2 is used than calculate the instantaneous incomming Radiation
    if Method_Radiation_inst==2:
        # calculate the transmissivity index for direct beam radiation
        Transm_corr = Transm_inst + 2e-5 * DEM_resh
        # Instantaneous incoming short wave radiation (W/m2):
        Rs_in_inst = Ra_inst * Transm_corr

    # Atmospheric emissivity, by Bastiaanssen (1995):
    Transm_corr[Transm_corr<0.001]=0.1
    Transm_corr[Transm_corr>1]=1
    atmos_emis = 0.85 * np.power(-np.log(Transm_corr), 0.09)

    # Instantaneous incoming longwave radiation:
    lw_in_inst = atmos_emis * SB_const * np.power(Temp_inst + 273.15, 4)

    # calculates the ground heat flux and the solar radiation
    Rn_24,rn_inst,g_inst,Rnl_24_FAO =Calc_Meteo(Rs_24,eact_24,Temp_24,Surf_albedo,dr,tir_emis,Surface_temp, water_mask,NDVI,Transm_24,SB_const,lw_in_inst,Rs_in_inst)


    print '---------------------------------------------------------'
    print '-------------------- Hot/Cold Pixels --------------------'
    print '---------------------------------------------------------'

    # Temperature at sea level corrected for elevation: ??
    ts_dem,air_dens,Temp_corr=Correct_Surface_Temp(Surface_temp,Temp_lapse_rate,DEM_resh,Pair,dr,Transm_corr,cos_zn,Sun_elevation,deg2rad,QC_Map)

    # Selection of hot and cold pixels


    # Cold pixels vegetation
    ts_dem_cold_veg = Calc_Cold_Pixels_Veg(NDVI,QC_Map,ts_dem, Cold_Pixel_Constant)

    # Cold pixels water
    ts_dem_cold,cold_pixels,ts_dem_cold_mean = Calc_Cold_Pixels(ts_dem,water_mask,QC_Map,ts_dem_cold_veg,Cold_Pixel_Constant)
    if np.isnan(ts_dem_cold) == True:
        ts_dem_cold = Temp_inst

    # Hot pixels
    ts_dem_hot,hot_pixels = Calc_Hot_Pixels(ts_dem,QC_Map, water_mask, NDVI,Hot_Pixel_Constant, ts_dem_cold)

    print '---------------------------------------------------------'
    print '----------------- Sensible heat flux --------------------'
    print '---------------------------------------------------------'

    # Change the minimum windspeed to prevent high values in further calculations
    if Wind_inst_kind_of_data is 0:
        if Wind_inst<1.5:
            Wind_inst=1.5
    if Wind_24_kind_of_data is 0:
        if Wind_24<1.5:
            Wind_24=1.5

    # calculate windspeed at the blending height and the friction velocity by using the Raupach model or NDVI
    Surf_roughness,u_200,ustar_1=Calc_Wind_Speed_Friction(h_obst,Wind_inst,zx,LAI,NDVI,Surf_albedo,water_mask,surf_roughness_equation_used)

    # Computation of surface roughness for momentum transport
    k_vk = 0.41      # Von Karman constant

    # Sensible heat 1 (Step 5)
    # Corrected value for the aerodynamic resistance (eq 41 with psi2 = psi1):
    rah1 = np.log(2.0/0.01) / (k_vk * ustar_1)
    i=0
    L, psi_m200_stable, psi, psi_m200,h_inst,dT, slope_dt, offset_dt = sensible_heat(
            rah1, ustar_1, rn_inst, g_inst, ts_dem, ts_dem_hot, ts_dem_cold,
            air_dens, Surface_temp, k_vk,QC_Map, hot_pixels, slope)

    # do the calculation iteratively 10 times
    for i in range(1,10):
        L,psi,psi_m200,psi_m200_stable,h_inst,ustar_corr,rah_corr,dT, slope_dt, offset_dt = Iterate_Friction_Velocity(k_vk,u_200,Surf_roughness,g_inst,rn_inst, ts_dem, ts_dem_hot, ts_dem_cold,air_dens, Surface_temp,L,psi,psi_m200,psi_m200_stable,QC_Map, hot_pixels, slope)

    print '---------------------------------------------------------'
    print '-------------------- Evaporation ------------------------'
    print '---------------------------------------------------------'

    # calculate reference net radiation
    Rn_ref, Refl_rad_water, rah_grass=Calc_Rn_Ref(shape_lsc,water_mask,Rn_24,Ra_mountain_24,Transm_24,Rnl_24_FAO,Wind_24)

    # Calculate rah of PM for the ET act (dT after iteration) and ETpot (4 degrees)
    rah_pm_act=((np.log((2.0-0.0)/(Surf_roughness*0.1))*np.log((2.0-0.0)/(Surf_roughness)))/(k_vk*1.5**2))*((1-5*(-9.82*dT*(2.0-0.0))/((273.15+Temp_inst)*1.5**2))**(-0.75))
    rah_pm_act[rah_pm_act<25]=25

    rah_pm_pot=((np.log((2.0-0.0)/(Surf_roughness*0.1))*np.log((2.0-0.0)/(Surf_roughness)))/(k_vk*1.5**2))*((1-5*(-9.82*4.0*(2.0-0.0))/((273.15+Temp_inst)*1.5**2))**(-0.75))
    rah_pm_pot[rah_pm_pot<25]=25

    # calculate reference potential evaporation.
    ETpot_24,ETref_24,Lhv,rs_min=Calc_Ref_Pot_ET(LAI,Surface_temp,sl_es_24,Rn_ref,air_dens,esat_24,eact_24,rah_grass,Psychro_c,Rn_24,Refl_rad_water,rah_pm_pot)

    # Instantaneous evapotranspiration
    LE_inst = rn_inst - g_inst - h_inst

    # Evaporative fraction
    EF_inst=Calc_instantaneous_ET_fraction(LE_inst,rn_inst,g_inst)

    # Daily Evaporation and advection factor
    ETA_24, AF=Calc_ETact(esat_24,eact_24,EF_inst,Rn_24,Refl_rad_water,Lhv, Image_Type)
    # Bulk surface resistance (s/m):
    bulk_surf_resis_24=Calc_Bulk_surface_resistance(sl_es_24,Rn_24,Refl_rad_water,air_dens,esat_24,eact_24,rah_pm_act,ETA_24,Lhv,Psychro_c)

    # crop factor
    kc = ETA_24 / ETref_24  # Crop factor
    ETP_24 = np.where(ETpot_24 < ETA_24, ETA_24, ETpot_24)
    ET_24_deficit = ETP_24 - ETA_24
    kc_max = ETP_24 / ETref_24

    print '---------------------------------------------------------'
    print '-------------------- Soil Moisture ----------------------'
    print '---------------------------------------------------------'

    #  Calculate soil properties
    #SM_stress_trigger, total_soil_moisture, RZ_SM,moisture_stress_biomass,irrigation_needs,top_soil_moisture=Calc_Soil_Moisture(ETA_24,accum_prec_14d,accum_ETo_14d,EF_inst,water_mask,vegt_cover,Theta_sat,Theta_res)
    SM_stress_trigger, total_soil_moisture, root_zone_moisture_first, moisture_stress_biomass_first,top_soil_moisture,RZ_SM_NAN = Calc_Soil_Moisture(ETA_24,EF_inst,QC_Map,water_mask,vegt_cover,Theta_sat_top,Theta_sat_sub, Theta_res_top,Theta_res_sub, depl_factor,Field_Capacity,FPAR, Soil_moisture_wilting_point)

    # seperation of E and T
    Eact_24,Tpot_24,Tact_24,moisture_stress_biomass,T24_deficit,beneficial_fraction,root_zone_moisture_final,top_zone_moisture_final=Separate_E_T(LAI,ETP_24,Theta_res_top, Theta_res_sub,Theta_sat_top,Theta_sat_sub,top_soil_moisture,sl_es_24, Psychro_c,moisture_stress_biomass_first,vegt_cover,ETA_24,SM_stress_trigger,root_zone_moisture_first,total_soil_moisture)

    # Irrigation:
    irrigation_needs = Classify_Irrigation(moisture_stress_biomass, vegt_cover)


    print '---------------------------------------------------------'
    print '---------------------- Biomass --------------------------'
    print '---------------------------------------------------------'

    # calculate biomass production
    LUE,Biomass_prod,Biomass_wp,Biomass_deficit = Calc_Biomass_production(LAI,ETP_24,moisture_stress_biomass,ETA_24,Ra_mountain_24,Transm_24,FPAR,esat_24,eact_24,Temp_24,LUEmax)
    lsc=None


    print '...................................................................'
    print '............................DONE!..................................'
    print '...................................................................'

    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # FUNCTIONS
    #-------------------------------------------------------------------------

def Create_Buffer(Data_In):

   '''
   This function creates a 3D array which is used to apply the moving window
   '''
   Buffer_area = 2 # A block of 2 times Buffer_area + 1 will be 1 if there is the pixel in the middle is 1
   Data_Out=np.empty((len(Data_In),len(Data_In[1])))
   Data_Out[:,:] = Data_In
   for ypixel in range(1,Buffer_area + 1):
        Data_Out[:,0:-ypixel] += Data_In[:,ypixel:]
        Data_Out[:,ypixel:] += Data_In[:,:-ypixel]

        for xpixel in range(1,Buffer_area + 1):
           Data_Out[0:-xpixel,ypixel:] += Data_In[xpixel:,:-ypixel]
           Data_Out[xpixel:,ypixel:] += Data_In[:-xpixel,:-ypixel]
           Data_Out[0:-xpixel,0:-ypixel] += Data_In[xpixel:,ypixel:]
           Data_Out[xpixel:,0:-ypixel] += Data_In[:-xpixel,ypixel:]

           if ypixel==1:
                Data_Out[xpixel:,:] += Data_In[:-xpixel,:]
                Data_Out[:,0:-xpixel] += Data_In[:,xpixel:]


   Data_Out[Data_Out>0.1] = 1
   Data_Out[Data_Out<=0.1] = 0

   return(Data_Out)

def Calc_Biomass_production(LAI,ETP_24,moisture_stress_biomass,ETA_24,Ra_mountain_24,Transm_24,FPAR,esat_24,eact_24,Temp_24,LUEmax):
    """
    Function to calculate the biomass production and water productivity
    """

    # Biomass constants
    Th = 35.0                        # Upper limit of stomatal activity
    Kt = 23.0                        # Optimum conductance temperature (°C), range: [17 - 19]
    Tl = 0.0                         # Lower limit of stomatal activity

    Ksolar = Ra_mountain_24 * Transm_24

    # Incident Photosynthetically active radiation (PAR, MJ/m2) per time period
    PAR = 0.48 * Ksolar

    # Aborbed Photosynthetical Active Radiation (APAR) by the vegetation:
    APAR = FPAR * PAR

    vapor_stress = 0.88 - 0.183 * np.log(esat_24 - eact_24)
    vapor_stress_biomass = vapor_stress.clip(0.0, 1.0)
    Jarvis_coeff = (Th - Kt) / (Kt - Tl)
    heat_stress_biomass = ((Temp_24 - Tl) * np.power(Th - Temp_24, Jarvis_coeff) /
                           ((Kt - Tl) * np.power(Th - Kt, Jarvis_coeff)))
    print 'vapor stress biomass =', '%0.3f' % np.nanmean(vapor_stress_biomass)
    print 'heat stress biomass =', '%0.3f' % np.nanmean(heat_stress_biomass)

    # Light use efficiency, reduced below its potential value by low
    # temperature or water shortage:
    LUE = (LUEmax * heat_stress_biomass * vapor_stress_biomass * moisture_stress_biomass)

    # Dry matter production (kg/ha/d):
    Biomass_prod = APAR * LUE * 0.864             # C3 vegetation

    # Water productivity
    Biomass_wp = Biomass_prod/ (ETA_24 * 10)  # C3 vegetation
    Biomass_wp[ETA_24 == 0.0] = 0.0

    # Water deficit
    Biomass_deficit = (Biomass_prod / moisture_stress_biomass -
                          Biomass_prod)


    return(LUE,Biomass_prod,Biomass_wp,Biomass_deficit)

#------------------------------------------------------------------------------
def Classify_Irrigation(moisture_stress_biomass, vegt_cover):
   '''
   This function makes a classification with 4 categories which show the irrigation needs
   '''
   for_irrigation = np.copy(moisture_stress_biomass)

   # make a discreed irrigation needs map with the following categories
      # Irrigation needs:
       # 0: No need for irrigation
       # 1: Perhaps irrigate
       # 2: Irrigate
       # 3: Irrigate immediately
   irrigation_needs = np.copy(for_irrigation)
   irrigation_needs[np.where(irrigation_needs >= 1.0)] == 0.0
   irrigation_needs[np.logical_and(irrigation_needs >= 0.9, irrigation_needs < 1.0)] = 1.0
   irrigation_needs[np.where((irrigation_needs >= 0.8) & (irrigation_needs < 0.9))] = 2.0
   irrigation_needs[np.where(irrigation_needs < 0.8)] = 3.0
   irrigation_needs[vegt_cover <= 0.3] = 0.0
   return(irrigation_needs)

#------------------------------------------------------------------------------
def Separate_E_T(LAI,ETP_24,Theta_res_top,Theta_res_sub, Theta_sat_top, Theta_sat_sub, top_soil_moisture,sl_es_24, Psychro_c,moisture_stress_biomass_first,vegt_cover,ETA_24,SM_stress_trigger,root_zone_moisture_first,total_soil_moisture):
   '''
   Separate the Evapotranspiration into evaporation and Transpiration
   '''
   # constants
   Light_use_extinction_factor = 0.5    # Light use extinction factor for Bear's Law

   Tpot_24_estimate=(1-np.exp(-Light_use_extinction_factor*LAI))*ETP_24
   SE_top = (top_soil_moisture-Theta_res_top)/(Theta_sat_top-Theta_res_top)
   Eact_24_estimate=np.minimum(1,1 / np.power(SE_top + 0.1,-2.0))*(ETP_24-Tpot_24_estimate)
   #RS_soil = RS_soil_min * np.power(SE_top,-2.0)
   #Eact_24_estimate=(sl_es_24+Psychro_c*(1+RS_soil_min/Rah_PM))/(sl_es_24+Psychro_c*(1+RS_soil/Rah_PM))*(ETP_24-Tpot_24_estimate)
   n66_memory = moisture_stress_biomass_first * Tpot_24_estimate

   # calulate the first estimation of actual daily tranpiration
   Tact_24_estimate = np.copy(n66_memory)
   Tact_24_estimate[n66_memory > 0.99*ETA_24]=ETA_24[n66_memory > 0.99*ETA_24]
   Tact_24_estimate[vegt_cover == 0.0] = 0.0

   # calculate the second estimation and end estimation of the actual daily tranpiration
   Tact_24 = np.abs((Tact_24_estimate/(Tact_24_estimate + Eact_24_estimate))*ETA_24)

   # calculate the actual daily potential transpiration
   Tpot_24 = np.copy(Tpot_24_estimate)
   Tpot_24[Tpot_24_estimate < Tact_24] = Tact_24[Tpot_24_estimate < Tact_24]

   # calculate moisture stress biomass
   moisture_stress_biomass = Tact_24 / Tpot_24

   # Calculate root zone moisture final
   Se_Poly=2.23*np.power(moisture_stress_biomass,3)-3.35*np.power(moisture_stress_biomass,2)+1.98*moisture_stress_biomass+0.07
   root_zone_moisture1=Se_Poly*(SM_stress_trigger+0.02-Theta_res_sub)+Theta_res_sub
   root_zone_moisture_final=np.where(root_zone_moisture1>root_zone_moisture_first,root_zone_moisture1,root_zone_moisture_first)

   # Calculate top zone moisture final
   top_zone_moisture1=(total_soil_moisture-root_zone_moisture_final*vegt_cover)/(1-vegt_cover)
   top_zone_moisture_final=top_zone_moisture1.clip(Theta_res_top,Theta_sat_top)

   # calculate the actual daily evaporation
   Eact_24 = ETA_24 - Tact_24

   # calculate the Transpiration deficit
   T24_deficit = Tpot_24 - Tact_24

   # calculate the beneficial fraction
   beneficial_fraction=Tact_24 / ETA_24
   beneficial_fraction[ETA_24 == 0.0] = 0.0

   return(Eact_24,Tpot_24,Tact_24,moisture_stress_biomass,T24_deficit,beneficial_fraction,root_zone_moisture_final,top_zone_moisture_final)

#------------------------------------------------------------------------------
def Calc_Soil_Moisture(ETA_24,EF_inst,QC_Map, water_mask,vegt_cover,Theta_sat_top, Theta_sat_sub,Theta_res_top, Theta_res_sub,depl_factor,Field_Capacity,FPAR, Soil_moisture_wilting_point):
    """
    Function to calculate soil characteristics
    """
    # constants:
    Veg_Cover_Threshold_RZ = 0.9 # Threshold vegetation cover for root zone moisture

    # Average fraction of TAW that can be depleted from the root zone
    # before stress:
    p_factor = depl_factor + 0.04 * (5.0 - ETA_24)  # page 163 of FAO 56
    # The factor p differs from one crop to another. It normally varies from
    # 0.30 for shallow rooted plants at high rates of ETc (> 8 mm d-1)
    # to 0.70 for deep rooted plants at low rates of ETc (< 3 mm d-1)

    # Critical value under which plants get stressed:
    SM_stress_trigger = Field_Capacity - p_factor * (Field_Capacity - Soil_moisture_wilting_point)
    EF_inst[EF_inst >= 1.0] = 0.999

    # Total soil water content (cm3/cm3):
    total_soil_moisture = Theta_sat_sub * np.exp((EF_inst - 1.0) / 0.421)   # asce paper Scott et al. 2003
    total_soil_moisture[np.logical_or(water_mask == 1.0,QC_Map == 1.0)] = 1.0  # In water and snow is 1
    total_soil_moisture[QC_Map == 1.0] = np.nan  # Where clouds no data

    # Root zone soil moisture:
    RZ_SM = np.copy(total_soil_moisture)
    RZ_SM[vegt_cover <= Veg_Cover_Threshold_RZ] = np.nan
    if np.isnan(np.nanmean(RZ_SM)) == True:
        Veg_Cover_Threshold_RZ = np.nanpercentile(vegt_cover, 80)
        RZ_SM = np.copy(total_soil_moisture)
        RZ_SM[vegt_cover <= Veg_Cover_Threshold_RZ] = np.nan
        print 'No RZ_SM so the vegetation Threshold for RZ is adjusted from 0,9 to =', '%0.3f' % Veg_Cover_Threshold_RZ

    #RZ_SM = RZ_SM.clip(Theta_res, (0.85 * Theta_sat))
    #RZ_SM[np.logical_or(water_mask == 1.0, water_mask == 2.0)] = 1.0
    RZ_SM_NAN = np.copy(RZ_SM)
    RZ_SM_NAN[RZ_SM==0] = np.nan
    RZ_SM_min = np.nanmin(RZ_SM_NAN)
    RZ_SM_max = np.nanmax(RZ_SM_NAN)
    RZ_SM_mean = np.nanmean(RZ_SM_NAN)
    print 'Root Zone Soil moisture mean =', '%0.3f (cm3/cm3)' % RZ_SM_mean
    print 'Root Zone Soil moisture min =', '%0.3f (cm3/cm3)' % RZ_SM_min
    print 'Root Zone Soil moisture max =', '%0.3f (cm3/cm3)' % RZ_SM_max

    Max_moisture_RZ = vegt_cover * (RZ_SM_max - RZ_SM_mean) + RZ_SM_mean

    # Soil moisture in the top (temporary)
    top_soil_moisture_temp = np.copy(total_soil_moisture)
    top_soil_moisture_temp[np.logical_or(vegt_cover <= 0.02, vegt_cover >= 0.1)] = 0
    top_soil_moisture_temp[top_soil_moisture_temp == 0] = np.nan
    top_soil_moisture_std = np.nanstd(top_soil_moisture_temp)
    top_soil_moisture_mean = np.nanmean(top_soil_moisture_temp)
    print 'Top Soil moisture mean =', '%0.3f (cm3/cm3)' % top_soil_moisture_mean
    print 'Top Soil moisture Standard Deviation', '%0.3f (cm3/cm3)' % top_soil_moisture_std

    # calculate root zone moisture
    root_zone_moisture_temp = (total_soil_moisture - (top_soil_moisture_mean + top_soil_moisture_std) * (1-vegt_cover))/vegt_cover # total soil moisture = soil moisture no vegtatation *(1-vegt_cover)+soil moisture root zone * vegt_cover
    try:
        root_zone_moisture_temp[root_zone_moisture_temp <= Theta_res_sub] = Theta_res_sub[root_zone_moisture_temp <= Theta_res_sub]
    except:
        root_zone_moisture_temp[root_zone_moisture_temp <= Theta_res_sub] = Theta_res_sub


    root_zone_moisture_temp[root_zone_moisture_temp >= Max_moisture_RZ] = Max_moisture_RZ[root_zone_moisture_temp >= Max_moisture_RZ]

    root_zone_moisture_first = np.copy(root_zone_moisture_temp)
    root_zone_moisture_first[np.logical_or(QC_Map ==1.0 ,np.logical_or(water_mask == 1.0, vegt_cover < 0.0))] = 0

    # Normalized stress trigger:
    norm_trigger = (root_zone_moisture_first - Soil_moisture_wilting_point)/ (SM_stress_trigger + 0.02 - Soil_moisture_wilting_point)
    norm_trigger[norm_trigger > 1.0] = 1.0

    # moisture stress biomass:
    moisture_stress_biomass_first = norm_trigger - (np.sin(2 * np.pi * norm_trigger)) / (2 * np.pi)
    moisture_stress_biomass_first=np.where(moisture_stress_biomass_first<0.5*FPAR,0.5*FPAR,moisture_stress_biomass_first)
    moisture_stress_biomass_first[moisture_stress_biomass_first <= 0.0] = 0
    moisture_stress_biomass_first[moisture_stress_biomass_first > 1.0] = 1.0

     # Soil moisture in the top layer - Recalculated ??
    top_soil_moisture = ((total_soil_moisture - root_zone_moisture_first * vegt_cover) / (1.0 - vegt_cover))

    try:
        top_soil_moisture[top_soil_moisture > Theta_sat_top] = Theta_sat_top [top_soil_moisture > Theta_sat_top]
    except:
        top_soil_moisture[top_soil_moisture > Theta_sat_top] = Theta_sat_top

    top_soil_moisture[np.logical_or(water_mask == 1.0, QC_Map == 1.0)] = 1.0

    return(SM_stress_trigger, total_soil_moisture, root_zone_moisture_first, moisture_stress_biomass_first,top_soil_moisture,RZ_SM_NAN)

#------------------------------------------------------------------------------
def Calc_Bulk_surface_resistance(sl_es_24,Rn_24,Refl_rad_water,air_dens,esat_24,eact_24,rah_pm_act,ETA_24,Lhv,Psychro_c):
    """
    Function to calculate the bulk surface resistance
    """
    # Bulk surface resistance (s/m):
    bulk_surf_resis_24 = ((((sl_es_24 * (Rn_24 - Refl_rad_water) + air_dens *
                          1004 * (esat_24 - eact_24) / rah_pm_act) / (ETA_24 * Lhv / 86400) -
                          sl_es_24) / Psychro_c - 1.0) * rah_pm_act)
    bulk_surf_resis_24[ETA_24 <= 0.0] = 100000.0
    bulk_surf_resis_24 = bulk_surf_resis_24.clip(0.0, 100000.0)
    return(bulk_surf_resis_24)

#------------------------------------------------------------------------------
def Calc_ETact(esat_24, eact_24, EF_inst, Rn_24, Refl_rad_water, Lhv, Image_Type):
    """
    Function to calculate the daily evaporation
    """
    # Advection factor
    if Image_Type == 2:
        AF = np.ones(Rn_24.shape)
    else:
        AF = 1 + 0.985 * (np.exp((esat_24 - eact_24) * 0.08) - 1.0) * EF_inst
    # Daily evapotranspiration:
    ETA_24 = EF_inst * AF * (Rn_24 - Refl_rad_water) / (Lhv * 1000) * 86400000
    ETA_24=ETA_24.clip(0,15.0)
    return(ETA_24, AF)

#------------------------------------------------------------------------------
def Calc_instantaneous_ET_fraction(LE_inst,rn_inst,g_inst):
    """
    Function to calculate the evaporative fraction
    """
    EF_inst = LE_inst / (rn_inst - g_inst)    # Evaporative fraction
    EF_inst = EF_inst.clip(0.0, 1.8)
    EF_inst[LE_inst<0] = 0

    return(EF_inst)

#------------------------------------------------------------------------------
def Calc_Ref_Pot_ET(LAI,Surface_temp,sl_es_24,Rn_ref,air_dens,esat_24,eact_24,rah_grass,Psychro_c,Rn_24,Refl_rad_water,rah_pm_pot):
    """
    Function to calculate the reference potential evapotransporation and potential evaporation
    """

    # Constants
    rl = 130                         # Bulk stomatal resistance of the well-illuminated leaf (s/m)

    # Effective leaf area index involved, see Allen et al. (2006):
    LAI_eff = LAI / (0.3 * LAI + 1.2)
    rs_min = rl / LAI_eff  # Min (Bulk) surface resistance (s/m)
    # Latent heat of vaporization (J/kg):
    Lhv = (2.501 - 2.361e-3 * (Surface_temp - 273.15)) * 1E6

    # Reference evapotranspiration- grass
    # Penman-Monteith of the combination equation (eq 3 FAO 56) (J/s/m2)
    LET_ref_24 = ((sl_es_24 * Rn_ref + air_dens * 1004 * (esat_24 - eact_24) /
                  rah_grass) / (sl_es_24 + Psychro_c * (1 + 70.0/rah_grass)))
    # Reference evaportranspiration (mm/d):
    ETref_24 = LET_ref_24 / (Lhv * 1000) * 86400000

    # Potential evapotranspiration
    # Penman-Monteith of the combination equation (eq 3 FAO 56) (J/s/m2)
    LETpot_24 = ((sl_es_24 * (Rn_24 - Refl_rad_water) + air_dens * 1004 *
               (esat_24 - eact_24)/rah_pm_pot) / (sl_es_24 + Psychro_c * (1 + rs_min/rah_pm_pot)))
    # Potential evaportranspiration (mm/d)
    ETpot_24 = LETpot_24 / (Lhv * 1000) * 86400000
    ETpot_24[ETpot_24 > 15.0] = 15.0
    return(ETpot_24,ETref_24,Lhv,rs_min)

#------------------------------------------------------------------------------
def Calc_Rn_Ref(shape_lsc,water_mask,Rn_24,Ra_mountain_24,Transm_24,Rnl_24_FAO,Wind_24):
    """
    Function to calculate the net solar radiation
    """
    # constants:
    G24_water = 0.1  # G24 ratio for water - reflectivity?

    # Reflected radiation at water surface: ??
    Refl_rad_water = np.zeros((shape_lsc[1], shape_lsc[0]))
    Refl_rad_water = np.where(water_mask != 0.0, G24_water * Rn_24, 0.0)

    # Aerodynamic resistance (s/m) for grass surface:
    rah_grass = 208.0 / Wind_24
    print  'rah_grass=', '%0.3f (s/m)' % np.nanmean(rah_grass)
    # Net radiation for grass Rn_ref, eq 40, FAO56:
    Rn_ref = Ra_mountain_24 * Transm_24 * (1 - 0.23) - Rnl_24_FAO  # Rnl avg(fao-slob)?
    return(Rn_ref, Refl_rad_water,rah_grass)


#------------------------------------------------------------------------------
def Iterate_Friction_Velocity(k_vk,u_200,Surf_roughness,g_inst,rn_inst, ts_dem, ts_dem_hot, ts_dem_cold,air_dens, Surface_temp,L,psi,psi_m200,psi_m200_stable,QC_Map, hot_pixels, slope):
    """
    Function to correct the windspeed and aerodynamic resistance for the iterative process the output can be used as the new input for this model
    """
    # Sensible heat 2 (Step 6)
    # Corrected value for the friction velocity, unstable
    ustar_corr_unstable = (k_vk * u_200 / (np.log(200.0 / Surf_roughness) -
                        psi_m200))
    # Corrected value for the friction velocity, stable
    ustar_corr_stable = (k_vk * u_200 / (np.log(200.0 / Surf_roughness) -
                      psi_m200_stable))
    ustar_corr = np.where(L > 0.0, ustar_corr_stable, ustar_corr_unstable)
    ustar_corr[ustar_corr < 0.02] = 0.02

    rah_corr_unstable = (np.log(2.0/0.01) - psi) / (k_vk * ustar_corr)     # unstable
    rah_corr_stable = (np.log(2.0/0.01) - 0.0) / (k_vk * ustar_corr)       # stable
    rah_corr = np.where(L > 0.0, rah_corr_stable, rah_corr_unstable)
    L_corr, psi_m200_corr_stable, psi_corr, psi_m200_corr,h,dT, slope_dt, offset_dt = sensible_heat(
            rah_corr, ustar_corr, rn_inst, g_inst, ts_dem, ts_dem_hot, ts_dem_cold,
            air_dens, Surface_temp, k_vk,QC_Map, hot_pixels, slope)
    return(L_corr,psi_corr,psi_m200_corr,psi_m200_corr_stable,h,ustar_corr,rah_corr,dT,slope_dt, offset_dt)

#------------------------------------------------------------------------------
def Calc_Wind_Speed_Friction(h_obst,Wind_inst,zx,LAI,NDVI,Surf_albedo,water_mask,surf_roughness_equation_used):
    """
    Function to calculate the windspeed and friction by using the Raupach or NDVI model
    """

    # constants
    k_vk = 0.41      # Von Karman constant
    h_grass = 0.12   # Grass height (m)
    cd = 53          # Free parameter for displacement height, default = 20.6
    # 1) Raupach model
    zom_Raupach=Raupach_Model(h_obst,cd,LAI)

    # 2) NDVI model
    zom_NDVI=NDVI_Model(NDVI,Surf_albedo,water_mask)

    if surf_roughness_equation_used == 1:
        Surf_roughness = zom_NDVI
    else:
        Surf_roughness = zom_Raupach

    zom_grass = 0.123 * h_grass
    # Friction velocity for grass (m/s):
    ustar_grass = k_vk * Wind_inst / np.log(zx / zom_grass)
    print 'u*_grass = ', '%0.3f (m/s)' % np.mean(ustar_grass)
    # Wind speed (m/s) at the "blending height" (200m):
    u_200 = ustar_grass * np.log(200 / zom_grass) / k_vk
    print 'Wind speed at the blending height, u200 =', '%0.3f (m/s)' % np.mean(u_200)
    # Friction velocity (m/s):
    ustar_1 = k_vk * u_200 / np.log(200 / Surf_roughness)
    return(Surf_roughness,u_200,ustar_1)

#------------------------------------------------------------------------------
def Raupach_Model(h_obst,cd,LAI):
    """
    Function for the Raupach model to calculate the surface roughness (based on Raupach 1994)
    """
    # constants
    cw = 2.0
    LAIshelter = 2.5

    # calculate psi
    psi = np.log(cw) - 1 + np.power(2.0, -1)  # Vegetation influence function

    # Calculate Ustar divided by U
    ustar_u = np.power((0.003+0.3*LAI/2), 0.5)
    ustar_u[LAI<LAIshelter] = 0.3

    # calculate: 1 - d/hv
    inv_d_hv =np.power((1-np.exp(np.power(-1*(cd*LAI),0.5)))/(cd * LAI),0.5)

    # Calculate: surface roughness/hv
    zom_hv = inv_d_hv * np.exp(-0.41/ustar_u-psi)

    # Calculate: surface roughness
    zom_Raupach = zom_hv * h_obst

    return(zom_Raupach)

#------------------------------------------------------------------------------
def NDVI_Model(NDVI,Surf_albedo,water_mask):
    """
    Function for the NDVI model to calculate the surface roughness
    """
    zom_NDVI = np.exp(1.096 * NDVI / Surf_albedo - 5.307)
    zom_NDVI[water_mask == 1.0] = 0.001
    zom_NDVI[zom_NDVI > 10.0] = 10.0
    return(zom_NDVI)

#------------------------------------------------------------------------------
def Correct_Surface_Temp(Surface_temp,Temp_lapse_rate,DEM_resh,Pair,dr,Transm_corr,cos_zn,Sun_elevation,deg2rad,ClipLandsat):
    """
    Function to correct the surface temperature based on the DEM map
    """
    #constants:
    Gsc = 1367        # Solar constant (W / m2)

    cos_zenith_flat = np.cos((90 - Sun_elevation) * deg2rad)
    Temp_corr = Surface_temp + Temp_lapse_rate * DEM_resh  # rescale everything to sea level
    Temp_corr[Surface_temp == 350.0] = 0.0
    air_dens = 1000 * Pair / (1.01 * Surface_temp * 287)
    #
    ts_dem = (Temp_corr + (Gsc * dr * Transm_corr * cos_zn -
              Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))
    #(Temp_corr - (Gsc * dr * Transm_corr * cos_zn -
    #          Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))
    ts_dem[ClipLandsat==1]=np.nan
    ts_dem[ts_dem==0]=np.nan
    ts_dem[ts_dem<273]=np.nan
    ts_dem[ts_dem>350]=np.nan

    return(ts_dem,air_dens,Temp_corr)

#------------------------------------------------------------------------------
def Calc_Hot_Pixels(ts_dem,QC_Map, water_mask, NDVI,Hot_Pixel_Constant, ts_dem_cold):
    """
    Function to calculates the hot pixels based on the surface temperature and NDVI
    """

    # Limits for hot pixels selection
    NDVIhot_low = 0.03               # Lower NDVI treshold for hot pixels
    NDVIhot_high = 0.20              # Higher NDVI treshold for hot pixels

    for_hot = np.copy(ts_dem)
    for_hot[NDVI <= NDVIhot_low] = 0.0
    for_hot[NDVI >= NDVIhot_high] = 0.0
    for_hot[np.logical_or(water_mask != 0.0, QC_Map != 0.0)] = 0.0
    hot_pixels = np.copy(for_hot)
    hot_pixels[for_hot < ts_dem_cold] = np.nan
    ts_dem_hot_max = np.nanmax(hot_pixels)    # Max
    ts_dem_hot_mean = np.nanmean(hot_pixels)  # Mean
    ts_dem_hot_std = np.nanstd(hot_pixels)    # Standard deviation
    #ts_dem_hot = ts_dem_hot_max - 0.25 * ts_dem_hot_std
    #ts_dem_hot = (ts_dem_hot_max + ts_dem_hot_mean)/2


    ts_dem_hot=ts_dem_hot_mean + Hot_Pixel_Constant * ts_dem_hot_std


    print 'hot : max= %0.3f (Kelvin)' % ts_dem_hot_max, ', sd= %0.3f (Kelvin)' % ts_dem_hot_std, \
           ', mean= %0.3f (Kelvin)' % ts_dem_hot_mean, ', value= %0.3f (Kelvin)' % ts_dem_hot
    return(ts_dem_hot,hot_pixels)

#------------------------------------------------------------------------------
def Calc_Cold_Pixels(ts_dem,water_mask,QC_Map,ts_dem_cold_veg,Cold_Pixel_Constant):
    """
    Function to calculates the the cold pixels based on the surface temperature
    """
    for_cold = np.copy(ts_dem)
    for_cold[water_mask != 1.0] = 0.0
    for_cold[QC_Map != 0] = 0.0
    cold_pixels = np.copy(for_cold)
    cold_pixels[for_cold < 278.0] = np.nan
    cold_pixels[for_cold > 320.0] = np.nan
    # cold_pixels[for_cold < 285.0] = 285.0
    ts_dem_cold_std = np.nanstd(cold_pixels)     # Standard deviation
    ts_dem_cold_min = np.nanmin(cold_pixels)     # Min
    ts_dem_cold_mean = np.nanmean(cold_pixels)   # Mean

    # If average temperature is below zero or nan than use the vegetation cold pixel
    if ts_dem_cold_mean <= 0.0:
        ts_dem_cold = ts_dem_cold_veg + Cold_Pixel_Constant * ts_dem_cold_std
    if np.isnan(ts_dem_cold_mean) == True:
        ts_dem_cold = ts_dem_cold_veg + Cold_Pixel_Constant * ts_dem_cold_std
    else:
        ts_dem_cold = ts_dem_cold_mean + Cold_Pixel_Constant * ts_dem_cold_std

    if ts_dem_cold > ts_dem_cold_veg:
        ts_dem_cold = ts_dem_cold_veg

    print 'cold water: min=%0.3f (Kelvin)' %ts_dem_cold_min , ', sd= %0.3f (Kelvin)' % ts_dem_cold_std, \
           ', mean= %0.3f (Kelvin)' % ts_dem_cold_mean, ', value= %0.3f (Kelvin)' % ts_dem_cold
    return(ts_dem_cold,cold_pixels,ts_dem_cold_mean)

#------------------------------------------------------------------------------
def Calc_Cold_Pixels_Veg(NDVI,QC_Map,ts_dem, Cold_Pixel_Constant):
    """
    Function to calculates the the cold pixels based on vegetation
    """
    cold_pixels_vegetation = np.copy(ts_dem)
    NDVI_std = np.nanstd(NDVI)
    NDVI_max = np.nanmax(NDVI)
    cold_pixels_vegetation[np.logical_and(NDVI <= (NDVI_max-0.1*NDVI_std),QC_Map != 0.0)] = 0.0
    cold_pixels_vegetation[cold_pixels_vegetation==0.0] = np.nan
    ts_dem_cold_std_veg = np.nanstd(cold_pixels_vegetation)
    ts_dem_cold_mean_veg = np.nanmean(cold_pixels_vegetation)
    ts_dem_cold_veg = ts_dem_cold_mean_veg + Cold_Pixel_Constant * ts_dem_cold_std_veg

    return(ts_dem_cold_veg)

#------------------------------------------------------------------------------
def Calc_Meteo(Rs_24,eact_24,Temp_24,Surf_albedo,dr,tir_emis,Surface_temp,water_mask,NDVI,Transm_24,SB_const,lw_in_inst,Rs_in_inst):
    """
    Calculates the instantaneous Ground heat flux and solar radiation.
    """

    # Net shortwave radiation (W/m2):
    Rns_24 = Rs_24 * (1 - Surf_albedo)


    # Net outgoing longwave radiation (W/m2):
    Rnl_24_FAO = (SB_const * np.power(Temp_24 + 273.15, 4) * (0.34-0.14 *
                  np.power(eact_24, 0.5)) * (1.35 * Transm_24 / 0.8 - 0.35))

    Rnl_24_Slob = 110 * Transm_24
    print 'Mean Daily Net Radiation (Slob) = %0.3f (W/m2)' % np.nanmean(Rnl_24_Slob)

    # Net 24 hrs radiation (W/m2):
    Rn_24_FAO = Rns_24 - Rnl_24_FAO          # FAO equation
    Rn_24_Slob = Rns_24 - Rnl_24_Slob       # Slob equation
    Rn_24 = (Rn_24_FAO + Rn_24_Slob) / 2  # Average


    # Instantaneous outgoing longwave radiation:
    lw_out_inst = tir_emis * SB_const * np.power(Surface_temp, 4)

    # Instantaneous net radiation
    rn_inst = (Rs_in_inst * (1 - Surf_albedo) + lw_in_inst - lw_out_inst -
               (1 - tir_emis) * lw_in_inst)
    # Instantaneous Soil heat flux
    g_inst = np.where(water_mask != 0.0, 0.4 * rn_inst,
                      ((Surface_temp - 273.15) * (0.0038 + 0.0074 * Surf_albedo) *
                       (1 - 0.978 * np.power(NDVI, 4))) * rn_inst)
    return(Rn_24,rn_inst,g_inst,Rnl_24_FAO)

#------------------------------------------------------------------------------
def Calc_surface_water_temp(Temp_inst,Landsat_nr,Lmax,Lmin,therm_data,b10_emissivity,k1_c,k2_c,eact,shape_lsc,water_mask_temp,Bands_thermal,Rp,tau_sky,surf_temp_offset,Image_Type):
    """
    Calculates the surface temperature and create a water mask
    """

    # Spectral radiance for termal
    if Landsat_nr == 8:
        if Bands_thermal == 1:
            k1 = k1_c[0]
            k2 = k2_c[0]
            L_lambda_b10 = (Lmax[-1] - Lmin[-1]) / (65535-1) * therm_data[:, :, 0] + Lmin[-1]

            # Get Temperature
            Temp_TOA = Get_Thermal(L_lambda_b10,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

        elif Bands_thermal == 2:
            L_lambda_b10 = (Lmax[-2] - Lmin[-2]) / (65535-1) * therm_data[:, :, 0] + Lmin[-2]
            L_lambda_b11 = (Lmax[-1] - Lmin[-1]) / (65535-1) * therm_data[:, :, 1] + Lmin[-1]

            # Brightness temperature
            # From Band 10:
            Temp_TOA_10 = (k2_c[0] / np.log(k1_c[0] / L_lambda_b10 + 1.0))
            # From Band 11:
            Temp_TOA_11 = (k2_c[1] / np.log(k1_c[1] / L_lambda_b11 + 1.0))
            # Combined:
            Temp_TOA = (Temp_TOA_10 + 1.378 * (Temp_TOA_10 - Temp_TOA_11) +
                           0.183 * np.power(Temp_TOA_10 - Temp_TOA_11, 2) - 0.268 +
                           (54.30 - 2.238 * eact) * (1 - b10_emissivity))

    elif Landsat_nr == 7:
        k1=666.09
        k2=1282.71
        L_lambda_b6 = (Lmax[-1] - Lmin[-1]) / (256-1) * therm_data[:, :, 0] + Lmin[-1]

        # Brightness temperature - From Band 6:
        Temp_TOA = Get_Thermal(L_lambda_b6,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

    elif Landsat_nr == 5:
        k1=607.76
        k2=1260.56
        L_lambda_b6 = ((Lmax[-1] - Lmin[-1]) / (256-1) * therm_data[:, :, 0] +
                       Lmin[-1])

       # Brightness temperature - From Band 6:
        Temp_TOA = Get_Thermal(L_lambda_b6,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

    # Surface temperature
    Surface_temp = Temp_TOA
    Surface_temp = Surface_temp.clip(230.0, 360.0)

    # Cloud mask:
    temp_water = np.zeros((shape_lsc[1], shape_lsc[0]))
    temp_water = np.copy(Surface_temp)
    temp_water[water_mask_temp == 0.0] = np.nan
    temp_water_sd = np.nanstd(temp_water)     # Standard deviation
    temp_water_mean = np.nanmean(temp_water)  # Mean
    print 'Mean water temperature = ', '%0.3f (Kelvin)' % temp_water_mean
    print 'SD water temperature = ', '%0.3f (Kelvin)' % temp_water_sd
    cloud_mask = np.zeros((shape_lsc[1], shape_lsc[0]))
    cloud_mask[Surface_temp < np.minimum((temp_water_mean - 1.0 * temp_water_sd -
               surf_temp_offset),290)] = 1.0

    return(Surface_temp,cloud_mask)


#------------------------------------------------------------------------------
def Get_Thermal(lambda_b10,Rp,Temp_inst,tau_sky,TIR_Emissivity,k1,k2):

    # Narrow band downward thermal radiation from clear sky, rsky (W/m2/sr/µm)
    rsky = (1.807E-10 * np.power(Temp_inst + 273.15, 4) * (1 - 0.26 *
            np.exp(-7.77E-4 * np.power((-Temp_inst), -2))))
    print 'Rsky = ', '%0.3f (W/m2/sr/µm)' % np.nanmean(rsky)

    # Corrected thermal radiance from the surface, Wukelikc et al. (1989):
    correc_lambda_b10 = ((lambda_b10 - Rp) / tau_sky -
                               (1.0 - TIR_Emissivity) * rsky)
    # Brightness temperature - From Band 10:
    Temp_TOA = (k2 / np.log(TIR_Emissivity * k1 /
                       correc_lambda_b10 + 1.0))

    return(Temp_TOA)
#------------------------------------------------------------------------------
def Calc_vegt_para(NDVI,water_mask_temp,shape_lsc):
    """
    Calculates the Fraction of PAR, Thermal infrared emissivity, Nitrogen, Vegetation Cover, LAI, b10_emissivity
    """
    # Fraction of PAR absorbed by the vegetation canopy (FPAR):
    FPAR = -0.161 + 1.257 * NDVI
    FPAR[NDVI < 0.125] = 0.0

    # Termal infrared emissivity
    tir_emis = 1.009 + 0.047 * np.log(NDVI)
    tir_emis[np.logical_or(water_mask_temp == 1.0, water_mask_temp == 2.0)] = 1.0
    tir_emis[np.logical_and(NDVI < 0.125, water_mask_temp == 0.0)] = 0.92

    # Vegetation Index - Regression model from Bagheri et al. (2013)
    VI = 38.764 * np.square(NDVI) - 24.605 * NDVI + 5.8103

    # Nitrogen computation
    Nitrogen = np.copy(VI)
    Nitrogen[VI <= 0.0] = 0.0
    Nitrogen[NDVI <= 0.0] = 0.0

    # Vegetation cover:
    vegt_cover = 1 - np.power((0.8 - NDVI)/(0.8 - 0.125), 0.7)
    vegt_cover[NDVI < 0.125] = 0.0
    vegt_cover[NDVI > 0.8] = 0.99

    # Leaf Area Index (LAI)
    LAI_1 = np.log(-(vegt_cover - 1)) / -0.45
    LAI_1[LAI_1 > 8] = 8.0
    LAI_2 = (9.519 * np.power(NDVI, 3) + 0.104 * np.power(NDVI, 2) +
             1.236 * NDVI - 0.257)

    LAI = (LAI_1 + LAI_2) / 2.0  # Average LAI
    LAI[LAI < 0.001] = 0.001

    b10_emissivity = np.zeros((shape_lsc[1], shape_lsc[0]))
    b10_emissivity = np.where(LAI <= 3.0, 0.95 + 0.01 * LAI, 0.98)
    b10_emissivity[water_mask_temp != 0.0] = 1.0

    return(FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity)

#------------------------------------------------------------------------------
def Calc_Ra_Mountain(lon,DOY,hour,minutes,lon_proy,lat_proy,slope,aspect):
    """
    Calculates the extraterrestiral solar radiation by using the date, slope and aspect.
    """

    # Constants
    deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
    Min_cos_zn = 0.1  # Min value for cos zenith angle
    Max_cos_zn = 1.0  # Max value for cos zenith angle
    Gsc = 1367        # Solar constant (W / m2)
    try:
        Loc_time = float(hour) + float(minutes)/60  # Local time (hours)
    except:
        Loc_time = np.float_(hour) + np.float_(minutes)/60  # Local time (hours)
    # Rounded difference of the local time from Greenwich (GMT) (hours):
    offset_GTM = round(np.sign(lon[int(lon.shape[0])/2, int(lon.shape[1])/2]) * lon[int(lon.shape[0])/2,int(lon.shape[1])/2] * 24 / 360)

    print '  Local time: ', '%0.3f' % np.nanmean(Loc_time)
    print '  Difference of local time (LT) from Greenwich (GMT): ', offset_GTM

    # 1. Calculation of extraterrestrial solar radiation for slope and aspect
    # Computation of Hour Angle (HRA = w)
    B = 360./365 * (DOY-81)           # (degrees)
    # Computation of cos(theta), where theta is the solar incidence angle
    # relative to the normal to the land surface
    delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)
    phi = lat_proy * deg2rad                                     # latitude of the pixel (radians)
    s = slope * deg2rad                                          # Surface slope (radians)
    gamma = (aspect-180) * deg2rad                               # Surface aspect angle (radians)
    w=w_time(Loc_time, lon_proy, DOY)                            # Hour angle (radians)
    a,b,c = Constants(delta,s,gamma,phi)
    cos_zn= AngleSlope(a,b,c,w)
    cos_zn = cos_zn.clip(Min_cos_zn, Max_cos_zn)

    print 'Average Cos Zenith Angle: ', '%0.3f (Radians)' % np.nanmean(cos_zn)

    dr = 1 + 0.033 * cos(DOY*2*pi/365)  # Inverse relative distance Earth-Sun
    # Instant. extraterrestrial solar radiation (W/m2), Allen et al.(2006):
    Ra_inst = Gsc * cos_zn * dr

    # 24-hours extraterrestrial radiation
    # 1.) determine if there are one or two periods of sun
    # 2.) calculate the 24-hours extraterrestrial radiation if there are two periods of sun
    # 3.) calculate the 24-hours extraterrestrial radiation if there is one period of sun

    #1.) determine amount of sun periods
    Ra_24 = np.zeros(np.shape(lat_proy))*np.nan
    constant=Gsc*dr/(2*np.pi)
    TwoPeriod= TwoPeriods(delta,s,phi)  # all input in radians

    #2.) calculate the 24-hours extraterrestrial radiation (2 periods)
    ID = np.where(np.ravel(TwoPeriod==True))
    Ra_24.flat[ID]=TwoPeriodSun(constant,delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])

    #3.) calculate the 24-hours extraterrestrial radiation (1 period)
    ID = np.where(np.ravel(TwoPeriod==False))
    Ra_24.flat[ID]=OnePeriodSun(constant,delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])

    # Horizontal surface
    ws = np.arccos(-np.tan(delta) * np.tan(phi))  # Sunrise/sunset time angle

    # Extraterrestial radiation for a horizontal surface for 24-h period:
    Ra_hor_24 = (Gsc * dr / np.pi * (np.sin(delta) * np.sin(phi) * ws + np.cos(delta) * np.cos(phi) * np.sin(ws)))
    # cos_theta_flat = (np.sin(delta) * np.sin(phi) + np.cos(delta) * np.cos(phi) * np.cos(w))

    # Mountain radiation
    Ra_mountain_24 = np.where(Ra_24 > Min_cos_zn * Ra_hor_24, Ra_24 / np.cos(s),
                           Ra_hor_24)
    Ra_mountain_24[Ra_mountain_24 > 600.0] = 600.0

    return(Ra_mountain_24,Ra_inst,cos_zn,dr,phi,delta)

#------------------------------------------------------------------------------
def OnePeriodSun(constant,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006
    Calculate the 24-hours extraterrestrial radiation when there is one sun period
    '''
    sunrise,sunset = SunHours(delta,s,gamma,phi)
    Vals=IntegrateSlope(constant,sunrise,sunset,delta,s,gamma,phi)

    return(Vals)

#------------------------------------------------------------------------------
def TwoPeriodSun(constant,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006
    Calculate the 24-hours extraterrestrial radiation when there are two sun period
    '''
    A1, A2 = SunHours(delta,s,gamma,phi)
    a,b,c = Constants(delta,s,gamma,phi)
    riseSlope, setSlope = BoundsSlope(a,b,c)
    B1 = np.maximum(riseSlope,setSlope)
    B2 = np.minimum(riseSlope,setSlope)
    Angle_B1 = AngleSlope(a,b,c,B1)
    Angle_B2 = AngleSlope(a,b,c,B2)

    B1[abs(Angle_B1) > 0.001] = np.pi - B1[abs(Angle_B1) > 0.001]
    B2[abs(Angle_B2) > 0.001] = -np.pi - B2[abs(Angle_B2) > 0.001]

    # Check if two periods really exist
    ID = np.ravel_multi_index(np.where(np.logical_and(B2 >= A1, B1 >= A2) == True),a.shape)
    Val = IntegrateSlope(constant,B2.flat[ID],B1.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])
    ID = ID[Val < 0]

    # Finally calculate resulting values
    Vals = np.zeros(B1.shape)

    Vals.flat[ID] = (IntegrateSlope(constant,A1.flat[ID],B2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])  +
                   IntegrateSlope(constant,B1.flat[ID],A2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID]))
    ID = np.ravel_multi_index(np.where(Vals == 0),a.shape)
    Vals.flat[ID] = IntegrateSlope(constant,A1.flat[ID],A2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])

    return(Vals)

#------------------------------------------------------------------------------
def IntegrateSlope(constant,sunrise,sunset,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 equation 5
    Calculate the 24 hours extraterrestrial radiation
    '''
    # correct the sunset and sunrise angels for days that have no sunset or no sunrise
    SunOrNoSun = np.logical_or(((np.abs(delta + phi)) > (np.pi/2)),((np.abs(delta - phi)) > (np.pi/2)))
    integral=np.zeros(s.shape)
    ID = np.where(np.ravel(SunOrNoSun==True))

    # No sunset
    if abs(delta+phi.flat[ID])>(np.pi/2):
        sunset1=np.pi
        sunrise1=-np.pi
        integral.flat[ID] = constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(sunset1-sunrise1)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(sunset1-sunrise1)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(sunset1)-np.sin(sunrise1))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(sunset1)-np.sin(sunrise1))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(sunset1)-np.cos(sunrise1)))

    # No sunrise
    elif np.abs(delta-phi.flat[ID])>(np.pi/2):
        integral.flat[ID]=constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(0)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(0)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(0)-np.sin(0))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(0)-np.sin(0))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(0)-np.cos(0)))

    ID = np.where(np.ravel(SunOrNoSun==False))
    integral.flat[ID] = constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(sunset-sunrise)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(sunset-sunrise)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(sunset)-np.sin(sunrise))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(sunset)-np.sin(sunrise))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(sunset)-np.cos(sunrise)))

    return(integral)

#------------------------------------------------------------------------------
def TwoPeriods(delta,s,phi):
    '''
    Based on Richard G. Allen 2006
    Create a boolean map with True values for places with two sunsets
    '''
    TwoPeriods = (np.sin(s) > np.ones(s.shape)*np.sin(phi)*np.sin(delta)+np.cos(phi)*np.cos(delta))

    return(TwoPeriods)

#------------------------------------------------------------------------------
def SunHours(delta,slope,slopedir,lat):
    # Define sun hours in case of one sunlight period

    a,b,c = Constants(delta,slope,slopedir,lat)
    riseSlope, setSlope = BoundsSlope(a,b,c)
    bound = BoundsHorizontal(delta,lat)

    Calculated = np.zeros(slope.shape, dtype = bool)
    RiseFinal = np.zeros(slope.shape)
    SetFinal = np.zeros(slope.shape)

    # First check sunrise is not nan
    # This means that their is either no sunrise (whole day night) or no sunset (whole day light)
    # For whole day light, use the horizontal sunrise and whole day night a zero..
    Angle4 = AngleSlope(a,b,c,-bound)
    RiseFinal[np.logical_and(np.isnan(riseSlope),Angle4 >= 0)] = -bound[np.logical_and(np.isnan(riseSlope),Angle4 >= 0)]
    Calculated[np.isnan(riseSlope)] = True

    # Step 1 > 4
    Angle1 = AngleSlope(a,b,c,riseSlope)
    Angle2 = AngleSlope(a,b,c,-bound)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(Angle2 < Angle1+0.001 ,Angle1 < 0.001),Calculated == False) == True),a.shape)
    RiseFinal.flat[ID] = riseSlope.flat[ID]
    Calculated.flat[ID] = True
    # step 5 > 7
    Angle3 = AngleSlope(a,b,c,-np.pi - riseSlope)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(-bound<(-np.pi-riseSlope),Angle3 <= 0.001),Calculated == False) == True),a.shape)
    RiseFinal.flat[ID] = -np.pi -riseSlope.flat[ID]
    Calculated.flat[ID] = True

    # For all other values we use the horizontal sunset if it is positive, otherwise keep a zero
    RiseFinal[Calculated == False] = -bound[Calculated == False]

    # Then check sunset is not nan or < 0
    Calculated = np.zeros(slope.shape, dtype = bool)

    Angle4 = AngleSlope(a,b,c,bound)
    SetFinal[np.logical_and(np.isnan(setSlope),Angle4 >= 0)] = bound[np.logical_and(np.isnan(setSlope),Angle4 >= 0)]
    Calculated[np.isnan(setSlope)] = True

    # Step 1 > 4
    Angle1 = AngleSlope(a,b,c,setSlope)
    Angle2 = AngleSlope(a,b,c,bound)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(Angle2 < Angle1+0.001,Angle1 < 0.001),Calculated == False) == True),a.shape)
    SetFinal.flat[ID] = setSlope.flat[ID]
    Calculated.flat[ID] = True
    # step 5 > 7
    Angle3 = AngleSlope(a,b,c,np.pi - setSlope)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(bound>(np.pi-setSlope),Angle3 <= 0.001),Calculated == False) == True),a.shape)
    SetFinal.flat[ID] = np.pi - setSlope.flat[ID]
    Calculated.flat[ID] = True

    # For all other values we use the horizontal sunset if it is positive, otherwise keep a zero
    SetFinal[Calculated == False] = bound[Calculated == False]

    #    Angle4 = AngleSlope(a,b,c,bound)
    #    SetFinal[np.logical_and(Calculated == False,Angle4 >= 0)] = bound[np.logical_and(Calculated == False,Angle4 >= 0)]

    # If Sunrise is after Sunset there is no sunlight during the day
    SetFinal[SetFinal <= RiseFinal] = 0
    RiseFinal[SetFinal <= RiseFinal] = 0

    return(RiseFinal,SetFinal)

#------------------------------------------------------------------------------
def Constants(delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 equation 11
    determines constants for calculating the exterrestial solar radiation
    '''
    a = np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma) - np.sin(delta)*np.sin(phi)*np.cos(s)
    b = np.cos(delta)*np.cos(phi)*np.cos(s) + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)
    c = np.cos(delta)*np.sin(s)*np.sin(gamma)

    return(a,b,c)

#------------------------------------------------------------------------------
def BoundsSlope(a,b,c):
    '''
    Based on Richard G. Allen 2006 equation 13
    This function calculates candidate values for sunrise and sunset hour angles
    '''
    Div = (b**2+c**2)
    Div[Div <= 0] = 0.00001
    sinB = (a*c + b*np.sqrt(b**2+c**2-a**2)) / Div
    sinA = (a*c - b*np.sqrt(b**2+c**2-a**2)) / Div

    sinB[sinB < -1] = -1; sinB[sinB > 1] = 1    # Limits see appendix A.2.i
    sinA[sinA < -1] = -1; sinA[sinA > 1] = 1    # Limits see appendix A.2.i

    sunrise = np.arcsin(sinA)
    sunset = np.arcsin(sinB)

    return(sunrise,sunset)

#------------------------------------------------------------------------------
def BoundsHorizontal(delta,phi):
    ''''
    Based on Richard G. Allen 2006
    This function calculates sunrise hours based on earth inclination and latitude
    If there is no sunset or sunrise hours the values are either set to 0 (polar night) or pi (polar day)
    '''
    bound = np.arccos(-np.tan(delta)*np.tan(phi))
    bound[abs(delta+phi) > np.pi/2] = np.pi
    bound[abs(delta-phi) > np.pi/2] = 0

    return(bound)

#------------------------------------------------------------------------------
def AngleSlope(a,b,c,w):
    '''
    Based on Richard G. Allen 2006
    Calculate the cos zenith angle by using the hour angle and constants
    '''
    angle = -a + b*np.cos(w) + c*np.sin(w)

    return(angle)

#------------------------------------------------------------------------------
def Calc_Gradient(dataset,pixel_spacing):
    """
    This function calculates the slope and aspect of a DEM map.
    """
    # constants
    deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
    rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree

    # Calculate slope
    x, y = np.gradient(dataset)
    slope = np.arctan(np.sqrt(np.square(x/pixel_spacing) + np.square(y/pixel_spacing))) * rad2deg

    # calculate aspect
    aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
    aspect = 180 + aspect

    return(deg2rad,rad2deg,slope,aspect)

#------------------------------------------------------------------------------
def DEM_lat_lon(DEM_fileName,output_folder):
    """
    This function retrieves information about the latitude and longitude of the
    DEM map.

    """
    # name for output
    lat_fileName = os.path.join(output_folder, 'Output_radiation_balance','latitude.tif')
    lon_fileName = os.path.join(output_folder, 'Output_radiation_balance','longitude.tif')

    g = gdal.Open(DEM_fileName)     # Open DEM
    geo_t = g.GetGeoTransform()     # Get the Geotransform vector:
    x_size = g.RasterXSize          # Raster xsize - Columns
    y_size = g.RasterYSize          # Raster ysize - Rows

    # create a longitude and a latitude array
    lon = np.zeros((y_size, x_size))
    lat = np.zeros((y_size, x_size))
    for col in np.arange(x_size):
        lon[:, col] = geo_t[0] + col * geo_t[1] + geo_t[1]/2
        # ULx + col*(E-W pixel spacing) + E-W pixel spacing
    for row in np.arange(y_size):
        lat[row, :] = geo_t[3] + row * geo_t[5] + geo_t[5]/2
        # ULy + row*(N-S pixel spacing) + N-S pixel spacing,
        # negative as we will be counting from the UL corner

    # Define shape of the raster
    shape = [x_size, y_size]

    # Save lat and lon files in geo- coordinates
    save_GeoTiff_proy(g, lat, lat_fileName, shape, nband=1)
    save_GeoTiff_proy(g, lon, lon_fileName, shape, nband=1)

    return(lat,lon,lat_fileName,lon_fileName)

#------------------------------------------------------------------------------
def reproject_dataset(dataset, pixel_spacing, UTM_Zone):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection
    """

    # 1) Open the dataset
    g = gdal.Open(dataset)
    if g is None:
        print 'input folder does not exist'

     # Define the EPSG code...
    EPSG_code = '326%02d' % UTM_Zone
    epsg_to = int(EPSG_code)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    try:
        proj = g.GetProjection()
        Proj_in=proj.split('EPSG","')
        epsg_from=int((str(Proj_in[-1]).split(']')[0])[0:-1])
    except:
        epsg_from = int(4326)    # Get the Geotransform vector:
    geo_t = g.GetGeoTransform()

    # Vector components:
    # 0- The Upper Left easting coordinate (i.e., horizontal)
    # 1- The E-W pixel spacing
    # 2- The rotation (0 degrees if image is "North Up")
    # 3- The Upper left northing coordinate (i.e., vertical)
    # 4- The rotation (0 degrees)
    # 5- The N-S pixel spacing, negative as it is counted from the UL corner
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize

    epsg_to = int(epsg_to)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    inProj = Proj(init='epsg:%d' %epsg_from)
    outProj = Proj(init='epsg:%d' %epsg_to)

    nrow_skip = round((0.07*y_size)/2)
    ncol_skip = round((0.04*x_size)/2)

    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the to
    ulx, uly = transform(inProj,outProj,geo_t[0], geo_t[3] + nrow_skip * geo_t[5])
    lrx, lry = transform(inProj,outProj,geo_t[0] + geo_t[1] * (x_size-ncol_skip),
                                        geo_t[3] + geo_t[5] * (y_size-nrow_skip))

    # See how using 27700 and WGS84 introduces a z-value!
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')

    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    col = int((lrx - ulx)/pixel_spacing)
    rows = int((uly - lry)/pixel_spacing)

    # Re-define lr coordinates based on whole number or rows and columns
    (lrx, lry) = (ulx + col * pixel_spacing, uly -
                  rows * pixel_spacing)

    dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)

    if dest is None:
        print 'input folder to large for memory, clip input map'

   # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2], uly,
               geo_t[4], - pixel_spacing)

    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(osng.ExportToWkt())

    # Perform the projection/resampling
    gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(),gdal.GRA_Bilinear)

    return dest, ulx, lry, lrx, uly, epsg_to
#------------------------------------------------------------------------------
def reproject_dataset_example(dataset, dataset_example, method = 1):

    # open example dataset
    g_ex = gdal.Open(dataset_example)
    try:
        proj = g_ex.GetProjection()
        Proj=proj.split('EPSG","')
        epsg_to=int((str(Proj[-1]).split(']')[0])[0:-1])
    except:
        epsg_to = int(4326)

    Y_raster_size = g_ex.RasterYSize
    X_raster_size = g_ex.RasterXSize

    Geo = g_ex.GetGeoTransform()
    ulx = Geo[0]
    uly = Geo[3]
    lrx = ulx + X_raster_size * Geo[1]
    lry = uly + Y_raster_size * Geo[5]

    # open dataset that must be transformed
    g_in = gdal.Open(dataset)
    try:
        proj = g_in.GetProjection()
        Proj=proj.split('EPSG","')
        epsg_from=int((str(Proj[-1]).split(']')[0])[0:-1])
    except:
        epsg_from = int(4326)

    # Set the EPSG codes
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    # Create new raster
    mem_drv = gdal.GetDriverByName('MEM')
    dest1 = mem_drv.Create('', X_raster_size, Y_raster_size, 1, gdal.GDT_Float32)
    dest1.SetGeoTransform(Geo)
    dest1.SetProjection(osng.ExportToWkt())

    # Perform the projection/resampling
    if method == 1:
        gdal.ReprojectImage(g_in, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    if method == 2:
        gdal.ReprojectImage(g_in, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
    if method == 3:
        gdal.ReprojectImage(g_in, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Cubic)

    return(dest1, ulx, lry, lrx, uly, epsg_to)

#------------------------------------------------------------------------------
def save_GeoTiff_proy(src_dataset, dst_dataset_array, dst_fileName, shape_lsc, nband):
    """
    This function saves an array dataset in GeoTiff, using the parameters
    from the source dataset, in projected coordinates

    """
    dst_dataset_array	= np.float_(dst_dataset_array)
    dst_dataset_array[dst_dataset_array<-9999] = np.nan
    geotransform = src_dataset.GetGeoTransform()
    spatialreference = src_dataset.GetProjection()

    # create dataset for output
    fmt = 'GTiff'
    driver = gdal.GetDriverByName(fmt)
    dir_name = os.path.dirname(dst_fileName)

    # If the directory does not exist, make it.
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    dst_dataset = driver.Create(dst_fileName, shape_lsc[0], shape_lsc[1], nband,gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).SetNoDataValue(-9999)
    dst_dataset.GetRasterBand(1).WriteArray(dst_dataset_array)
    dst_dataset = None

#------------------------------------------------------------------------------
def w_time(LT,lon_proy, DOY):
    """
    This function computes the hour angle (radians) of an image given the
    local time, longitude, and day of the year.

    """
    nrow, ncol = lon_proy.shape

    # Difference of the local time (LT) from Greenwich Mean Time (GMT) (hours):
    delta_GTM = np.sign(lon_proy[nrow/2, ncol/2]) * lon_proy[nrow/2, ncol/2] * 24 / 360
    if np.isnan(delta_GTM) == True:
         delta_GTM = np.nanmean(lon_proy) * np.nanmean(lon_proy)  * 24 / 360


    # Local Standard Time Meridian (degrees):
    LSTM = 15 * delta_GTM

    # Ecuation of time (EoT, minutes):
    B = 360./365 * (DOY-81)  # (degrees)
    EoT = 9.87*sin(np.deg2rad(2*B))-7.53*cos(np.deg2rad(B))-1.5*sin(np.deg2rad(B))

    # Net Time Correction Factor (minutes) at the center of the image:
    TC = 4 * (lon_proy - LSTM) + EoT     # Difference in time over the longitude
    LST = LT + delta_GTM + TC/60         # Local solar time (hours)
    HRA = 15 * (LST-12)                  # Hour angle HRA (degrees)
    deg2rad = np.pi / 180.0              # Factor to transform from degree to rad
    w = HRA * deg2rad                    # Hour angle HRA (radians)
    return w

#------------------------------------------------------------------------------
def sensible_heat(rah, ustar, rn_inst, g_inst, ts_dem, ts_dem_hot, ts_dem_cold,
                  air_dens, Surf_temp, k_vk, QC_Map, hot_pixels = None, slope = None, offset_dt = None, slope_dt = None):
    """
    This function computes the instantaneous sensible heat given the
    instantaneous net radiation, ground heat flux, and other parameters.

    """

    #dT_hot_fileName = os.path.join(output_folder, 'Output_cloud_masked','test.tif')
    #save_GeoTiff_proy(dest, dT_hot, dT_hot_fileName,shape, nband=1)

    if (hot_pixels is not None and slope is not None):

        # Near surface temperature difference (dT):
        dT_ini = (rn_inst - g_inst) * rah / (air_dens * 1004)
        dT_hot = np.copy(dT_ini)

        # dT for hot pixels - hot, (dry) agricultural fields with no green veget.:
        dT_hot[ts_dem <= (ts_dem_hot - 0.5)] = np.nan
        dT_hot[QC_Map == 1] = np.nan
        dT_hot[dT_hot == 0] = np.nan
        if np.all(np.isnan(dT_hot)) == True:
           dT_hot = np.copy(dT_ini)
           ts_dem_hot = np.nanpercentile(hot_pixels, 99.5)
           dT_hot[ts_dem <= (ts_dem_hot - 0.5)] = np.nan
           dT_hot[dT_hot == 0] = np.nan

        dT_hot=np.float32(dT_hot)
        dT_hot[slope > 10]=np.nan

        dT_hot_mean = np.nanmean(dT_hot)

        # Compute slope and offset of linear relationship dT = b + a * Ts
        slope_dt = (dT_hot_mean - 0.0) / (ts_dem_hot - ts_dem_cold)  # EThot = 0.0
        offset_dt = dT_hot_mean - slope_dt * ts_dem_hot
    elif (slope_dt is not None and offset_dt is not None):
        pass

    dT = offset_dt + slope_dt * ts_dem

    # Sensible heat flux:
    h = air_dens * 1004 * dT / rah
    h[QC_Map == 1] = np.nan
    h[h==0]=np.nan
    h[QC_Map != 0] = np.nan

    # Monin-Obukhov length (m):
    L_MO = ((-1.0 * air_dens * 1004 * np.power(ustar, 3) * Surf_temp) /
            (k_vk * 9.81 * h))
    L_MO[L_MO < -1000] = -1000

    # Stability correction for momentum, stable conditions (L_MO >= 0):
    psi_200_stable = -0.05 * 200 / L_MO

    # Stability correction for momentum and heat transport, unstable
    # conditions (L_MO < 0):
    x2 = np.power((1.0 - 16.0 * (2.0/L_MO)), 0.25)  # x at 2m
    x200 = np.power(1.0 - 16.0 * (200/L_MO), 0.25)  # x at 200m
    psi_h = 2 * np.log((1 + np.power(x2, 2))/2)
    psi_m200 = (2 * np.log((1 + x200) / 2) + np.log((1 + np.power(x200, 2)) /
                2) - 2 * np.arctan(x200) + 0.5*np.pi)


    return L_MO, psi_200_stable, psi_h, psi_m200, h, dT, slope_dt, offset_dt

#------------------------------------------------------------------------------

def Reshape_Reproject_Input_data(input_File_Name, Example_extend_fileName):

   # Reproject the dataset based on the example
   data_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset_example(
       input_File_Name, Example_extend_fileName)

   # Get the array information from the new created map
   band_data = data_rep.GetRasterBand(1) # Get the reprojected dem band
   ncol_data = data_rep.RasterXSize
   nrow_data = data_rep.RasterYSize

   # Save new dataset
   #stats = band.GetStatistics(0, 1)
   data = band_data.ReadAsArray(0, 0, ncol_data, nrow_data)

   return(data)

#------------------------------------------------------------------------------

def Run_command_window(argument):
    """
    This function runs the argument in the command window without showing cmd window

    Keyword Arguments:
    argument -- string, name of the adf file
    """
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

    process = subprocess.Popen(argument, startupinfo=startupinfo, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()

    return()
