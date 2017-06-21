# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 15:09:49 2016

@author: tih
"""
import numpy as np
import os
import gdal
from math import sin, cos, pi
import re
import subprocess
from openpyxl import load_workbook
import time
import osr
import shutil
from datetime import datetime, timedelta
import glob
import pandas as pd
from pyproj import Proj, transform
import numpy.polynomial.polynomial as poly	

def main():
    ############################## INPUT ##########################################		
        # Input for preSEBAL.py
    for number in range(2,3):                                                                                   # Number defines the column of the inputExcel
        inputExcel=r'J:\SEBAL_Tadla\Excel\InputEXCEL_v3_3_6.xlsx'               # The excel with all the SEBAL input data
        VegetationExcel = r'J:\SEBAL_Tadla\Excel\Vegetation height model.xlsx'  # This excel defines the p and c factor and vegetation height.
        output_folder = r'J:\SEBAL_Tadla\Preprocessing_Output'         # Output folder
        LU_data_FileName = r'J:\SEBAL_Tadla\LandCover\LU_map.tif'        # Path to Land Use map
    
        # optional paramater 
        DSSF_Folder=r'J:\SEBAL_Tadla\LANDSAF'
    
    ######################## Load Excels ##########################################	
        # Open Excel workbook for SEBAL inputs
        wb = load_workbook(inputExcel)
				
        # Open Excel workbook used for Vegetation c and p factor conversions				
        wb_veg = load_workbook(VegetationExcel, data_only=True)			

    ############################## Define General info ############################ 

        Rp = 0.91                        # Path radiance in the 10.4-12.5 µm band (W/m2/sr/µm)
        tau_sky = 0.866                  # Narrow band transmissivity of air, range: [10.4-12.5 µm]
        surf_temp_offset = 3             # Surface temperature offset for water 
        
    ######################## Open General info from SEBAL Excel ################### 

        # Open the General_Input sheet			
        ws = wb['General_Input']
    
        # Extract the input and output folder, and Image type from the excel file			
        input_folder = str(ws['B%d' % number].value)                             
        Image_Type = int(ws['D%d' % number].value)                               # Type of Image (1=Landsat & 2 = VIIRS & GLOBA-V)     
        output_folder = os.path.join(input_folder,'Preprocessing_output')
        output_folder_temp = os.path.join(input_folder,'Preprocessing_output','Temp')

        # Create or empty output folder		
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
								
        NDVI_outfolder = os.path.join(output_folder,'NDVI') 				
        SAVI_outfolder = os.path.join(output_folder,'SAVI') 
        Albedo_outfolder = os.path.join(output_folder,'Albedo') 								
        LAI_outfolder = os.path.join(output_folder,'LAI') 
        Surface_Temperature_Sharp_outfolder = os.path.join(output_folder,'Surface_Temperature_Sharpenend') 
        TRANS_outfolder = os.path.join(output_folder,'Transmissivity') 
        Surface_Temperature_outfolder = os.path.join(output_folder,'Surface_Temperature') 
        
        if not os.path.exists(NDVI_outfolder):
            os.makedirs(NDVI_outfolder)
        if not os.path.exists(SAVI_outfolder):
            os.makedirs(SAVI_outfolder)								
        if not os.path.exists(Albedo_outfolder):
            os.makedirs(Albedo_outfolder)
        if not os.path.exists(LAI_outfolder):
            os.makedirs(LAI_outfolder)
        if not os.path.exists(Surface_Temperature_outfolder):
            os.makedirs(Surface_Temperature_outfolder)
        if not os.path.exists(Surface_Temperature_Sharp_outfolder):
            os.makedirs(Surface_Temperature_Sharp_outfolder)
        if not os.path.exists(TRANS_outfolder):
            os.makedirs(TRANS_outfolder)	
            
        # Extract the Path to the DEM map from the excel file
        DEM_fileName = '%s' %str(ws['E%d' % number].value) #'DEM_HydroShed_m'  

        # Open DEM and create Latitude and longitude files
        lat,lon,lat_fileName,lon_fileName=DEM_lat_lon(DEM_fileName,output_folder_temp)
    
  
    ######################## Extract general data for Landsat ##########################################	
        if Image_Type == 1:	
       
            # Open the Landsat_Input sheet				
            ws = wb['Landsat_Input']		
      
            # Extract Landsat name, number and amount of thermal bands from excel file 
            Name_Landsat_Image = str(ws['B%d' % number].value)    # From glovis.usgs.gov
            Landsat_nr = int(ws['C%d' % number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)
            Bands_thermal = int(ws['D%d' %number].value)         # Number of LS bands to use to retrieve land surface 
 
            # Pixel size of the model
            pixel_spacing=int(30) 
            
            # the path to the MTL file of landsat				
            Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' % Name_Landsat_Image)
            
            # read out the general info out of the MTL file in Greenwich Time
            year, DOY, hour, minutes, UTM_Zone, Sun_elevation = info_general_metadata(Landsat_meta_fileName) # call definition info_general_metadata
            date=datetime.strptime('%s %s'%(year,DOY), '%Y %j')
            month = date.month
            day = date.day
            
            # define the kind of sensor and resolution of the sensor
            sensor1 = 'L%d' % Landsat_nr
            sensor2 = 'L%d' % Landsat_nr
            sensor3 = 'L%d' % Landsat_nr
            res1 = '30m'
            res2 = '%sm' %int(pixel_spacing)
            res3 = '30m'
            
    ######################## Extract general data for VIIRS-PROBAV ##########################################	
        if Image_Type == 2:	
            
            # Open the VIIRS_PROBAV_Input sheet					
            ws = wb['VIIRS_PROBAV_Input']
            
            # Extract the name of the thermal and quality VIIRS image from the excel file	
            Name_VIIRS_Image_TB = '%s' %str(ws['B%d' % number].value)
            
            # Extract the name to the PROBA-V image from the excel file	
            Name_PROBAV_Image = '%s' %str(ws['D%d' % number].value)    # Must be a tiff file 
                 
            # Pixel size of the model
            pixel_spacing=int(100) 	
 		
            # UTM Zone of the end results					
            UTM_Zone = float(ws['G%d' % number].value)

            #Get time from the VIIRS dataset name (IMPORTANT TO KEEP THE TEMPLATE OF THE VIIRS NAME CORRECT example: npp_viirs_i05_20150701_124752_wgs84_fit.tif)
            Total_Day_VIIRS = Name_VIIRS_Image_TB.split('_')[3]
            Total_Time_VIIRS = Name_VIIRS_Image_TB.split('_')[4]
				
            # Get the information out of the VIIRS name in GMT (Greenwich time)
            year = int(Total_Day_VIIRS[1:5])
            month = int(Total_Day_VIIRS[5:7])
            day = int(Total_Day_VIIRS[7:9])
            Startdate = '%d-%02d-%02d' % (year,month,day)
            DOY=datetime.strptime(Startdate,'%Y-%m-%d').timetuple().tm_yday
            hour = int(Total_Time_VIIRS[1:3])
            minutes = int(Total_Time_VIIRS[3:5])
            
            # Rounded difference of the local time from Greenwich (GMT) (hours):
            delta_GTM = round(np.sign(lon[int(np.shape(lon)[0]/2), int(np.shape(lon)[1]/2)]) * lon[int(np.shape(lon)[0]/2), int(np.shape(lon)[1]/2)] * 24 / 360)
            if np.isnan(delta_GTM) == True:
                delta_GTM = round(np.nanmean(lon) * np.nanmean(lon)  * 24 / 360)

            # Calculate local time 
            hour += delta_GTM
            if hour < 0.0:
                day -= 1
                hour += 24
            if hour >= 24:
                day += 1
                hour -= 24         
 
            # define the kind of sensor and resolution of the sensor	
            sensor1 = 'PROBAV'
            sensor2 = 'VIIRS'
            res1 = '375m'
            res2 = '%sm' %int(pixel_spacing)
            res3 = '30m'
        
    ######################## Extract general data from DEM file and create Slope map ##########################################	
        # Variable date name
        Var_name = '%d%02d%02d' %(year, month, day)
    								
        # Reproject from Geog Coord Syst to UTM -
        # 1) DEM - Original DEM coordinates is Geographic: lat, lon
        dest, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(
                       DEM_fileName, pixel_spacing, UTM_Zone=UTM_Zone)
        band = dest.GetRasterBand(1)   # Get the reprojected dem band
        ncol = dest.RasterXSize        # Get the reprojected dem column size
        nrow = dest.RasterYSize        # Get the reprojected dem row size
        shape=[ncol, nrow]
       
        # Read out the DEM band and print the DEM properties
        data_DEM = band.ReadAsArray(0, 0, ncol, nrow)

        # 2) Latitude file - reprojection    
        # reproject latitude to the landsat projection and save as tiff file																																
        lat_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(
                        lat_fileName, pixel_spacing, UTM_Zone=UTM_Zone)
                    
        # Get the reprojected latitude data															
        lat_proy = lat_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)
     
        # 3) Longitude file - reprojection
        # reproject longitude to the landsat projection	 and save as tiff file	
        lon_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset(lon_fileName, pixel_spacing, UTM_Zone=UTM_Zone)

        # Get the reprojected longitude data	
        lon_proy = lon_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)
        
        lon_fileName = os.path.join(output_folder_temp,'lon_resh.tif')
        save_GeoTiff_proy(dest, lon_proy, lon_fileName, shape, nband=1)	
				
        # Calculate slope and aspect from the reprojected DEM
        deg2rad,rad2deg,slope,aspect=Calc_Gradient(data_DEM, pixel_spacing)

        # calculate the coz zenith angle
        Ra_mountain_24, Ra_inst, cos_zn_resh, dr, phi, delta = Calc_Ra_Mountain(lon,DOY,hour,minutes,lon_proy,lat_proy,slope,aspect) 
        cos_zn_fileName = os.path.join(output_folder_temp,'cos_zn.tif')
        save_GeoTiff_proy(dest, cos_zn_resh, cos_zn_fileName, shape, nband=1)
        
        # Save the Ra
        Ra_inst_fileName = os.path.join(output_folder_temp,'Ra_inst.tif')
        save_GeoTiff_proy(dest, Ra_inst, Ra_inst_fileName, shape, nband=1)	
        Ra_mountain_24_fileName = os.path.join(output_folder_temp,'Ra_mountain_24.tif')
        save_GeoTiff_proy(dest, Ra_mountain_24, Ra_mountain_24_fileName, shape, nband=1)	
        
        #################### Calculate Transmissivity ##########################################	
        
        # Open the General_Input sheet			
        ws = wb['Meteo_Input']

        # Extract the method radiation value	
        Value_Method_Radiation_inst = '%s' %str(ws['L%d' % number].value)

        # Values to check if data is created
        Check_Trans_inst = 0
        Check_Trans_24 = 0
    
        # Extract the data to the method of radiation
        if int(Value_Method_Radiation_inst) == 2:
            Field_Radiation_inst = '%s' %str(ws['N%d' % number].value)        
                
            if Field_Radiation_inst == 'None':

                # Create directory for transmissivity
                if not os.path.exists(TRANS_outfolder):
                    os.makedirs(TRANS_outfolder)

                # Instantanious Transmissivity files must be created
                Check_Trans_inst = 1
                
                # Calculate Transmissivity
                quarters_hours = np.ceil(minutes/30.) * 30
                hours_GMT = hour - delta_GTM
                if quarters_hours >= 60:
                    hours_GMT += 1
                    quarters_hours = 0
                
                # Define the instantanious LANDSAF file
                name_Landsaf_inst = 'HDF5_LSASAF_MSG_DSSF_MSG-Disk_%d%02d%02d%02d%02d.tif' %(year, month,day, hours_GMT, quarters_hours)
                file_Landsaf_inst = os.path.join(DSSF_Folder,name_Landsaf_inst)  
                             
                # Reproject the Ra_inst data to match the LANDSAF data
                Ra_inst_3Km_dest = reproject_dataset_example(Ra_inst_fileName, file_Landsaf_inst, method = 2)   
                Ra_inst_3Km = Ra_inst_3Km_dest.GetRasterBand(1).ReadAsArray()
                Ra_inst_3Km[Ra_inst_3Km==0] =np.nan
                
                # Open the Rs LANDSAF data           
                dest_Rs_inst_3Km = gdal.Open(file_Landsaf_inst)
                Rs_inst_3Km = dest_Rs_inst_3Km.GetRasterBand(1).ReadAsArray()
                Rs_inst_3Km = np.float_(Rs_inst_3Km)/10
                Rs_inst_3Km[Rs_inst_3Km<0]=np.nan
                
                # Get shape LANDSAF data
                shape_trans=[dest_Rs_inst_3Km.RasterXSize , dest_Rs_inst_3Km.RasterYSize ]
                
                # Calculate Transmissivity 3Km
                Transmissivity_3Km = Rs_inst_3Km/Ra_inst_3Km
                Transmissivity_3Km_fileName = os.path.join(output_folder_temp,'Transmissivity_3Km.tif')
                save_GeoTiff_proy(Ra_inst_3Km_dest, Transmissivity_3Km, Transmissivity_3Km_fileName, shape_trans, nband=1)	
                
                # Reproject Transmissivity to match DEM (now this is done by using the nearest neighbour method)
                Transmissivity_inst_dest = reproject_dataset_example(Transmissivity_3Km_fileName, cos_zn_fileName)  
                Transmissivity_inst = Transmissivity_inst_dest.GetRasterBand(1).ReadAsArray()
                Transmissivity_inst_fileName = os.path.join(TRANS_outfolder,'Transmissivity_inst_%s.tif' %Var_name)                
                save_GeoTiff_proy(Transmissivity_inst_dest, Transmissivity_inst, Transmissivity_inst_fileName, shape, nband=1)	
                
        # Extract the method radiation value	
        Value_Method_Radiation_24 = '%s' %str(ws['I%d' % number].value)

        # Extract the data to the method of radiation
        if Value_Method_Radiation_24 == 2:
            Field_Radiation_24 = '%s' %str(ws['K%d' % number].value)        
                
            if Field_Radiation_24 == 'None':
                
                # Daily Transmissivity files must be created
                Check_Trans_24 = 1

                # Create directory for transmissivity
                if not os.path.exists(TRANS_outfolder):
                    os.makedirs(TRANS_outfolder)

                # Create times that are needed to calculate daily Rs (LANDSAF)
                Starttime_GMT = datetime.strptime(Startdate,'%Y-%m-%d') + timedelta(hours=-delta_GTM)
                Endtime_GMT = Starttime_GMT + timedelta(days=1)
                Times = pd.date_range(Starttime_GMT, Endtime_GMT,freq = '30min')
 
                for Time in Times[:-1]:
                    year_LANDSAF = Time.year
                    month_LANDSAF = Time.month
                    day_LANDSAF = Time.day
                    hour_LANDSAF = Time.hour
                    min_LANDSAF = Time.minute
 
                    # Define the instantanious LANDSAF file
                    name_Landsaf_inst = 'HDF5_LSASAF_MSG_DSSF_MSG-Disk_%d%02d%02d%02d%02d.tif' %(year_LANDSAF, month_LANDSAF,day_LANDSAF, hour_LANDSAF, min_LANDSAF)
                    file_Landsaf_inst = os.path.join(DSSF_Folder,name_Landsaf_inst)  

                    # Open the Rs LANDSAF data           
                    dest_Rs_inst_3Km = gdal.Open(file_Landsaf_inst)
                    Rs_one_3Km = dest_Rs_inst_3Km.GetRasterBand(1).ReadAsArray()
                    Rs_one_3Km = np.float_(Rs_one_3Km)/10
                    Rs_one_3Km[Rs_one_3Km < 0]=np.nan
                  
                    if Time == Times[0]:
                        Rs_24_3Km_tot = Rs_one_3Km
                    else:
                        Rs_24_3Km_tot += Rs_one_3Km
                
                Rs_24_3Km = Rs_24_3Km_tot / len(Times[:-1])  
                 
                # Reproject the Ra_inst data to match the LANDSAF data
                Ra_24_3Km_dest = reproject_dataset_example(Ra_mountain_24_fileName, file_Landsaf_inst, method = 2)   
                Ra_24_3Km = Ra_24_3Km_dest.GetRasterBand(1).ReadAsArray()
                Ra_24_3Km[Ra_24_3Km==0] = np.nan  
                 
                # Get shape LANDSAF data
                shape_trans=[dest_Rs_inst_3Km.RasterXSize , dest_Rs_inst_3Km.RasterYSize ]
                
                # Calculate Transmissivity 3Km
                Transmissivity_24_3Km = Rs_24_3Km/Ra_24_3Km
            
                Transmissivity_24_3Km_fileName = os.path.join(output_folder_temp,'Transmissivity_24_3Km.tif')
                save_GeoTiff_proy(Ra_24_3Km_dest, Transmissivity_24_3Km, Transmissivity_24_3Km_fileName, shape_trans, nband=1)	
                
                # Reproject Transmissivity to match DEM (now this is done by using the nearest neighbour method)
                Transmissivity_24_dest = reproject_dataset_example(Transmissivity_24_3Km_fileName, cos_zn_fileName)  
                Transmissivity_24 = Transmissivity_24_dest.GetRasterBand(1).ReadAsArray()
                Transmissivity_24_fileName = os.path.join(TRANS_outfolder,'Transmissivity_24_%s.tif' %Var_name)                
                save_GeoTiff_proy(Transmissivity_24_dest, Transmissivity_24, Transmissivity_24_fileName, shape, nband=1)	

    #################### Calculate NDVI and SAVI for LANDSAT ##########################################	

        if Image_Type == 1:	

            # Define bands used for each Landsat number
            if Landsat_nr == 5 or Landsat_nr == 7:
                Bands = np.array([1, 2, 3, 4, 5, 7, 6])
            elif Landsat_nr == 8:
                Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])  
            else:
                print 'Landsat image not supported, use Landsat 7 or 8'

            # Open MTL landsat and get the correction parameters   
            Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' % Name_Landsat_Image)
            Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)
           
            # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
            # for the different Landsat images (L5, L7, or L8)
            ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
            ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
            ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])
    
            # Open one band - To get the metadata of the landsat images only once (to get the extend)
            src_FileName = os.path.join(input_folder, '%s_B2.TIF' % Name_Landsat_Image)  # before 10!
            ls,band_data,ulx,uly,lrx,lry,x_size_ls,y_size_ls = Get_Extend_Landsat(src_FileName)
         
            # Crop the Landsat images to the DEM extent -
            dst_FileName = os.path.join(output_folder_temp,'cropped_LS_b2.tif')  # Before 10 !!
     							
            # Clip the landsat image to match the DEM map											
            lsc = reproject_dataset_example(src_FileName, cos_zn_fileName)	
            data_LS = lsc.GetRasterBand(1).ReadAsArray()
            save_GeoTiff_proy(dest, data_LS, dst_FileName, shape, nband=1)	
    
            # Get the extend of the remaining landsat file after clipping based on the DEM file	
            lsc,band_data,ulx,uly,lrx,lry,x_size_lsc,y_size_lsc = Get_Extend_Landsat(dst_FileName)

            # Create the corrected signals of Landsat in 1 array
            Reflect = Landsat_Reflect(Bands,input_folder,Name_Landsat_Image,output_folder,shape,Lmax,Lmin,ESUN_L5,ESUN_L7,ESUN_L8,cos_zn_resh,dr,Landsat_nr, cos_zn_fileName)

            # Calculate temporal water mask
            water_mask_temp=Water_Mask(shape,Reflect)       

            # Calculate the NDVI and SAVI								
            NDVI,SAVI,albedo = Calc_NDVI_SAVI_albedo(Reflect)

            NDVI_FileName = os.path.join(NDVI_outfolder,'NDVI_LS_%s.tif'%Var_name)
            save_GeoTiff_proy(dest, NDVI, NDVI_FileName, shape, nband=1)
								
            SAVI_FileName = os.path.join(SAVI_outfolder,'SAVI_LS_%s.tif'%Var_name)
            save_GeoTiff_proy(dest, SAVI, SAVI_FileName, shape, nband=1)
								
            albedo_FileName = os.path.join(Albedo_outfolder,'Albedo_LS_%s.tif'%Var_name)
            save_GeoTiff_proy(dest, albedo, albedo_FileName, shape, nband=1)

    ################### Extract Meteo data for Landsat days from SEBAL Excel ##################
        # Open the Meteo_Input sheet	
        ws = wb['Meteo_Input']	
        # ---------------------------- Instantaneous Air Temperature ------------
        # Open meteo data, first try to open as value, otherwise as string (path)	  
        try:
           Temp_inst = float(ws['B%d' %number].value)                # Instantaneous Air Temperature (°C)

        # if the data is not a value, than open as a string	
        except:
            Temp_inst_name = '%s' %str(ws['B%d' %number].value) 
            Temp_inst_fileName = os.path.join(output_folder, 'Temp', 'Temp_inst_input.tif')
            Temp_inst = Reshape_Reproject_Input_data(Temp_inst_name, Temp_inst_fileName, cos_zn_fileName)

        try:
            RH_inst = float(ws['D%d' %number].value)                # Instantaneous Relative humidity (%)
 
        # if the data is not a value, than open as a string							
        except:
            RH_inst_name = '%s' %str(ws['D%d' %number].value) 
            RH_inst_fileName = os.path.join(output_folder, 'Temp', 'RH_inst_input.tif')
            RH_inst = Reshape_Reproject_Input_data(RH_inst_name, RH_inst_fileName, cos_zn_fileName)			
         
        esat_inst = 0.6108 * np.exp(17.27 * Temp_inst / (Temp_inst + 237.3)) 													
        eact_inst = RH_inst * esat_inst / 100
        
     #################### Calculate NDVI and SAVI for VIIRS-PROBAV ##########################################	
     
        if Image_Type == 2:	
    
            # Define the bands that will be used
            bands=['SM', 'B1', 'B2', 'B3', 'B4']  #'SM', 'BLUE', 'RED', 'NIR', 'SWIR'

            # Set the index number at 0
            index=0
							
            # create a zero array with the shape of the reprojected DEM file						
            data_PROBAV=np.zeros((shape[1], shape[0]))
            spectral_reflectance_PROBAV=np.zeros([shape[1], shape[0], 5])
       
            # constants
            n188_float=248       # Now it is 248, but we do not exactly know what this really means and if this is for constant for all images.
   
            # write the data one by one to the spectral_reflectance_PROBAV       
            for bandnmr in bands:

                # Translate the PROBA-V names to the Landsat band names                
                Band_number={'SM':7,'B1':8,'B2':10,'B3':9,'B4':11}

                # Open the .hdf file              
                Band_PROBAVhdf_fileName = os.path.join(input_folder, '%s.HDF5' % (Name_PROBAV_Image))   
                    
                # Open the dataset        
                g=gdal.Open(Band_PROBAVhdf_fileName)
        
                # open the subdataset to get the projection
                sds_b3 = gdal.Open(g.GetSubDatasets()[Band_number[bandnmr]][0])
                Data = sds_b3.GetRasterBand(1).ReadAsArray()
 
                # Define the x and y spacing
                Meta_data = g.GetMetadata()
                Lat_Bottom = float(Meta_data['LEVEL3_GEOMETRY_BOTTOM_LEFT_LATITUDE'])
                Lat_Top = float(Meta_data['LEVEL3_GEOMETRY_TOP_RIGHT_LATITUDE'])
                Lon_Left = float(Meta_data['LEVEL3_GEOMETRY_BOTTOM_LEFT_LONGITUDE'])
                Lon_Right = float(Meta_data['LEVEL3_GEOMETRY_TOP_RIGHT_LONGITUDE'])           
                Pixel_size = float((Meta_data['LEVEL3_GEOMETRY_VNIR_VAA_MAPPING']).split(' ')[-3])
                
                # Define the georeference of the PROBA-V data
                geo_PROBAV=[Lon_Left-0.5*Pixel_size, Pixel_size, 0, Lat_Top+0.5*Pixel_size, 0, -Pixel_size] #0.000992063492063
		 									
                # Define the name of the output file											
                PROBAV_data_name=os.path.join(input_folder, '%s_%s.tif' % (Name_PROBAV_Image,bandnmr)) 									
                dst_fileName=os.path.join(input_folder, PROBAV_data_name)
                
                # create gtiff output with the PROBA-V band
                fmt = 'GTiff'
                driver = gdal.GetDriverByName(fmt)

                dst_dataset = driver.Create(dst_fileName, int(Data.shape[1]), int(Data.shape[0]), 1,gdal.GDT_Float32)
                dst_dataset.SetGeoTransform(geo_PROBAV)
            
                # set the reference info
                srs = osr.SpatialReference()
                srs.SetWellKnownGeogCS("WGS84")
                dst_dataset.SetProjection(srs.ExportToWkt())
           
                # write the array in the geotiff band
                dst_dataset.GetRasterBand(1).WriteArray(Data)
                dst_dataset = None
 
                # Open the PROBA-V band in SEBAL											
                g=gdal.Open(PROBAV_data_name.replace("\\","/")) 
	 										
                # If the data cannot be opened, change the extension											
                if g is None:
                    PROBAV_data_name=os.path.join(input_folder, '%s_%s.tiff' % (Name_PROBAV_Image,bandnmr))  
  
                # Reproject the PROBA-V band  to match DEM's resolution          
                PROBAV = reproject_dataset_example(
                                  PROBAV_data_name, lon_fileName)

                # Open the reprojected PROBA-V band data                         
                data_PROBAV_DN = PROBAV.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)
	 										
                # Define the filename to store the cropped Landsat image
                dst_FileName = os.path.join(output_folder, 'Output_PROBAV','proy_PROBAV_%s.tif' % bandnmr)
		 									
                # close the PROBA-V 											
                g=None
                                   
                # If the band data is not SM change the DN values into PROBA-V values and write into the spectral_reflectance_PROBAV                      
                if bandnmr is not 'SM':  
                    data_PROBAV[:, :]=data_PROBAV_DN/2000                           
                    spectral_reflectance_PROBAV[:, :, index]=data_PROBAV[:, :]
                
                # If the band data is the SM band than write the data into the spectral_reflectance_PROBAV  and create cloud mask            
                else:
                    data_PROBAV[:, :]=data_PROBAV_DN
                    Cloud_Mask_PROBAV=np.zeros((shape[1], shape[0]))
                    Cloud_Mask_PROBAV[data_PROBAV[:,:]!=n188_float]=1
                    spectral_reflectance_PROBAV[:, :, index]=Cloud_Mask_PROBAV
  
                # Change the spectral reflectance to meet certain limits                               
                spectral_reflectance_PROBAV[:, :, index]=np.where(spectral_reflectance_PROBAV[:, :, index]<=0,np.nan,spectral_reflectance_PROBAV[:, :, index])   
                spectral_reflectance_PROBAV[:, :, index]=np.where(spectral_reflectance_PROBAV[:, :, index]>=150,np.nan,spectral_reflectance_PROBAV[:, :, index])   
  										
                # Go to the next index 									
                index=index+1
 
            # Bands in PROBAV spectral reflectance
            # 0 = MS
            # 1 = BLUE
            # 2 = NIR
            # 3 = RED
            # 4 = SWIR
       
            # Calculate SAVI based on PROBA-V
            L = 0.5							
            SAVI = (1+L)*(spectral_reflectance_PROBAV[:, :, 3]-spectral_reflectance_PROBAV[:, :, 2])/(L+spectral_reflectance_PROBAV[:, :, 2]+spectral_reflectance_PROBAV[:, :, 3])
 
            # Calculate surface albedo based on PROBA-V
            Surface_Albedo_PROBAV = 0.219 * spectral_reflectance_PROBAV[:, :, 1] + 0.361 * spectral_reflectance_PROBAV[:, :, 2] + 0.379 * spectral_reflectance_PROBAV[:, :, 3] + 0.041 * spectral_reflectance_PROBAV[:, :, 4]
  
            # Create Water mask based on PROBA-V             
            water_mask_temp = np.zeros((shape[1], shape[0])) 
            water_mask_temp[np.logical_and(spectral_reflectance_PROBAV[:, :, 2] >= spectral_reflectance_PROBAV[:, :, 3],data_DEM>0)]=1

            # Reproject the Veg_height to the LAI projection
            # dest = reproject_dataset3(Veg_Height_proj_FileName, LAI_FileName)			
            Albedo_FileName = os.path.join(Albedo_outfolder,'Albedo_PROBAV_%s.tif' %Var_name) 
								
            save_GeoTiff_proy(dest, Surface_Albedo_PROBAV, Albedo_FileName, shape, nband=1)	  

            # Calculate the NDVI based on PROBA-V     
            n218_memory = spectral_reflectance_PROBAV[:, :, 2] + spectral_reflectance_PROBAV[:, :, 3]
            NDVI = np.zeros((shape[1], shape[0]))
            NDVI[n218_memory != 0] =  ( spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] - spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] )/ ( spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] + spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] )

            NDVI_FileName = os.path.join(NDVI_outfolder,'NDVI_PROBAV_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, NDVI, NDVI_FileName, shape, nband=1)	

            SAVI_FileName = os.path.join(SAVI_outfolder,'SAVI_PROBAV_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, SAVI, SAVI_FileName, shape, nband=1)															
				

     #################### Calculate vegetation parameters ##########################################	

        FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity = Calc_vegt_para(NDVI,SAVI,water_mask_temp,shape)

        if Image_Type == 1:

            LAI_FileName = os.path.join(LAI_outfolder,'LAI_LS_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, LAI, LAI_FileName, shape, nband=1)	

        if Image_Type == 2:

            LAI_FileName = os.path.join(LAI_outfolder,'LAI_PROBAV_%s.tif' %Var_name) 				
            save_GeoTiff_proy(dest, LAI, LAI_FileName, shape, nband=1)	
            
     #################### Calculate thermal for Landsat ##########################################	

        if Image_Type == 1:
								
            therm_data = Landsat_therm_data(Bands,input_folder,Name_Landsat_Image,output_folder,ulx_dem,lry_dem,lrx_dem,uly_dem,shape)          
            Surface_temp=Calc_surface_water_temp(Temp_inst,Landsat_nr,Lmax,Lmin,therm_data,b10_emissivity,k1_c,k2_c,eact_inst,shape,water_mask_temp,Bands_thermal,Rp,tau_sky,surf_temp_offset,Image_Type)
            therm_data_FileName = os.path.join(Surface_Temperature_outfolder,'Surface_Temperature_LS_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, Surface_temp, therm_data_FileName, shape, nband=1)	

     ################################ Apply thermal sharpening #######################################

            # Thermal Sharpening of Thermal data Landsat

            # Upscale DEM to 90m
            pixel_spacing_upscale=90

            dest_90, ulx_dem_90, lry_dem_90, lrx_dem_90, uly_dem_90, epsg_to = reproject_dataset(
                                                      DEM_fileName, pixel_spacing_upscale, UTM_Zone = UTM_Zone)
                                                      
            DEM_90 = dest_90.GetRasterBand(1).ReadAsArray()
            Y_raster_size_90 = dest_90.RasterYSize				
            X_raster_size_90 = dest_90.RasterXSize
            shape_90=([X_raster_size_90, Y_raster_size_90])
						
            proyDEM_fileName_90 = os.path.join(output_folder_temp, 'DEM_90m.tif')      
      
            save_GeoTiff_proy(dest_90, DEM_90, proyDEM_fileName_90, shape_90, nband=1)
	  	
            # save landsat surface temperature
            surf_temp_fileName = os.path.join(output_folder_temp, '%s_%s_surface_temp_%s_%s_%s.tif' %(sensor1, sensor2, res2, year, DOY))
            save_GeoTiff_proy(lsc, Surface_temp, surf_temp_fileName, shape, nband=1)							

            # Upscale NDVI data																																	
            dest_up = reproject_dataset_example(
                                       NDVI_FileName, proyDEM_fileName_90)
  
            NDVI_Landsat_up = dest_up.GetRasterBand(1).ReadAsArray()

            # Upscale Thermal data
            dest_up = reproject_dataset_example(
                                       surf_temp_fileName, proyDEM_fileName_90)
            surface_temp_up = dest_up.GetRasterBand(1).ReadAsArray()
 
            # Define the width of the moving window box
            Box=7	
  
            temp_surface_sharpened_fileName = os.path.join(output_folder_temp,'LS_surface_temp_Sharpened_%s.tif' %Var_name)
  
            # Apply thermal sharpening           
            temp_surface_sharpened = Thermal_Sharpening(surface_temp_up, NDVI_Landsat_up, NDVI, Box, dest_up, output_folder, NDVI_FileName, shape, lsc, temp_surface_sharpened_fileName)	
	  
  
            # Divide temporal watermask in snow and water mask by using surface temperature
            snow_mask, water_mask, ts_moist_veg_min, NDVI_max, NDVI_std = CalculateSnowWaterMask(NDVI,shape,water_mask_temp,temp_surface_sharpened)
       
            # Replace water values by old values before shapening            
            temp_surface_sharpened[water_mask == 1] = Surface_temp[water_mask == 1]
            temp_surface_sharpened = np.where(np.isnan(temp_surface_sharpened),Surface_temp,temp_surface_sharpened)
            
		   # save landsat surface temperature
            surf_temp_fileName = os.path.join(Surface_Temperature_Sharp_outfolder ,'LS_surface_temp_sharpened_%s.tif' %Var_name)
            save_GeoTiff_proy(lsc, temp_surface_sharpened, surf_temp_fileName, shape, nband=1)	


    ################################## Calculate VIIRS surface temperature ########################

        if Image_Type == 2:

            # Define the VIIRS thermal data name
            VIIRS_data_name=os.path.join(input_folder, '%s' % (Name_VIIRS_Image_TB))
							
            # Reproject VIIRS thermal data								
            VIIRS = reproject_dataset_example(VIIRS_data_name, LAI_FileName)
																
            # Open VIIRS thermal data																		
            data_VIIRS = VIIRS.GetRasterBand(1).ReadAsArray()    
				
            # Define the thermal VIIRS output name
            proyVIIRS_fileName = os.path.join(output_folder_temp, 'Surface_Temp_VIIRS_%s.tif' %Var_name)
	 											
            # Save the thermal VIIRS data 												
            save_GeoTiff_proy(dest, data_VIIRS, proyVIIRS_fileName, shape, nband=1)	

            # Set the conditions for the brightness temperature (100m)
            brightness_temp=np.where(data_VIIRS>=250, data_VIIRS, np.nan)
 
            # Constants
            k1=606.399172
            k2=1258.78
            L_lambda_b10_100=((2*6.63e-34*(3.0e8)**2)/((11.45e-6)**5*(np.exp((6.63e-34*3e8)/(1.38e-23*(11.45e-6)*brightness_temp))-1)))*1e-6
         
            # Get Temperature for 100 and 375m resolution
            Temp_TOA_100 = Get_Thermal(L_lambda_b10_100,Rp,Temp_inst,tau_sky,tir_emis,k1,k2) 
       
            # Conditions for surface temperature (100m)
            n120_surface_temp=Temp_TOA_100.clip(250, 450)

            # Save the surface temperature of the VIIRS in 100m resolution
            temp_surface_100_fileName_beforeTS = os.path.join(Surface_Temperature_outfolder,'Surface_Temperature_VIIRS_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, n120_surface_temp, temp_surface_100_fileName_beforeTS, shape, nband=1)     
        
        
            print '---------------------------------------------------------'
            print '-------------------- Downscale VIIRS --------------------'
            print '---------------------------------------------------------'

        ################################ Thermal Sharpening #####################################################

            # Upscale DEM to 90m
            pixel_spacing_upscale=400

            dest_400, ulx_dem_400, lry_dem_400, lrx_dem_400, uly_dem_400, epsg_to = reproject_dataset(
                                                      DEM_fileName, pixel_spacing_upscale, UTM_Zone = UTM_Zone)
                                                      
            DEM_400 = dest_400.GetRasterBand(1).ReadAsArray()
            Y_raster_size_400 = dest_400.RasterYSize				
            X_raster_size_400 = dest_400.RasterXSize
            shape_400=([X_raster_size_400, Y_raster_size_400])
						
            proyDEM_fileName_400 = os.path.join(output_folder_temp, 'DEM_400m.tif')      
      
            save_GeoTiff_proy(dest_400, DEM_400, proyDEM_fileName_400, shape_400, nband=1)
    																		
            # Upscale thermal band VIIRS from 100m to 400m
            VIIRS_Upscale = reproject_dataset_example(temp_surface_100_fileName_beforeTS, proyDEM_fileName_400)
            data_Temp_Surf_400 = VIIRS_Upscale.GetRasterBand(1).ReadAsArray()
 
            # Upscale PROBA-V NDVI from 100m to 400m       
            NDVI_PROBAV_Upscale = reproject_dataset_example(
                              NDVI_FileName, proyDEM_fileName_400)
            data_NDVI_400 = NDVI_PROBAV_Upscale.GetRasterBand(1).ReadAsArray()

            # Define the width of the moving window box
            Box=9
           
            # The surface temperature temporary file name
            temp_surface_sharpened_fileName = os.path.join(output_folder_temp,'VIIRS_surface_temp_Sharpened_%s.tif' %Var_name)
  
            # Apply the surface temperature sharpening        
            temp_surface_sharpened = Thermal_Sharpening(data_Temp_Surf_400, data_NDVI_400, NDVI, Box, NDVI_PROBAV_Upscale, output_folder, NDVI_FileName, shape, dest, temp_surface_sharpened_fileName)	
      
            # Divide temporal watermask in snow and water mask by using surface temperature
            snow_mask, water_mask, ts_moist_veg_min, NDVI_max, NDVI_std = CalculateSnowWaterMask(NDVI,shape,water_mask_temp,temp_surface_sharpened)
       
            # Replace water values by old values before shapening            
            temp_surface_sharpened[water_mask == 1] = n120_surface_temp[water_mask == 1]
            temp_surface_sharpened = np.where(np.isnan(temp_surface_sharpened),n120_surface_temp,temp_surface_sharpened)  
            
		   # save landsat surface temperature
            surf_temp_fileName = os.path.join(Surface_Temperature_Sharp_outfolder ,'VIIRS_surface_temp_sharpened_%s.tif' %Var_name)
            save_GeoTiff_proy(dest, temp_surface_sharpened, surf_temp_fileName, shape, nband=1)	        

    ################################## Set path to output maps to excel file ########################

        # Open meteo sheet
        xfile = load_workbook(inputExcel)
        sheet = xfile.get_sheet_by_name('Meteo_Input')

        # If instantanious Transmissivity is calculated in PreSEBAL        
        if Check_Trans_inst == 1:

            sheet['N%d'%(number)] = str(Transmissivity_inst_fileName)
            xfile.save(inputExcel)

        # If daily Transmissivity is calculated in PreSEBAL                
        if Check_Trans_24 == 1:  

            sheet['K%d'%(number)] = str(Transmissivity_24_fileName)
            xfile.save(inputExcel)            
            
################################################### HANTS #######################################################





''' 	

################################## All input is now calculated, so preprosessing can start ########################

    # Open preprosessing excel the Vegetation_Height sheet				
    ws_veg = wb_veg['Vegetation_Height'] 

    # Define output name for the LandUse map 
    dst_FileName = os.path.join(output_folder,'LU_%s.tif' %Var_name) 

    # Open LU data
    LU_dest = gdal.Open(LU_data_FileName)
    LU_data = LU_dest.GetRasterBand(1).ReadAsArray() 
 
    # Reproject the LAI to the same projection as LU
    dest1 = reproject_dataset_example(LAI_FileName, LU_data_FileName)	 ## input after HANTS
    LAI_proj = dest1.GetRasterBand(1).ReadAsArray() 

    # Read out the excel file coefficient numbers			
    Array = np.zeros([ws_veg.max_row-1,4])
    for j in ['A','C','D','E']:
        j_number={'A' : 0, 'C' : 1, 'D' : 2, 'E' : 3}					
        for i in range(2,ws_veg.max_row+1):											
	        Value = (ws_veg['%s%s' %(j,i)].value)  																
	        Array[i-2, j_number[j]] = Value										

    # Create maps with the coefficient numbers for the right land cover
    coeff = np.zeros([int(np.shape(LU_data)[0]),int(np.shape(LU_data)[1]),3])
    for coeff_nmbr in range(0,3):				
        for Class in range(0,len(Array)):
	        coeff[LU_data==Array[Class,0],coeff_nmbr] = Array[Class,coeff_nmbr+1]

    # Get some dimensions of the projected dataset 
    band_data = dest1.GetRasterBand(1) 
    ncol_data = dest1.RasterXSize
    nrow_data = dest1.RasterYSize
    shape_data=[ncol_data, nrow_data]

    # Calculate the vegetation height in the LU projection
    Veg_Height_proj = coeff[:,:,0] * np.power(LAI_proj,2) + coeff[:,:,1] * LAI_proj + coeff[:,:,2]
    Veg_Height_proj = np.clip(Veg_Height_proj, 0, 600)

    # Save the vegetation height in the lU projection in the temporary directory
    Veg_Height_proj_FileName = os.path.join(output_folder_temp,'Veg_Height_proj.tif') 				
    save_GeoTiff_proy(dest1, Veg_Height_proj, Veg_Height_proj_FileName, shape_data, nband=1)	
				
    # Reproject the Veg_height to the LAI projection
    dest = reproject_dataset_example(Veg_Height_proj_FileName, LAI_FileName)			

    # Get some dimensions of the original dataset 
    band_data = dest.GetRasterBand(1)
    ncol_data = dest.RasterXSize
    nrow_data = dest.RasterYSize

    # Open the Veg_height with the same projection as LAI				
    Veg_Height = band_data.ReadAsArray(0, 0, ncol_data, nrow_data)
    Veg_Height[Veg_Height == 0] = np.nan				

    # Save Vegetation Height in the end folder				
    dst_FileName = os.path.join(output_folder,'Vegetation_Height_%s.tif' %Var_name) 	
    save_GeoTiff_proy(dest, Veg_Height, dst_FileName, shape, nband=1)			

######################## calculate p-factor by using the Landuse map #########################
    ws_p = wb_veg['p-factor'] 
			
    Array_P = np.zeros([ws_p.max_row-1,2])
    for j in ['A','C']:
        j_number={'A' : 0, 'C' : 1}					
        for i in range(2,ws_p.max_row+1):											
            Value = (ws_p['%s%s' %(j,i)].value)  																
            Array_P[i-2, j_number[j]] = Value	

    p_factor = np.zeros([int(np.shape(LU_data)[0]),int(np.shape(LU_data)[1])])		
    for Class in range(0,len(Array_P)):
	    p_factor[LU_data==Array_P[Class,0]] = Array_P[Class,1]

    p_factor[p_factor == 0] = np.nan

    dst_FileName = os.path.join(output_folder_temp,'p-factor_proj.tif') 	
    save_GeoTiff_proy(dest1, p_factor, dst_FileName, shape_data, nband=1)

    dest = reproject_dataset_example(dst_FileName, LAI_FileName)	

    band_data = dest.GetRasterBand(1) # Get the reprojected dem band	
    ncol_data = dest.RasterXSize
    nrow_data = dest.RasterYSize				
    p_factor = band_data.ReadAsArray(0, 0, ncol_data, nrow_data)
    p_factor[p_factor == 0] = np.nan

    dst_pfactor_FileName = os.path.join(output_folder,'p-factor_%s.tif' %Var_name) 
    save_GeoTiff_proy(dest, p_factor, dst_pfactor_FileName, shape, nband=1)	

######################## calculate c-factor by using the Landuse map #########################

    ws_c = wb_veg['C-factor'] 
			
    Array_C = np.zeros([ws_c.max_row-1,2])
    for j in ['A','C']:
        j_number={'A' : 0, 'C' : 1}					
        for i in range(2,ws_c.max_row+1):											
            Value = (ws_c['%s%s' %(j,i)].value)  																
            Array_C[i-2, j_number[j]] = Value	

    c_factor = np.zeros([int(np.shape(LU_data)[0]),int(np.shape(LU_data)[1])])		
    for Class in range(0,len(Array_C)):
	    c_factor[LU_data==Array_C[Class,0]] = Array_C[Class,1]

    c_factor[np.logical_and(c_factor != 3.0, c_factor != 4.0)] = np.nan

    LUE_max = np.zeros([int(np.shape(LU_data)[0]),int(np.shape(LU_data)[1])])	
    LUE_max[c_factor == 3] = 2.5
    LUE_max[c_factor == 4] = 4.5
    LUE_max[LUE_max == 0] = np.nan

    dst_FileName = os.path.join(output_folder_temp,'LUE_max_proj.tif') 	
    save_GeoTiff_proy(dest1, LUE_max, dst_FileName, shape_data, nband=1)

    dest = reproject_dataset_example(dst_FileName, LAI_FileName)	

    band_data = dest.GetRasterBand(1) # Get the reprojected dem band	
    ncol_data = dest.RasterXSize
    nrow_data = dest.RasterYSize				
    LUE_max = band_data.ReadAsArray(0, 0, ncol_data, nrow_data)
    LUE_max[LUE_max == 0] = np.nan

    dst_LUEmax_FileName = os.path.join(output_folder,'LUE_max_%s.tif' %Var_name) 
    save_GeoTiff_proy(dest, LUE_max, dst_LUEmax_FileName, shape, nband=1)	

############################# delete temporary directory ########################
    shutil.rmtree(output_folder_temp)
#################################################################################
''' 			
				
				
# Functions
#################################################################################   
def DEM_lat_lon(DEM_fileName,output_folder_temp):
    """
    This function retrieves information about the latitude and longitude of the
    DEM map. 
    
    """
    # name for output
    lat_fileName = os.path.join(output_folder_temp,'latitude.tif')
    lon_fileName = os.path.join(output_folder_temp,'longitude.tif')
				
				
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

    # Save lat and lon files in geo- coordinates
    save_GeoTiff_geo(g, lat, lat_fileName, x_size, y_size, nband=1)
    save_GeoTiff_geo(g, lon, lon_fileName, x_size, y_size, nband=1)
    
    return(lat,lon,lat_fileName,lon_fileName)

#------------------------------------------------------------------------------
def Calc_Gradient(dataset,pixel_spacing):
    """
    This function calculates the slope and aspect of a DEM map.
    """
    # constants
    deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
    rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree
    
    # calulate slope from DEM map
    x, y = np.gradient(dataset)
    slope = np.arctan(np.sqrt(np.square(x/pixel_spacing) + np.square(y/pixel_spacing))) * rad2deg
    
    # calculate aspect                  
    aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
    aspect = 180 + aspect

    return(deg2rad,rad2deg,slope,aspect)
  
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
    Loc_time = float(hour) + float(minutes)/60  # Local time (hours)
    
    # 1. Calculation of extraterrestrial solar radiation for slope and aspect
    # Computation of Hour Angle (HRA = w)
    B = 360./365 * (DOY-81)           # (degrees)
    # Computation of cos(theta), where theta is the solar incidence angle
    # relative to the normal to the land surface
    delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)
    phi = lat_proy * deg2rad                                     # latitude of the pixel (radians)
    s = slope * deg2rad                                          # Surface slope (radians)
    gamma = (aspect-180) * deg2rad                               # Surface aspect angle (radians)
    w = w_time(Loc_time, lon_proy, DOY)                            # Hour angle (radians)
    a,b,c = Constants(delta,s,gamma,phi)
    cos_zn= AngleSlope(a,b,c,w)
    cos_zn = cos_zn.clip(Min_cos_zn, Max_cos_zn)

    print '  Local time: ', '%0.3f' % Loc_time 
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
    return(Ra_mountain_24, Ra_inst, cos_zn, dr, phi, delta)

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
def OnePeriodSun(constant,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006
    Calculate the 24-hours extraterrestrial radiation when there is one sun period
    '''
    sunrise,sunset = SunHours(delta,s,gamma,phi)    
    Vals=IntegrateSlope(constant,sunrise,sunset,delta,s,gamma,phi)
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
def SunHours(delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 Appendix A
    refines integration limits (hour angels) so that they are correctly applied
    for all combinations of slope, aspect, and latitude
    '''
    a,b,c = Constants(delta,s,gamma,phi)
    riseSlope, setSlope = BoundsSlope(a,b,c)
    bound = BoundsHorizontal(delta,phi)
    
    Calculated = np.zeros(s.shape, dtype = bool)    
    RiseFinal = np.zeros(s.shape)    
    SetFinal = np.zeros(s.shape)

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
    
    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(-bound<(-np.pi-riseSlope),Angle3 <= 0.001),Calculated == False) == True),a.shape) # Limits see appendix A.2.iv 
    RiseFinal.flat[ID] = -np.pi -riseSlope.flat[ID] # Limits see appendix A.2.v
    Calculated.flat[ID] = True
    
    # For all other values we use the horizontal sunset if it is positive, otherwise keep a zero   
    RiseFinal[Calculated == False] = -bound[Calculated == False]      
        
    # Then check sunset is not nan or < 0 
    Calculated = np.zeros(s.shape, dtype = bool)     

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
def TwoPeriods(delta,s,phi):
    '''
    Based on Richard G. Allen 2006
    Create a boolean map with True values for places with two sunsets
    '''
    TwoPeriods = (np.sin(s) > np.ones(s.shape)*np.sin(phi)*np.sin(delta)+np.cos(phi)*np.cos(delta))
    return(TwoPeriods)  

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
def AngleSlope(a,b,c,w):
    '''
    Based on Richard G. Allen 2006
    Calculate the cos zenith angle by using the hour angle and constants
    '''
    angle = -a + b*np.cos(w) + c*np.sin(w)
    return(angle)    

#------------------------------------------------------------------------------
def Landsat_therm_data(Bands,input_folder,Name_Landsat_Image,output_folder,ulx_dem,lry_dem,lrx_dem,uly_dem,shape_lsc):          
    """
    This function calculates and returns the thermal data from the landsat image.
    """                             
    therm_data = np.zeros((shape_lsc[1], shape_lsc[0], len(Bands)-6))
    for band in Bands[-(len(Bands)-6):]:
        # Open original Landsat image for the band number
        src_FileName = os.path.join(input_folder, '%s_B%1d.TIF'
                                    % (Name_Landsat_Image, band))
        if not os.path.exists(src_FileName):
             src_FileName = os.path.join(input_folder, '%s_B%1d_VCID_2.TIF'
                                    % (Name_Landsat_Image, band))		    																																
																																				
        # Define the filename to store the cropped Landsat image
        dst_FileName = os.path.join(output_folder, 'Temp',
                                    'cropped_LS_b%1d.tif' % band)

        ls_data=Open_landsat(src_FileName,dst_FileName,ulx_dem,lry_dem,lrx_dem,uly_dem,shape_lsc) 																																	

        index = np.where(Bands[:] == band)[0][0] - 6
        therm_data[:, :, index] = ls_data
								
    return(therm_data)

#------------------------------------------------------------------------------    
def Landsat_Reflect(Bands,input_folder,Name_Landsat_Image,output_folder,shape_lsc,Lmax,Lmin,ESUN_L5,ESUN_L7,ESUN_L8,cos_zn_resh,dr,Landsat_nr, proyDEM_fileName):
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
        # stats = band_data.GetStatistics(0, 1)

        index = np.where(Bands[:-(len(Bands)-6)] == band)[0][0]
        if Landsat_nr == 8:
            # Spectral radiance for each band:
            L_lambda = Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L8, index, cos_zn_resh, dr)
        elif Landsat_nr == 7:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L7, index, cos_zn_resh, dr)
        elif Landsat_nr == 5:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda =Landsat_rho_lambda(L_lambda, ESUN_L5, index, cos_zn_resh, dr)
        else:
            print 'Landsat image not supported, use Landsat 5, 7 or 8'

        Spec_Rad[:, :, index] = L_lambda
        Reflect[:, :, index] = rho_lambda
    Reflect = Reflect.clip(0.0, 1.0)
    return(Reflect,Spec_Rad)
     

#------------------------------------------------------------------------------   
def Open_landsat(src_FileName, proyDEM_fileName):
    """
    This function opens a landsat image and returns the data array of a specific landsat band.
    """                           
    # crop band to the DEM extent
    ls, ulx, uly, lrx, lry, epsg_to = reproject_dataset_example(src_FileName, proyDEM_fileName)											
    
    # Open the cropped Landsat image for the band number
    ls_data = ls.GetRasterBand(1).ReadAsArray()
    return(ls_data) 

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
def Landsat_rho_lambda(L_lambda,ESUN,index,cos_zn_resh,dr):
    """
    Calculates the lambda from landsat
    """
    rho_lambda = np.pi * L_lambda / (ESUN[index] * cos_zn_resh * dr)
    return(rho_lambda)
   				
#------------------------------------------------------------------------------
def Calc_NDVI_SAVI_albedo(Reflect):
    """
    This function calculates and returns the Surface albedo, NDVI, and SAVI by using the refectance from the landsat image.
    """
    # Computation of Normalized Difference Vegetation Index (NDVI)
    # and Soil Adjusted Vegetation Index (SAVI):
    L = 0.5				
    NDVI = ((Reflect[:, :, 3] - Reflect[:, :, 2]) /
            (Reflect[:, :, 3] + Reflect[:, :, 2]))
    SAVI = (1 + L) * ((Reflect[:, :, 3] - Reflect[:, :, 2]) /
                      (L + Reflect[:, :, 3] + Reflect[:, :, 2]))

    Apparent_atmosf_transm = 0.89    # This value is used for atmospheric correction of broad band albedo. This value is used for now, would be better to use tsw.
    path_radiance = 0.03             # Recommended, Range: [0.025 - 0.04], based on Bastiaanssen (2000).
																						
    # Surface albedo:
    Surf_albedo = (0.254 * Reflect[:, :, 0] + 0.149 * Reflect[:, :, 1] +
                   0.147 * Reflect[:, :, 2] + 0.311 * Reflect[:, :, 3] +
                   0.103 * Reflect[:, :, 4] + 0.036 * Reflect[:, :, 5] -
                   path_radiance) / np.power(Apparent_atmosf_transm, 2)

    # Better tsw instead of Apparent_atmosf_transm ??
    Surf_albedo = Surf_albedo.clip(0.0, 0.6)																						
																						
    return(NDVI,SAVI,Surf_albedo)

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

    return(Surface_temp)    

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
def Calc_vegt_para(NDVI,SAVI,water_mask_temp,shape_lsc):
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
    VI_NDVI = 38.764 * np.square(NDVI) - 24.605 * NDVI + 5.8103
    VI_SAVI = 6.3707 * np.square(SAVI) - 2.8503 * SAVI + 1.6335
    VI = (VI_NDVI + VI_SAVI) / 2.0  # Average of computed from NDVI and SAVI

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
    LAI_3 = 11.0 * np.power(SAVI, 3)
    LAI_3[SAVI >= 0.817] = 6.0
    LAI_4 = -np.log((0.69 - SAVI) / 0.59) / 0.91  # For South. Idaho, empirical
    LAI_4[SAVI < 0.0] = 0.0
    LAI_4[SAVI >= 0.689] = 6.0

    LAI = (LAI_1 + LAI_2 + LAI_3 + LAI_4) / 4.0  # Average LAI
    LAI[LAI < 0.001] = 0.001

    b10_emissivity = np.zeros((shape_lsc[1], shape_lsc[0]))
    b10_emissivity = np.where(LAI <= 3.0, 0.95 + 0.01 * LAI, 0.98)
    b10_emissivity[water_mask_temp != 0.0] = 1.0
    return(FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity)

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

    epsg_from = Get_epsg(g)	
   
    # Get the Geotransform vector:
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
    DEM_UL_lat = geo_t[3]

    # Define the EPSG code...
    if DEM_UL_lat > 0:
        EPSG_code = '326%02d' % UTM_Zone
    else:
        EPSG_code = '326%02d' % UTM_Zone
        UTM_Zone = - UTM_Zone
    epsg_to = int(EPSG_code)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    inProj = Proj(init='epsg:%d' %epsg_from)
    outProj = Proj(init='epsg:%d' %epsg_to)
				
    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the to

    # 3) Work out the boundaries of the new dataset in the target projection
    #   Skip some rows and columns in the border to avoid null values due to
    #   reprojection - rectangle to parallelogram
    nrow_skip = round((0.10*y_size)/2)
    ncol_skip = round((0.10*x_size)/2)
    

    ulx, uly = transform(inProj,outProj,geo_t[0]+ncol_skip*geo_t[1], geo_t[3] +
                       nrow_skip * geo_t[5])
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
    (ulx, uly) = (int(ulx), int(uly))
    (lrx, lry) = (int(ulx) + col * pixel_spacing, int(uly) -
                  rows * pixel_spacing)
    dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    if dest is None:
        print 'input folder to large for memory, clip input map'
     
   # Calculate the new geotransform
    new_geo = (int(ulx), pixel_spacing, geo_t[2], int(uly),
               geo_t[4], - pixel_spacing)
    
    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(osng.ExportToWkt())
      
    # Perform the projection/resampling
    gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(),gdal.GRA_Bilinear)
				
    return dest, ulx, lry, lrx, uly, epsg_to
				
#------------------------------------------------------------------------------
def reproject_dataset_example(dataset, dataset_example,method = 1):

    # open dataset that must be transformed    
    g = gdal.Open(dataset)
    epsg_from = Get_epsg(g)	   

    # open dataset that is used for transforming the dataset
    gland=gdal.Open(dataset_example) 
    epsg_to = Get_epsg(gland)	

    # Set the EPSG codes
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    # Get shape and geo transform from example				
    geo_land = gland.GetGeoTransform()			
    col=gland.RasterXSize
    rows=gland.RasterYSize

    # Create new raster			
    mem_drv = gdal.GetDriverByName('MEM')
    dest1 = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    dest1.SetGeoTransform(geo_land)
    dest1.SetProjection(osng.ExportToWkt())
    
    # Perform the projection/resampling
    if method == 1:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    if method == 2:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)

    return(dest1)		
				
#------------------------------------------------------------------------------
def Get_epsg(g):				
			
    try:
        # Get info of the dataset that is used for transforming     
        gland_proj = g.GetProjection()
        Projection=gland_proj.split('EPSG","')
        epsg_to=int((str(Projection[-1]).split(']')[0])[0:-1])				      
    except:
        epsg_to=4326	
        print 'Was not able to get the projection, so WGS84 is assumed'							
    return(epsg_to)	
				
#------------------------------------------------------------------------------				
def Get_Extend_Landsat(src_FileName):
    """
    This function gets the extend of the landsat image
    """
    ls = gdal.Open(src_FileName)       # Open Landsat image
    print 'Original LANDSAT Image - '
    geo_t_ls = ls.GetGeoTransform()    # Get the Geotransform vector
    x_size_ls = ls.RasterXSize         # Raster xsize - Columns
    y_size_ls = ls.RasterYSize         # Raster ysize - Rows
    print '  Size :', x_size_ls, y_size_ls
    (ulx, uly) = geo_t_ls[0], geo_t_ls[3]
    (lrx, lry) = (geo_t_ls[0] + geo_t_ls[1] * x_size_ls,
                  geo_t_ls[3] + geo_t_ls[5] * y_size_ls)
    band_data = ls.GetRasterBand(1)
    
    return(ls,band_data,ulx,uly,lrx,lry,x_size_ls,y_size_ls)
    				
#------------------------------------------------------------------------------
def save_GeoTiff_geo(src_dataset, dst_dataset_array, dst_fileName, ncol, nrow,
                     nband):
    """
    This function saves an array dataset in GeoTiff, using the parameters
    from the source dataset, in geographical coordinates

    """
    geotransform = src_dataset.GetGeoTransform()
    # create dataset for output
    fmt = 'GTiff'
    driver = gdal.GetDriverByName(fmt)
    dir_name = os.path.dirname(dst_fileName)
    # If the directory does not exist, make it.
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    dst_dataset = driver.Create(dst_fileName, ncol, nrow, nband,gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    # set the reference info
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84")
    dst_dataset.SetProjection(srs.ExportToWkt())
    # write the array in the geotiff band
    dst_dataset.GetRasterBand(1).WriteArray(dst_dataset_array)
    # stats = dst_dataset.GetRasterBand(1).GetStatistics(0, 1)
    dst_dataset = None				
				
#------------------------------------------------------------------------------
def Reshape_Reproject_Input_data(input_File_Name, output_File_Name, Example_extend_fileName):
       
   data_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = reproject_dataset_example(
       input_File_Name, Example_extend_fileName)
   band_data = data_rep.GetRasterBand(1) # Get the reprojected dem band
   ncol_data = data_rep.RasterXSize
   nrow_data = data_rep.RasterYSize
   shape_data=[ncol_data, nrow_data]
   
   #stats = band.GetStatistics(0, 1)
   data = band_data.ReadAsArray(0, 0, ncol_data, nrow_data)
   save_GeoTiff_proy(data_rep, data, output_File_Name, shape_data, nband=1)
   return(data)

#------------------------------------------------------------------------------				
def save_GeoTiff_proy(src_dataset, dst_dataset_array, dst_fileName, shape_lsc,nband):
    """
    This function saves an array dataset in GeoTiff, using the parameters
    from the source dataset, in projected coordinates

    """
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
    dst_dataset.GetRasterBand(1).WriteArray(dst_dataset_array)
    dst_dataset = None
    
#------------------------------------------------------------------------------   
def Thermal_Sharpening(surface_temp_up, NDVI_up, NDVI, Box, dest_up, output_folder, ndvi_fileName, shape_down, dest_down, temp_surface_sharpened_fileName):

    # Creating arrays to store the coefficients						
    CoefA=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
    CoefB=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
    CoefC=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
  
    # Fit a second polynominal fit to the NDVI and Thermal data and save the coeffiecents for each pixel  
    # NOW USING FOR LOOPS PROBABLY NOT THE FASTEST METHOD     
    for i in range(0,len(surface_temp_up)):
        for j in range(0,len(surface_temp_up[1])):
            if np.isnan(np.sum(surface_temp_up[i,j]))==False and np.isnan(np.sum(NDVI_up[i,j]))==False:
                x_data=NDVI_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)][np.logical_and(np.logical_not(np.isnan(NDVI_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)])),np.logical_not(np.isnan(surface_temp_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)])))]
                y_data=surface_temp_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)][np.logical_and(np.logical_not(np.isnan(NDVI_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)])),np.logical_not(np.isnan(surface_temp_up[np.maximum(0,i-(Box-1)/2):np.minimum(len(surface_temp_up),i+(Box-1)/2+1),np.maximum(0,j-(Box-1)/2):np.minimum(len(surface_temp_up[1]),j+(Box-1)/2+1)])))]
                if len(x_data)>6:
                    coefs = poly.polyfit(x_data, y_data, 2)
                    CoefA[i,j] = coefs[2]
                    CoefB[i,j] = coefs[1]
                    CoefC[i,j] = coefs[0]
                else:
                    CoefA[i,j] = np.nan
                    CoefB[i,j] = np.nan
                    CoefC[i,j] = np.nan
            else:
                CoefA[i,j] = np.nan
                CoefB[i,j] = np.nan
                CoefC[i,j] = np.nan
            
    # Define the shape of the surface temperature with the resolution of 400m      
    shape_up=[len(surface_temp_up[1]),len(surface_temp_up)]
           
    # Save the coefficients
    CoefA_fileName_Optie2 = os.path.join(output_folder, 'Output_temporary','coef_A.tif')
    save_GeoTiff_proy(dest_up,CoefA, CoefA_fileName_Optie2,shape_up, nband=1)
        
    CoefB_fileName_Optie2 = os.path.join(output_folder, 'Output_temporary','coef_B.tif')
    save_GeoTiff_proy(dest_up,CoefB, CoefB_fileName_Optie2,shape_up, nband=1)
 
    CoefC_fileName_Optie2 = os.path.join(output_folder, 'Output_temporary','coef_C.tif')
    save_GeoTiff_proy(dest_up,CoefC, CoefC_fileName_Optie2,shape_up, nband=1)
       
    # Downscale the fitted coefficients
    CoefA_Downscale = reproject_dataset_example(
                                  CoefA_fileName_Optie2, ndvi_fileName)
    CoefA = CoefA_Downscale.GetRasterBand(1).ReadAsArray()
       
    CoefB_Downscale = reproject_dataset_example(
                                  CoefB_fileName_Optie2, ndvi_fileName)
    CoefB = CoefB_Downscale.GetRasterBand(1).ReadAsArray()
        
    CoefC_downscale = reproject_dataset_example(
                                CoefC_fileName_Optie2, ndvi_fileName)
    CoefC = CoefC_downscale.GetRasterBand(1).ReadAsArray()

    # Calculate the surface temperature based on the fitted coefficents and NDVI
    temp_surface_sharpened=CoefA*NDVI**2+CoefB*NDVI+CoefC
    temp_surface_sharpened[temp_surface_sharpened < 250] = np.nan					
    temp_surface_sharpened[temp_surface_sharpened > 400] = np.nan	
				
    save_GeoTiff_proy(dest_down,temp_surface_sharpened, temp_surface_sharpened_fileName,shape_down, nband=1)	
    return(temp_surface_sharpened)		
    
#------------------------------------------------------------------------------  
def CalculateSnowWaterMask(NDVI,shape_lsc,water_mask_temp,Surface_temp):
   '''
   Devides the temporaly water mask into a snow and water mask by using the surface temperature
   '''
   NDVI_nan=np.copy(NDVI)
   NDVI_nan[NDVI==0]=np.nan
   NDVI_nan=np.float32(NDVI_nan)
   NDVI_std=np.nanstd(NDVI_nan)
   NDVI_max=np.nanmax(NDVI_nan)
   NDVI_treshold_cold_pixels=NDVI_max-0.1*NDVI_std
   print 'NDVI treshold for cold pixels = ', '%0.3f' % NDVI_treshold_cold_pixels
   ts_moist_veg_min=np.nanmin(Surface_temp[NDVI>NDVI_treshold_cold_pixels])
   
   # calculate new water mask 
   mask=np.zeros((shape_lsc[1], shape_lsc[0]))
   mask[np.logical_and(np.logical_and(water_mask_temp==1, Surface_temp <= 275),NDVI>=0.3)]=1   
   snow_mask=np.copy(mask)
   
   # calculate new water mask 
   mask=np.zeros((shape_lsc[1], shape_lsc[0]))
   mask[np.logical_and(water_mask_temp==1, Surface_temp > 273)]=1   
   water_mask=np.copy(mask)
   
   return(snow_mask,water_mask,ts_moist_veg_min, NDVI_max, NDVI_std)
     
#------------------------------------------------------------------------------
def Water_Mask(shape_lsc,Reflect):
    """
    Calculates the water mask
    """
    mask = np.zeros((shape_lsc[1], shape_lsc[0]))
    mask[np.logical_and(Reflect[:, :, 3] < Reflect[:, :, 2],
                        Reflect[:, :, 4] < Reflect[:, :, 1])] = 1.0
    water_mask_temp = np.copy(mask)
    return(water_mask_temp)
				
#------------------------------------------------------------------------------				
if __name__ == '__main__':
    main()