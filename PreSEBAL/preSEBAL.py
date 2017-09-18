# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 15:09:49 2016

@author: tih
"""
import numpy as np
import os
import gdal
from openpyxl import load_workbook
import osr
from datetime import datetime, timedelta
import pandas as pd

import SEBAL

def main():
    ############################## INPUT ##########################################		
        # Input for preSEBAL.py
    for number in range(2,3):                                                                                   # Number defines the column of the inputExcel
        inputExcel=r'G:\SEBAL_Tadla\Excel\InputEXCEL_v3_3_6.xlsx'               # The excel with all the SEBAL input data
        VegetationExcel = r'G:\SEBAL_Tadla\Excel\Vegetation height model.xlsx'  # This excel defines the p and c factor and vegetation height.
        output_folder = r'G:\SEBAL_Tadla\Preprocessing_Output'         # Output folder
        LU_data_FileName = r'G:\SEBAL_Tadla\LandCover\LU_map.tif'        # Path to Land Use map
    
        # optional paramater 
        DSSF_Folder=r'G:\SEBAL_Tadla\LANDSAF'
    
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
        lat,lon,lat_fileName,lon_fileName=SEBAL.DEM_lat_lon(DEM_fileName,output_folder_temp)

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
            year, DOY, hour, minutes, UTM_Zone, Sun_elevation = SEBAL.info_general_metadata(Landsat_meta_fileName) # call definition info_general_metadata
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
        dest, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset(
                       DEM_fileName, pixel_spacing, UTM_Zone=UTM_Zone)
        band = dest.GetRasterBand(1)   # Get the reprojected dem band
        ncol = dest.RasterXSize        # Get the reprojected dem column size
        nrow = dest.RasterYSize        # Get the reprojected dem row size
        shape=[ncol, nrow]
       
        # Read out the DEM band and print the DEM properties
        data_DEM = band.ReadAsArray(0, 0, ncol, nrow)

        # 2) Latitude file - reprojection    
        # reproject latitude to the landsat projection and save as tiff file																																
        lat_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset(
                        lat_fileName, pixel_spacing, UTM_Zone=UTM_Zone)
                    
        # Get the reprojected latitude data															
        lat_proy = lat_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)
     
        # 3) Longitude file - reprojection
        # reproject longitude to the landsat projection	 and save as tiff file	
        lon_rep, ulx_dem, lry_dem, lrx_dem, uly_dem, epsg_to = SEBAL.reproject_dataset(lon_fileName, pixel_spacing, UTM_Zone=UTM_Zone)

        # Get the reprojected longitude data	
        lon_proy = lon_rep.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow)
        
        lon_fileName = os.path.join(output_folder_temp,'lon_resh.tif')
        SEBAL.save_GeoTiff_proy(dest, lon_proy, lon_fileName, shape, nband=1)	
				
        # Calculate slope and aspect from the reprojected DEM
        deg2rad,rad2deg,slope,aspect=SEBAL.Calc_Gradient(data_DEM, pixel_spacing)

        # calculate the coz zenith angle
        Ra_mountain_24, Ra_inst, cos_zn_resh, dr, phi, delta = SEBAL.Calc_Ra_Mountain(lon,DOY,hour,minutes,lon_proy,lat_proy,slope,aspect) 
        cos_zn_fileName = os.path.join(output_folder_temp,'cos_zn.tif')
        SEBAL.save_GeoTiff_proy(dest, cos_zn_resh, cos_zn_fileName, shape, nband=1)
        
        # Save the Ra
        Ra_inst_fileName = os.path.join(output_folder_temp,'Ra_inst.tif')
        SEBAL.save_GeoTiff_proy(dest, Ra_inst, Ra_inst_fileName, shape, nband=1)	
        Ra_mountain_24_fileName = os.path.join(output_folder_temp,'Ra_mountain_24.tif')
        SEBAL.save_GeoTiff_proy(dest, Ra_mountain_24, Ra_mountain_24_fileName, shape, nband=1)	
        
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
                Ra_inst_3Km_dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(Ra_inst_fileName, file_Landsaf_inst, method = 2)   
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
                SEBAL.save_GeoTiff_proy(Ra_inst_3Km_dest, Transmissivity_3Km, Transmissivity_3Km_fileName, shape_trans, nband=1)	
                
                # Reproject Transmissivity to match DEM (now this is done by using the nearest neighbour method)
                Transmissivity_inst_dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(Transmissivity_3Km_fileName, cos_zn_fileName)  
                Transmissivity_inst = Transmissivity_inst_dest.GetRasterBand(1).ReadAsArray()
                Transmissivity_inst_fileName = os.path.join(TRANS_outfolder,'Transmissivity_inst_%s.tif' %Var_name)                
                SEBAL.save_GeoTiff_proy(Transmissivity_inst_dest, Transmissivity_inst, Transmissivity_inst_fileName, shape, nband=1)	
                
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
                Ra_24_3Km_dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(Ra_mountain_24_fileName, file_Landsaf_inst, method = 2)   
                Ra_24_3Km = Ra_24_3Km_dest.GetRasterBand(1).ReadAsArray()
                Ra_24_3Km[Ra_24_3Km==0] = np.nan  
                 
                # Get shape LANDSAF data
                shape_trans=[dest_Rs_inst_3Km.RasterXSize , dest_Rs_inst_3Km.RasterYSize ]
                
                # Calculate Transmissivity 3Km
                Transmissivity_24_3Km = Rs_24_3Km/Ra_24_3Km
            
                Transmissivity_24_3Km_fileName = os.path.join(output_folder_temp,'Transmissivity_24_3Km.tif')
                SEBAL.save_GeoTiff_proy(Ra_24_3Km_dest, Transmissivity_24_3Km, Transmissivity_24_3Km_fileName, shape_trans, nband=1)	
                
                # Reproject Transmissivity to match DEM (now this is done by using the nearest neighbour method)
                Transmissivity_24_dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(Transmissivity_24_3Km_fileName, cos_zn_fileName)  
                Transmissivity_24 = Transmissivity_24_dest.GetRasterBand(1).ReadAsArray()
                Transmissivity_24_fileName = os.path.join(TRANS_outfolder,'Transmissivity_24_%s.tif' %Var_name)                
                SEBAL.save_GeoTiff_proy(Transmissivity_24_dest, Transmissivity_24, Transmissivity_24_fileName, shape, nband=1)	

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
            Lmin, Lmax, k1_c, k2_c = SEBAL.info_band_metadata(Landsat_meta_fileName, Bands)
           
            # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
            # for the different Landsat images (L5, L7, or L8)
            ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
            ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
            ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])
    
            # Open one band - To get the metadata of the landsat images only once (to get the extend)
            src_FileName = os.path.join(input_folder, '%s_B2.TIF' % Name_Landsat_Image)  # before 10!
            ls,band_data,ulx,uly,lrx,lry,x_size_ls,y_size_ls = SEBAL.Get_Extend_Landsat(src_FileName)
         
            # Crop the Landsat images to the DEM extent -
            dst_FileName = os.path.join(output_folder_temp,'cropped_LS_b2.tif')  # Before 10 !!
     							
            # Clip the landsat image to match the DEM map											
            lsc, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(src_FileName, cos_zn_fileName)	
            data_LS = lsc.GetRasterBand(1).ReadAsArray()
            SEBAL.save_GeoTiff_proy(dest, data_LS, dst_FileName, shape, nband=1)	
    
            # Get the extend of the remaining landsat file after clipping based on the DEM file	
            lsc,band_data,ulx,uly,lrx,lry,x_size_lsc,y_size_lsc = SEBAL.Get_Extend_Landsat(dst_FileName)

            # Create the corrected signals of Landsat in 1 array
            Reflect = SEBAL.Landsat_Reflect(Bands,input_folder,Name_Landsat_Image,output_folder,shape,Lmax,Lmin,ESUN_L5,ESUN_L7,ESUN_L8,cos_zn_resh,dr,Landsat_nr, cos_zn_fileName)

            # Calculate temporal water mask
            water_mask_temp=SEBAL.Water_Mask(shape,Reflect)       

            # Calculate SAVI
            L_SAVI = 0.5
            SAVI = SEBAL.Calc_SAVI(Reflect, L_SAVI)			

            # Calculate NDVI
            NDVI = SEBAL.Calc_NDVI(Reflect)

            # Calculate albedo
            albedo = SEBAL.Calc_albedo(Reflect)

            NDVI_FileName = os.path.join(NDVI_outfolder,'NDVI_LS_%s.tif'%Var_name)
            SEBAL.save_GeoTiff_proy(dest, NDVI, NDVI_FileName, shape, nband=1)
								
            SAVI_FileName = os.path.join(SAVI_outfolder,'SAVI_LS_%s.tif'%Var_name)
            SEBAL.save_GeoTiff_proy(dest, SAVI, SAVI_FileName, shape, nband=1)
								
            albedo_FileName = os.path.join(Albedo_outfolder,'Albedo_LS_%s.tif'%Var_name)
            SEBAL.save_GeoTiff_proy(dest, albedo, albedo_FileName, shape, nband=1)

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
            Temp_inst = SEBAL.Reshape_Reproject_Input_data(Temp_inst_name, Temp_inst_fileName, cos_zn_fileName)

        try:
            RH_inst = float(ws['D%d' %number].value)                # Instantaneous Relative humidity (%)
 
        # if the data is not a value, than open as a string							
        except:
            RH_inst_name = '%s' %str(ws['D%d' %number].value) 
            RH_inst_fileName = os.path.join(output_folder, 'Temp', 'RH_inst_input.tif')
            RH_inst = SEBAL.Reshape_Reproject_Input_data(RH_inst_name, RH_inst_fileName, cos_zn_fileName)			
         
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
                PROBAV, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(
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
								
            SEBAL.save_GeoTiff_proy(dest, Surface_Albedo_PROBAV, Albedo_FileName, shape, nband=1)	  

            # Calculate the NDVI based on PROBA-V     
            n218_memory = spectral_reflectance_PROBAV[:, :, 2] + spectral_reflectance_PROBAV[:, :, 3]
            NDVI = np.zeros((shape[1], shape[0]))
            NDVI[n218_memory != 0] =  ( spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] - spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] )/ ( spectral_reflectance_PROBAV[:, :, 2][n218_memory != 0] + spectral_reflectance_PROBAV[:, :, 3][n218_memory != 0] )

            NDVI_FileName = os.path.join(NDVI_outfolder,'NDVI_PROBAV_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, NDVI, NDVI_FileName, shape, nband=1)	

            SAVI_FileName = os.path.join(SAVI_outfolder,'SAVI_PROBAV_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, SAVI, SAVI_FileName, shape, nband=1)															
				

     #################### Calculate vegetation parameters ##########################################	

        FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity = SEBAL.Calc_vegt_para(NDVI,SAVI,water_mask_temp,shape)

        if Image_Type == 1:

            LAI_FileName = os.path.join(LAI_outfolder,'LAI_LS_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, LAI, LAI_FileName, shape, nband=1)	

        if Image_Type == 2:

            LAI_FileName = os.path.join(LAI_outfolder,'LAI_PROBAV_%s.tif' %Var_name) 				
            SEBAL.save_GeoTiff_proy(dest, LAI, LAI_FileName, shape, nband=1)	
            
     #################### Calculate thermal for Landsat ##########################################	

        if Image_Type == 1:
								
            therm_data = SEBAL.Landsat_therm_data(Bands,input_folder,Name_Landsat_Image,output_folder,ulx_dem,lry_dem,lrx_dem,uly_dem,shape)          
            Surface_temp=SEBAL.Calc_surface_water_temp(Temp_inst,Landsat_nr,Lmax,Lmin,therm_data,b10_emissivity,k1_c,k2_c,eact_inst,shape,water_mask_temp,Bands_thermal,Rp,tau_sky,surf_temp_offset,Image_Type)
            therm_data_FileName = os.path.join(Surface_Temperature_outfolder,'Surface_Temperature_LS_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, Surface_temp, therm_data_FileName, shape, nband=1)	

     ################################ Apply thermal sharpening #######################################

            # Thermal Sharpening of Thermal data Landsat

            # Upscale DEM to 90m
            pixel_spacing_upscale=90

            dest_90, ulx_dem_90, lry_dem_90, lrx_dem_90, uly_dem_90, epsg_to = SEBAL.reproject_dataset(
                                                      DEM_fileName, pixel_spacing_upscale, UTM_Zone = UTM_Zone)
                                                      
            DEM_90 = dest_90.GetRasterBand(1).ReadAsArray()
            Y_raster_size_90 = dest_90.RasterYSize				
            X_raster_size_90 = dest_90.RasterXSize
            shape_90=([X_raster_size_90, Y_raster_size_90])
						
            proyDEM_fileName_90 = os.path.join(output_folder_temp, 'DEM_90m.tif')      
      
            SEBAL.save_GeoTiff_proy(dest_90, DEM_90, proyDEM_fileName_90, shape_90, nband=1)
	  	
            # save landsat surface temperature
            surf_temp_fileName = os.path.join(output_folder_temp, '%s_%s_surface_temp_%s_%s_%s.tif' %(sensor1, sensor2, res2, year, DOY))
            SEBAL.save_GeoTiff_proy(lsc, Surface_temp, surf_temp_fileName, shape, nband=1)							

            # Upscale NDVI data																																	
            dest_up, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(
                                       NDVI_FileName, proyDEM_fileName_90)
  
            NDVI_Landsat_up = dest_up.GetRasterBand(1).ReadAsArray()

            # Upscale Thermal data
            dest_up, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(
                                       surf_temp_fileName, proyDEM_fileName_90)
            surface_temp_up = dest_up.GetRasterBand(1).ReadAsArray()
 
            # Define the width of the moving window box
            Box=7	
  
            temp_surface_sharpened_fileName = os.path.join(output_folder_temp,'LS_surface_temp_Sharpened_%s.tif' %Var_name)
  
            # Apply thermal sharpening           
            temp_surface_sharpened = SEBAL.Thermal_Sharpening(surface_temp_up, NDVI_Landsat_up, NDVI, Box, dest_up, output_folder, NDVI_FileName, shape, lsc, temp_surface_sharpened_fileName)	
	  
  
            # Divide temporal watermask in snow and water mask by using surface temperature
            snow_mask, water_mask, ts_moist_veg_min, NDVI_max, NDVI_std = SEBAL.CalculateSnowWaterMask(NDVI,shape,water_mask_temp,temp_surface_sharpened)
       
            # Replace water values by old values before shapening            
            temp_surface_sharpened[water_mask == 1] = Surface_temp[water_mask == 1]
            temp_surface_sharpened = np.where(np.isnan(temp_surface_sharpened),Surface_temp,temp_surface_sharpened)
            
		      # save landsat surface temperature
            surf_temp_fileName = os.path.join(Surface_Temperature_Sharp_outfolder ,'LS_surface_temp_sharpened_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(lsc, temp_surface_sharpened, surf_temp_fileName, shape, nband=1)	


    ################################## Calculate VIIRS surface temperature ########################

        if Image_Type == 2:

            # Define the VIIRS thermal data name
            VIIRS_data_name=os.path.join(input_folder, '%s' % (Name_VIIRS_Image_TB))
							
            # Reproject VIIRS thermal data								
            VIIRS, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(VIIRS_data_name, LAI_FileName)
																
            # Open VIIRS thermal data																		
            data_VIIRS = VIIRS.GetRasterBand(1).ReadAsArray()    
				
            # Define the thermal VIIRS output name
            proyVIIRS_fileName = os.path.join(output_folder_temp, 'Surface_Temp_VIIRS_%s.tif' %Var_name)
	 											
            # Save the thermal VIIRS data 												
            SEBAL.save_GeoTiff_proy(dest, data_VIIRS, proyVIIRS_fileName, shape, nband=1)	

            # Set the conditions for the brightness temperature (100m)
            brightness_temp=np.where(data_VIIRS>=250, data_VIIRS, np.nan)
 
            # Constants
            k1=606.399172
            k2=1258.78
            L_lambda_b10_100=((2*6.63e-34*(3.0e8)**2)/((11.45e-6)**5*(np.exp((6.63e-34*3e8)/(1.38e-23*(11.45e-6)*brightness_temp))-1)))*1e-6
         
            # Get Temperature for 100 and 375m resolution
            Temp_TOA_100 = SEBAL.Get_Thermal(L_lambda_b10_100,Rp,Temp_inst,tau_sky,tir_emis,k1,k2) 
       
            # Conditions for surface temperature (100m)
            n120_surface_temp=Temp_TOA_100.clip(250, 450)

            # Save the surface temperature of the VIIRS in 100m resolution
            temp_surface_100_fileName_beforeTS = os.path.join(Surface_Temperature_outfolder,'Surface_Temperature_VIIRS_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, n120_surface_temp, temp_surface_100_fileName_beforeTS, shape, nband=1)     
        
        
            print '---------------------------------------------------------'
            print '-------------------- Downscale VIIRS --------------------'
            print '---------------------------------------------------------'

        ################################ Thermal Sharpening #####################################################

            # Upscale DEM to 90m
            pixel_spacing_upscale=400

            dest_400, ulx_dem_400, lry_dem_400, lrx_dem_400, uly_dem_400, epsg_to = SEBAL.reproject_dataset(
                                                      DEM_fileName, pixel_spacing_upscale, UTM_Zone = UTM_Zone)
                                                      
            DEM_400 = dest_400.GetRasterBand(1).ReadAsArray()
            Y_raster_size_400 = dest_400.RasterYSize				
            X_raster_size_400 = dest_400.RasterXSize
            shape_400=([X_raster_size_400, Y_raster_size_400])
						
            proyDEM_fileName_400 = os.path.join(output_folder_temp, 'DEM_400m.tif')      
      
            SEBAL.save_GeoTiff_proy(dest_400, DEM_400, proyDEM_fileName_400, shape_400, nband=1)
    																		
            # Upscale thermal band VIIRS from 100m to 400m
            VIIRS_Upscale, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(temp_surface_100_fileName_beforeTS, proyDEM_fileName_400)
            data_Temp_Surf_400 = VIIRS_Upscale.GetRasterBand(1).ReadAsArray()
 
            # Upscale PROBA-V NDVI from 100m to 400m       
            NDVI_PROBAV_Upscale, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(
                              NDVI_FileName, proyDEM_fileName_400)
            data_NDVI_400 = NDVI_PROBAV_Upscale.GetRasterBand(1).ReadAsArray()

            # Define the width of the moving window box
            Box=9
           
            # The surface temperature temporary file name
            temp_surface_sharpened_fileName = os.path.join(output_folder_temp,'VIIRS_surface_temp_Sharpened_%s.tif' %Var_name)
  
            # Apply the surface temperature sharpening        
            temp_surface_sharpened = SEBAL.Thermal_Sharpening(data_Temp_Surf_400, data_NDVI_400, NDVI, Box, NDVI_PROBAV_Upscale, output_folder, NDVI_FileName, shape, dest, temp_surface_sharpened_fileName)	
      
            # Divide temporal watermask in snow and water mask by using surface temperature
            snow_mask, water_mask, ts_moist_veg_min, NDVI_max, NDVI_std = SEBAL.CalculateSnowWaterMask(NDVI,shape,water_mask_temp,temp_surface_sharpened)
       
            # Replace water values by old values before shapening            
            temp_surface_sharpened[water_mask == 1] = n120_surface_temp[water_mask == 1]
            temp_surface_sharpened = np.where(np.isnan(temp_surface_sharpened),n120_surface_temp,temp_surface_sharpened)  
            
		   # save landsat surface temperature
            surf_temp_fileName = os.path.join(Surface_Temperature_Sharp_outfolder ,'VIIRS_surface_temp_sharpened_%s.tif' %Var_name)
            SEBAL.save_GeoTiff_proy(dest, temp_surface_sharpened, surf_temp_fileName, shape, nband=1)	        

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
    dest1, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(LAI_FileName, LU_data_FileName)	 ## input after HANTS
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
    dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(Veg_Height_proj_FileName, LAI_FileName)			

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

    dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(dst_FileName, LAI_FileName)	

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

    dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset_example(dst_FileName, LAI_FileName)	

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
if __name__ == '__main__':
    main()