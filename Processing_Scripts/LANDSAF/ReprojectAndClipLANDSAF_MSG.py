# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:32:01 2017
Tim Hessels

This script reprojects LANDSAF data (whole MSG disk) and clips to user extend

import is downloaded from: https://landsaf.ipma.pt/ and use the MSG extend

The data name needs to look similar as: HDF5_LSASAF_MSG_DSSF_MSG-Disk_201601101215.bz2
or HDF5_LSASAF_MSG_DSSF_MSG-Disk_201601101215

"""

import os
import shutil
import glob
import gdal
import SEBAL.pySEBAL.pySEBAL_code as SEBAL



# User parameters
latlim = [31, 33]
lonlim = [34, 36] 
input_folder = r"K:\Project_Jain\Radiation2" 

# Get some general data
temp_folder = os.path.join(input_folder, 'Temporary')
if not os.path.exists(temp_folder):
    os.mkdir(temp_folder)
    
os.chdir(input_folder)
for file_name in [e for e in glob.glob('HDF5_*') if not e.endswith('.tif')]:

    print('Processing: ', file_name)
    extension = os.path.splitext(file_name)[1]
    file_name_only = os.path.splitext(file_name)[0]
    file_name_tiff =file_name_only + '.tif'
    out_name = os.path.join(input_folder, file_name_tiff)

    if not os.path.exists(out_name):


        # unzip bunzipfile if needed 
        if str(extension) == '.bz2':
            fullCmd = '7z e %s -o%s' %(file_name,temp_folder)
            SEBAL.Run_command_window(fullCmd)
            input_folder_landsaf = temp_folder
        else:
            input_folder_landsaf = input_folder        
      
        # Get environmental variable
        SEBAL_env_paths = os.environ["SEBAL"].split(';')
        GDAL_env_path = SEBAL_env_paths[0]
        GDAL_TRANSLATE = os.path.join(GDAL_env_path, 'gdal_translate.exe')
        GDALWARP = os.path.join(GDAL_env_path, 'gdalwarp.exe')
    
        # Set projection of GOES/LANDSAF data
        output_GOES_projected = os.path.join(temp_folder,'GOES_projected_LANDSAF.tif')    
        fullCmd = '"%s" -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5570248.832537 5570248.832537 5570248.832537 -5570248.832537 HDF5:"%s"://DSSF "%s"' %(GDAL_TRANSLATE, os.path.join(input_folder_landsaf, file_name_only), output_GOES_projected)
        SEBAL.Run_command_window(fullCmd)
        
        # Reproject GOES/LANDSAF data and clip data
        fullCmd = '"%s" -overwrite -s_srs "+proj=geos +lon_0=0 +h=35785831 +x_0=0 +y_0=0 +a=6378169 +b=6356583.8 +units=m +no_defs" -t_srs EPSG:4326 -te %d %d %d %d -of GTiff "%s" "%s"' %(GDALWARP, lonlim[0], latlim[0], lonlim[1], latlim[1], output_GOES_projected, out_name)
        SEBAL.Run_command_window(fullCmd)

        dest = gdal.Open(out_name)
        sizeX = dest.RasterXSize
        sizeY = dest.RasterYSize        
        Array = dest.GetRasterBand(1).ReadAsArray()
        Array = Array * 0.1
        Array[Array < 0] = -9999

        
        SEBAL.save_GeoTiff_proy(dest,Array,out_name,[Array.shape[1],Array.shape[0]],1)


# Save data as Tiff
if os.path.exists(temp_folder):
    shutil.rmtree(temp_folder)    









