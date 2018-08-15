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


    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # FUNCTIONS
    #-------------------------------------------------------------------------
