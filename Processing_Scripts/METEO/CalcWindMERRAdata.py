# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 06:33:10 2022

@author: IrriWatch
"""


import os
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC
import numpy as np
import pandas as pd

Startdate ="2020-09-29"
Enddate ="2021-10-31"
# Temp_folder = r"K:\Project_Jain\Weather_Data\Model\GLDAS\three_hourly\tair_f_inst\Tair_GLDAS-NOAH_C_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
# Pres_folder = r"K:\Project_Jain\Weather_Data\Model\GLDAS\three_hourly\psurf_f_inst\P_GLDAS-NOAH_kpa_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
# Hum_folder = r"K:\Project_Jain\Weather_Data\Model\GLDAS\three_hourly\qair_f_inst\Hum_GLDAS-NOAH_kg-kg_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
# out_folder = r"K:\Project_Jain\Weather_Data\Model\GLDAS\three_hourly\relative_humidity_inst\Humidity_GLDAS-NOAH_Percentage_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"

wind_y_format = r"D:\Project_WAKAS\Input_Data\MERRA\Northward_Wind\hourly_MERRA2\v2m_MERRA_m-s-1_hourly_{yyyy}.{mm:02d}.{dd:02d}_H05.M00.tif"
wind_x_format = r"D:\Project_WAKAS\Input_Data\MERRA\Eastward_Wind\hourly_MERRA2\u2m_MERRA_m-s-1_hourly_{yyyy}.{mm:02d}.{dd:02d}_H05.M00.tif"
output_format = r"D:\Project_WAKAS\Input_Data\MERRA\Wind\hourly_MERRA2\Wind_MERRA_m-s-1_hourly_{yyyy}.{mm:02d}.{dd:02d}_H05.M00.tif"






def Calc_Humidity(Temp_format, P_format, Hum_format, output_format, Startdate, Enddate, freq ="D"):

    folder_dir_out = os.path.dirname(output_format)
    
    if not os.path.exists(folder_dir_out):
        os.makedirs(folder_dir_out)

    Dates = pd.date_range(Startdate, Enddate, freq = "D")
    
    for Date in Dates:
        
        print(Date)
        
        Day = Date.day
        Month = Date.month
        Year = Date.year
        
        Windyfile_one = wind_y_format.format(yyyy = Year, mm = Month, dd = Day)
        Windxfile_one = wind_x_format.format(yyyy = Year, mm = Month, dd = Day)
        out_folder_one = output_format.format(yyyy = Year, mm = Month, dd = Day)
    
        geo_out, proj, size_X, size_Y = RC.Open_array_info(Windxfile_one)
        Windxdata = RC.Open_tiff_array(Windxfile_one)
        Windxdata[Windxdata<-900]=-9999

        Windydata = RC.Open_tiff_array(Windyfile_one)
        Windydata[Windxdata<-900]=-9999
        
        # gapfilling
        Windydata = RC.gap_filling(Windydata,-9999)
        Windxdata = RC.gap_filling(Windxdata,-9999)

        Wind = np.sqrt(Windxdata**2 + Windydata **2)
                        
        DC.Save_as_tiff(out_folder_one,Wind,geo_out,"WGS84") 

    return()               