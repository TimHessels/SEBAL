# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:09:38 2017

@author: tih
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

Temp_format = r"D:\Project_WAKAS\Input_Data\MERRA\Air_Temperature\daily_MERRA2\t2m_MERRA_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif"
P_format = r"D:\Project_WAKAS\Input_Data\MERRA\Surface_Pressure\daily_MERRA2\ps_MERRA_kpa_daily_{yyyy}.{mm:02d}.{dd:02d}.tif"
Hum_format = r"D:\Project_WAKAS\Input_Data\MERRA\Specific_Humidity\daily_MERRA2\q2m_MERRA_kg-kg-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif"
output_format = r"D:\Project_WAKAS\Input_Data\MERRA\relative_humidity\Humidity_MERRA_kg-kg-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif"






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
        
        Tempfile_one = Temp_format.format(yyyy = Year, mm = Month, dd = Day)
        Presfile_one = P_format.format(yyyy = Year, mm = Month, dd = Day)
        Humfile_one = Hum_format.format(yyyy = Year, mm = Month, dd = Day)
        out_folder_one = output_format.format(yyyy = Year, mm = Month, dd = Day)
    
        geo_out, proj, size_X, size_Y = RC.Open_array_info(Tempfile_one)
        Tdata = RC.Open_tiff_array(Tempfile_one)
        if "MERRA_K" in Temp_format:
            Tdata = Tdata - 273.15
        
        Tdata[Tdata<-900]=-9999
        Pdata = RC.Open_tiff_array(Presfile_one)
        Hdata = RC.Open_tiff_array(Humfile_one)
        Pdata[Pdata<0]=-9999
        Hdata[Hdata<0]=-9999
        
        # gapfilling
        Tdata = RC.gap_filling(Tdata,-9999)
        Pdata = RC.gap_filling(Pdata,-9999)
        Hdata = RC.gap_filling(Hdata,-9999)
    
        Esdata = 0.6108*np.exp((17.27*Tdata)/(Tdata+237.3))
        HumData = np.minimum((1.6077717*Hdata*Pdata/Esdata),1)*100
        HumData = HumData.clip(0,100)
                        
        DC.Save_as_tiff(out_folder_one,HumData,geo_out,"WGS84") 

    return()               