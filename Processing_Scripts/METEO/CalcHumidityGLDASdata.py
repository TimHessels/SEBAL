# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:09:38 2017

@author: tih
"""

Tfile = r"J:\Tyler\Input\Meteo\daily\avgsurft_inst\mean\T_GLDAS-NOAH_C_daily_2016.06.15.tif"
Pfile = r"J:\Tyler\Input\Meteo\daily\psurf_f_inst\mean\P_GLDAS-NOAH_kpa_daily_2016.06.15.tif"
Hfile = r"J:\Tyler\Input\Meteo\daily\qair_f_inst\mean\Hum_GLDAS-NOAH_kg-kg_daily_2016.06.15.tif"
Outfilename = r"J:\Tyler\Input\Meteo\daily\Hum_Calculated\Humidity_percentage_Calculated_daily.tif"

import gdal
import os
import wa.General.raster_conversions as RC
import wa.General.data_conversions as DC
import numpy as np


geo_out, proj, size_X, size_Y = RC.Open_array_info(Tfile)
Tdata = RC.Open_tiff_array(Tfile)
Tdata[Tdata<-900]=np.nan
Pdata = RC.Open_tiff_array(Pfile)
Hdata = RC.Open_tiff_array(Hfile)


Esdata = 0.6108*np.exp((17.27*Tdata)/(Tdata+237.3))
HumData = np.minimum((1.6077717*Hdata*Pdata/Esdata),1)*100
                
DC.Save_as_tiff(Outfilename,HumData,geo_out,"WGS84")                