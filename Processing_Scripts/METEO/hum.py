# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:23:59 2017

@author: joo
"""
import os
import numpy as np
import wa.General.raster_conversions as RC
import wa.General.data_conversions as DC

tmean = r"D:\Weather\Lebanon\Bekaa-ETref\2017\Weather_Data\Model\GLDAS\daily\tair_f_inst\mean/"
tmeans = [os.path.join(tmean,fn) for fn in next(os.walk(tmean))[2] if fn[-4:] == '.tif']

humid = r"D:\Weather\Lebanon\Bekaa-ETref\2017\Weather_Data\Model\GLDAS\daily\qair_f_inst\mean/"
humids = [os.path.join(humid,fn) for fn in next(os.walk(humid))[2] if fn[-4:] == '.tif']

press = r"D:\Weather\Lebanon\Bekaa-ETref\2017\Weather_Data\Model\GLDAS\daily\psurf_f_inst\mean/"
presss = [os.path.join(press,fn) for fn in next(os.walk(press))[2] if fn[-4:] == '.tif']

Outfilepath = r"D:\Weather\Lebanon\Bekaa-ETref\2017\Weather_Data\Model\GLDAS\daily\humid/"

for i in range (0,len(tmeans)):
    geo_out, proj, size_X, size_Y = RC.Open_array_info(tmeans[i])
    Tdata = RC.Open_tiff_array(tmeans[i])
    Tdata[Tdata<-900]=np.nan
    Pdata = RC.Open_tiff_array(presss[i])
    Hdata = RC.Open_tiff_array(humids[i])
    
    Esdata = 0.6108*np.exp((17.27*Tdata)/(Tdata+237.3))
    HumData = np.minimum((1.6077717*Hdata*Pdata/Esdata),1)*100

    datestamp = tmeans[i][-14:-4]    
    
    Outfilename = os.path.join(Outfilepath,datestamp + '.tif')    
    
    DC.Save_as_tiff(Outfilename,HumData,geo_out,"WGS84")  
