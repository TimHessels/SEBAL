# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:09:38 2017

@author: tih
"""

import os
import gdal
import osr
import scipy
import numpy as np
import pandas as pd

Startdate ="2017-01-01"
Enddate ="2017-21-21"
Temp_folder = r"K:\Weather_Data\Model\GLDAS\three_hourly\tair_f_inst\Tair_GLDAS-NOAH_C_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
Pres_folder = r"K:\Weather_Data\Model\GLDAS\three_hourly\psurf_f_inst\P_GLDAS-NOAH_kpa_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
Hum_folder = r"K:\Weather_Data\Model\GLDAS\three_hourly\qair_f_inst\Hum_GLDAS-NOAH_kg-kg_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"
out_folder = r"K:\Weather_Data\Model\GLDAS\three_hourly\relative_humidity_inst\Humidity_GLDAS-NOAH_Percentage_3hour_{yyyy}.{mm:02d}.{dd:02d}_4.tif"

folder_dir_out = os.path.dirname(out_folder)

if not os.path.exists(folder_dir_out):
    os.makedirs(folder_dir_out)

Dates = pd.date_range(Startdate, Enddate, freq = "D")




def Open_array_info(filename=''):
    """
    Opening a tiff info, for example size of array, projection and transform matrix.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file

    """
    try:
        if filename.split('.')[-1] == 'tif':
            f = gdal.Open(r"%s" %filename)
        else:
            f = filename
    except:
            f = filename       
    try:
        geo_out = f.GetGeoTransform()
        proj = f.GetProjection()
        size_X = f.RasterXSize
        size_Y = f.RasterYSize
        f = None
    except:
        print('%s does not exists' %filename)
        
    return(geo_out, proj, size_X, size_Y)
           
    
def Save_as_tiff(name='', data='', geo='', projection=''):
    """
    This function save the array as a geotiff

    Keyword arguments:
    name -- string, directory name
    data -- [array], dataset of the geotiff
    geo -- [minimum lon, pixelsize, rotation, maximum lat, rotation,
            pixelsize], (geospatial dataset)
    projection -- integer, the EPSG code
    """
    
    dir_name = os.path.dirname(name)
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        
    # save as a geotiff
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name, int(data.shape[1]), int(data.shape[0]), 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    srse = osr.SpatialReference()
    if projection == '':
        srse.SetWellKnownGeogCS("WGS84")

    else:
        try:
            if not srse.SetWellKnownGeogCS(projection) == 6:
                srse.SetWellKnownGeogCS(projection)
            else:
                try:
                    srse.ImportFromEPSG(int(projection))
                except:
                    srse.ImportFromWkt(projection)
        except:
            try:
                srse.ImportFromEPSG(int(projection))
            except:
                srse.ImportFromWkt(projection)

    dst_ds.SetProjection(srse.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds = None
    return()
    
def gap_filling(dataset, NoDataValue, method = 1):
    """
    This function fills the no data gaps in a numpy array

    Keyword arguments:
    dataset -- 'C:/'  path to the source data (dataset that must be filled)
    NoDataValue -- Value that must be filled
    """
    import watertools.General.data_conversions as DC

    try:
        if dataset.split('.')[-1] == 'tif':
            # Open the numpy array
            data = Open_tiff_array(dataset)
            Save_as_tiff = 1
        else:
            data = dataset
            Save_as_tiff = 0
    except:
        data = dataset
        Save_as_tiff = 0

    # fill the no data values
    if NoDataValue is np.nan:
        mask = ~(np.isnan(data))
    else:
        mask = ~(data==NoDataValue)
    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    xym = np.vstack( (np.ravel(xx[mask]), np.ravel(yy[mask])) ).T
    data0 = np.ravel( data[:,:][mask] )

    if method == 1:
        interp0 = scipy.interpolate.NearestNDInterpolator( xym, data0 )
        data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

    if method == 2:
        interp0 = scipy.interpolate.LinearNDInterpolator( xym, data0 )
        data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

    if Save_as_tiff == 1:
        EndProduct=dataset[:-4] + '_GF.tif'

        # collect the geoinformation
        geo_out, proj, size_X, size_Y = Open_array_info(dataset)

        # Save the filled array as geotiff
        DC.Save_as_tiff(name=EndProduct, data=data_end, geo=geo_out, projection=proj)

    else:
        EndProduct = data_end

    return (EndProduct) 
   
def Open_tiff_array(filename='', band=''):
    """
    Opening a tiff array.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    band -- integer
        Defines the band of the tiff that must be opened.
    """
    f = gdal.Open(filename)
    if f is None:
        print('%s does not exists' %filename)
    else:
        if band == '':
            band = 1
        Data = f.GetRasterBand(band).ReadAsArray()
    return(Data)

for Date in Dates:
    
    Day = Date.day
    Month = Date.month
    Year = Date.year
    
    Tempfile_one = Temp_folder.format(yyyy = Year, mm = Month, dd = Day)
    Presfile_one = Pres_folder.format(yyyy = Year, mm = Month, dd = Day)
    Humfile_one = Hum_folder.format(yyyy = Year, mm = Month, dd = Day)
    out_folder_one = out_folder.format(yyyy = Year, mm = Month, dd = Day)

    geo_out, proj, size_X, size_Y = Open_array_info(Tempfile_one)
    Tdata = Open_tiff_array(Tempfile_one)
    Tdata[Tdata<-900]=-9999
    Pdata = Open_tiff_array(Presfile_one)
    Hdata = Open_tiff_array(Humfile_one)
    Pdata[Pdata<0]=-9999
    Hdata[Hdata<0]=-9999
    
    # gapfilling
    Tdata = gap_filling(Tdata,-9999)
    Pdata = gap_filling(Pdata,-9999)
    Hdata = gap_filling(Hdata,-9999)

    
    Esdata = 0.6108*np.exp((17.27*Tdata)/(Tdata+237.3))
    HumData = np.minimum((1.6077717*Hdata*Pdata/Esdata),1)*100
    HumData = HumData.clip(0,100)
                    
    Save_as_tiff(out_folder_one,HumData,geo_out,"WGS84")     
