#!/usr/bin/env python

import os
import numpy as np
from scipy import interpolate
from osgeo import gdal, osr

#functions

def GetGeoInfo(filename, Subdataset = 0):
    """Gives geo-information derived from a georeferenced map
    
    filename: file to be scrutinized
    subdataset: layer to be used in case of HDF4 format
    """
    SourceDS = gdal.Open(filename, gdal.GA_ReadOnly)
    Type = SourceDS.GetDriver().ShortName
    if Type == 'HDF4' or Type == 'netCDF':
        SourceDS = gdal.Open(SourceDS.GetSubDatasets()[Subdataset][0])
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    #DataType = SourceDS.GetRasterBand(1).DataType
    #DataType = gdal.GetDataTypeName(DataType)
    driver = gdal.GetDriverByName(Type)
    return driver, NDV, xsize, ysize, GeoT, Projection#, DataType

    
def OpenAsArray(filename, Bandnumber = 1, dtype = 'float32', nan_values = False):
    """Function opens a geo-map as an numpy array. 
    
    filename: map to be opened
    Bandnumber: band to open as array (or layer in case of HDF4)
    dtype: datatype of output array
    nan_values: if True, the NDV values in the output array are substituted 
    with Numpy NaN values. NOTE: dtype has to be float in that case.
    """
    datatypes = {"uint8": np.uint8, "int8": np.int8, "uint16": np.uint16, "int16":  np.int16, "Int16":  np.int16, "uint32": np.uint32,
    "int32": np.int32, "float32": np.float32, "float64": np.float64, "complex64": np.complex64, "complex128": np.complex128,
    "Int32": np.int32, "Float32": np.float32, "Float64": np.float64, "Complex64": np.complex64, "Complex128": np.complex128,}
    DataSet = gdal.Open(filename, gdal.GA_ReadOnly)
    Type = DataSet.GetDriver().ShortName
    if Type == 'GTiff':
        Subdataset = DataSet.GetRasterBand(Bandnumber)
        NDV = Subdataset.GetNoDataValue()
    if Type == 'HDF4':
        Subdataset = gdal.Open(DataSet.GetSubDatasets()[Bandnumber][0])
        NDV = int(Subdataset.GetMetadata()['_FillValue'])
    Array = Subdataset.ReadAsArray().astype(datatypes[dtype])
    if nan_values:
        Array[Array == NDV] = np.nan
    #missing value filling
    x = np.arange(0, Array.shape[1])
    y = np.arange(0, Array.shape[0])
    #mask invalid values
    Array = np.ma.masked_invalid(Array)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~Array.mask]
    y1 = yy[~Array.mask]
    newarr = Array[~Array.mask]
    Array = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='nearest')
    return Array

def CreateGeoTiff(Name, Array, driver, NDV, xsize, ysize, GeoT, Projection,output_folder,subdir):
    """Function creates a geotiff file from a Numpy Array.
    Name: output file
    Array: numpy array to save
    """
    datatypes = {"uint8": 1, "int8": 1, "uint16": 2, "int16": 3, "Int16": 3, "uint32": 4,
    "int32": 5, "float32": 6, "float64": 7, "complex64": 10, "complex128": 11,
    "Int32": 5, "Float32": 6, "Float64": 7, "Complex64": 10, "Complex128": 11,}
    
    file_name = output_folder + subdir + Name + '.tif'
    if not os.path.isdir(output_folder + subdir):
        os.makedirs(output_folder + subdir)
        
    DataSet = driver.Create(file_name,xsize,ysize,1,datatypes[Array.dtype.name])
    if NDV is None:
        NDV = -9999
    Array[np.isnan(Array)] = NDV
    DataSet.GetRasterBand(1).SetNoDataValue(NDV)
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projection.ExportToWkt())
    DataSet.GetRasterBand(1).WriteArray(Array)
    DataSet = None