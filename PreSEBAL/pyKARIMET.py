from __future__ import division
import os
import netCDF4
import numpy as np
import osr
import numpy
import scipy.interpolate
import gdal
import glob
from openpyxl import load_workbook

def main(VegetationExcel, nc_path_lst,nc_path_meteo, CT):
  
    # Open Excelfile
    wb_veg = load_workbook(VegetationExcel, data_only=True)
    
    # Open worksheets
    ws_veg = wb_veg['General_Input']
    
    # input parameters
    input_folder = r"%s" %str(ws_veg['B7'].value)

    # Create input datasets for HANTS
    output_folder = os.path.join(input_folder, "PreSEBAL_SEBAL_out")
    
    # Define End name
    name_end = 'KARIMET.nc'

    nc_file = netCDF4.Dataset(nc_path_lst, 'r')
    latx = nc_file.variables['latitude'][:]
    lonx = nc_file.variables['longitude'][:]

    X_tiles = int(np.ceil(len(lonx)/100))
    Y_tiles = int(np.ceil(len(latx)/100))

    for x_tile in range(0,X_tiles):
        x_start = x_tile * 100
        x_end = x_start + 100
        if x_end > len(lonx):
                x_end = len(lonx)

        for y_tile in range(0,Y_tiles):
            y_start = y_tile * 100
            y_end = y_start + 100
            if y_end > len(latx):
                y_end = len(latx)
            name='KARIMET_X%d_Y%d.nc' %(int(x_tile + 1), int(y_tile + 1))
            KARIMET(output_folder, nc_path_lst, nc_path_meteo, CT, x_start, x_end, y_start, y_end, name)

    Merge_NC_files(output_folder, nc_path_lst, nc_path_meteo, name_end)
    return(name_end)

def Merge_NC_files(output_folder, nc_path_lst, nc_path_meteo, name_end):

    nc_file = netCDF4.Dataset(nc_path_lst, 'r')
    latx = nc_file.variables['latitude'][:]
    lonx = nc_file.variables['longitude'][:]

    nc_file2 = netCDF4.Dataset(nc_path_meteo, 'r')
    time = nc_file2.variables['time'][:]

    ###################### Create netcdf ##########################################

    epsg=4326
    fill_val=-9999.0
    spa_ref = Spatial_Reference(epsg)
    Scaling_factor = 0.0001

    fill_val = -9999
    nc_fileo = netCDF4.Dataset(os.path.join(output_folder,name_end), 'w', format="NETCDF4")
    lat_n = len(latx)
    lon_n = len(lonx)
    nc_fileo.createDimension('latitude', lat_n)
    nc_fileo.createDimension('longitude', lon_n)
    nc_fileo.createDimension('time', None)

    crs_var = nc_fileo.createVariable('crs', 'i4')
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.crs_wkt = spa_ref

    timevar = nc_fileo.createVariable('time', 'f4', ('time',))
    timevar.units = 'Daily'
    timevar.standard_name = 'time'

    lat_var = nc_fileo.createVariable('latitude', 'f8', ('latitude'),
                                         fill_value=fill_val)
    lat_var.units = 'degrees_north'
    lat_var.standard_name = 'latitude'

    lon_var = nc_fileo.createVariable('longitude', 'f8', ('longitude'),
                                         fill_value=fill_val)
    lon_var.units = 'degrees_east'
    lon_var.standard_name = 'longitude'

    cloud = nc_fileo.createVariable('cloud_mask', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    cloud.long_name = 'cloud mask'
    cloud.add_offset = 0.00
    cloud.grid_mapping = 'crs'
    cloud.scale_factor = Scaling_factor
    cloud.set_auto_maskandscale(False)
    reclin = nc_fileo.createVariable('reconst_lin', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    reclin.long_name = 'reconst brut'
    reclin.add_offset = 0.00
    reclin.grid_mapping = 'crs'
    reclin.scale_factor = Scaling_factor
    reclin.set_auto_maskandscale(False)
    recbrt = nc_fileo.createVariable('reconst_brt', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    recbrt.long_name = 'reconst brut'
    recbrt.add_offset = 0.00
    recbrt.grid_mapping = 'crs'
    recbrt.scale_factor = Scaling_factor
    recbrt.set_auto_maskandscale(False)
    rmsebrt = nc_fileo.createVariable('rmse_brt', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmsebrt.long_name = 'rmse brut'
    rmsebrt.add_offset = 0.00
    rmsebrt.grid_mapping = 'crs'
    rmsebrt.scale_factor = Scaling_factor
    rmsebrt.set_auto_maskandscale(False)
    rmselin = nc_fileo.createVariable('rmse_lin', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmselin.long_name = 'rmse lin'
    rmselin.add_offset = 0.00
    rmselin.grid_mapping = 'crs'
    rmselin.scale_factor = Scaling_factor
    rmselin.set_auto_maskandscale(False)
    rmsehants = nc_fileo.createVariable('rmse_hants', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmsehants.long_name = 'rmse hants'
    rmsehants.add_offset = 0.00
    rmsehants.grid_mapping = 'crs'
    rmsehants.scale_factor = Scaling_factor
    rmsehants.set_auto_maskandscale(False)

    timevar[:] = np.float_(time)
    lat_var[:] = latx
    lon_var[:] = lonx

    os.chdir(output_folder)
    files_nc = glob.glob('KARIMET_X*_Y*.nc')


    parameter = 'cloud_mask'
    Array = Get_Array(output_folder, nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    cloud[:,:,:] = Array
    del Array


    parameter = 'reconst_lin'
    Array = Get_Array(output_folder, nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    reclin[:,:,:] = Array
    del Array


    parameter = 'reconst_brt'
    Array = Get_Array(output_folder, nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    recbrt[:,:,:] = Array
    del Array

    parameter = 'rmse_brt'
    Array = Get_Array(output_folder, nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    rmsebrt[:,:] = Array
    del Array

    parameter = 'rmse_lin'
    Array = Get_Array(output_folder, nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    rmselin[:,:] = Array
    del Array

    parameter = 'rmse_hants'
    Array = Get_Array(output_folder,  nc_path_lst, files_nc, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    rmsehants[:,:] = Array
    del Array

    nc_fileo.close()

def Get_Array(output_folder, nc_path, files_nc, parameter):

    nc_file = netCDF4.Dataset(nc_path, 'r')
    latx = nc_file.variables['latitude'][:]
    lonx = nc_file.variables['longitude'][:]
    timeo = nc_file.variables['time'][:]

    if parameter[0:4] == 'rmse':
        Array = np.ones([len(latx), len(lonx)]) * np.nan
    else:
        Array = np.ones([len(timeo), len(latx), len(lonx)]) * np.nan

    for file_nc in files_nc:
        nc_file_part = netCDF4.Dataset(os.path.join(output_folder,file_nc))

        latx_part = nc_file_part.variables['latitude'][:]
        lonx_part = nc_file_part.variables['longitude'][:]
        start_y = np.argwhere(latx == latx_part[0])[0][0]
        start_x = np.argwhere(lonx == lonx_part[0])[0][0]
        if parameter[0:4] == 'rmse':
            Array_part = nc_file_part.variables[parameter][:,:]
            Array[start_y: start_y + Array_part.shape[0],start_x: start_x + Array_part.shape[1]] = Array_part
        else:
            Array_part = nc_file_part.variables[parameter][:,:,:]
            Array[:, start_y: start_y + Array_part.shape[1],start_x: start_x + Array_part.shape[2]] = Array_part

    return(Array)


def KARIMET(output_folder, nc_path, nc_path_GLDAS, CT, x_start, x_end, y_start, y_end, name):

    # HANTS output modis
    nc_file = netCDF4.Dataset(nc_path, 'r')
    latx = nc_file.variables['latitude'][y_start:y_end]
    lonx = nc_file.variables['longitude'][x_start:x_end]
    cloud=numpy.empty((len(latx),len(lonx)))


    nc_file2 = netCDF4.Dataset(nc_path_GLDAS, 'r')
    lat = nc_file2.variables['latitude'][:]
    lon = nc_file2.variables['longitude'][:]
    time = nc_file2.variables['time'][:]

    ########################## Open arrays ##########################################

    # open MODIS data
    values_modis_o_all = nc_file.variables['original_values'][:, y_start:y_end,x_start:x_end]
    values_modis_h_all = nc_file.variables['hants_values'][:, y_start:y_end,x_start:x_end]
    values_modis_o_all[values_modis_o_all == -9999] = np.nan
    values_modis_h_all[values_modis_h_all == -9999] = np.nan
    
    # open GLDAS data and fill missing values
    values_gldas_o_all = nc_file2.variables['original_values'][:,:,:]
    values_gldas_o_all[np.isnan(values_gldas_o_all)] =-9999
    for i in range(0,int(np.shape(values_gldas_o_all)[2])):
        values_gldas_o_all[i,:,:] = gap_filling(values_gldas_o_all[i,:,:],-9999)

    values_gldas_h_all = nc_file2.variables['hants_values'][:,:,:]
    values_gldas_h_all.mask = False
    for i in range(0,int(np.shape(values_gldas_h_all)[2])):
        values_gldas_h_all[i,:,:] = gap_filling(values_gldas_h_all[i,:,:],-9999)

    values_modis_o_all = values_modis_o_all - 273.15
    values_modis_h_all = values_modis_h_all - 273.15

    dT_modis = values_modis_o_all - values_modis_h_all
    val_dmax = values_gldas_o_all - values_gldas_h_all

    ################# resize val_dmax array ##########################################
    proj = "WGS84"
    size_x = (lon[0] - lon[-1])/(len(lon) - 1)
    size_y = (lat[0] - lat[-1])/(len(lat) - 1)
    geo = tuple([lon[0] + 0.5 * size_x, -size_x, 0 ,lat[0] - 0.5 * size_y,0,-size_y])

    size_x_example = (lonx[0] - lonx[-1])/(len(lonx) - 1)
    size_y_example = (latx[0] - latx[-1])/(len(latx) - 1)
    geo_example = tuple([lonx[0] + 0.5 * size_x_example, -size_x_example, 0 ,latx[-1] - 0.5 * size_y_example,0,size_y_example])
    dest_example = Save_as_MEM(values_modis_o_all[0,:,:], geo_example, proj)

    val_dmax_resized = np.ones(values_modis_o_all.shape) * np.nan

    for i in range(0,int(np.shape(values_gldas_h_all)[0])):

        dest_in_h = Save_as_MEM(val_dmax[i,:,:], geo, proj)
        dest_out_h = reproject_dataset_example(dest_in_h, dest_example)
        Array_out_h = dest_out_h.GetRasterBand(1).ReadAsArray()
        val_dmax_resized[i,:,:]=Array_out_h

    ###################### Create netcdf ##########################################

    epsg=4326
    fill_val=-9999.0
    spa_ref = Spatial_Reference(epsg)
    Scaling_factor = 0.0001

    fill_val = -9999
    nc_fileo = netCDF4.Dataset(os.path.join(output_folder,name), 'w', format="NETCDF4")
    lat_n = len(latx)
    lon_n = len(lonx)
    nc_fileo.createDimension('latitude', lat_n)
    nc_fileo.createDimension('longitude', lon_n)
    nc_fileo.createDimension('time', None)

    crs_var = nc_fileo.createVariable('crs', 'i4')
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.crs_wkt = spa_ref

    timeo = nc_fileo.createVariable('time', 'f4', ('time',))
    timeo.units = 'Daily'
    timeo.standard_name = 'time'

    lat_var = nc_fileo.createVariable('latitude', 'f8', ('latitude'),
                                         fill_value=fill_val)
    lat_var.units = 'degrees_north'
    lat_var.standard_name = 'latitude'
    lon_var = nc_fileo.createVariable('longitude', 'f8', ('longitude'),
                                         fill_value=fill_val)
    lon_var.units = 'degrees_east'
    lon_var.standard_name = 'longitude'
    cloud = nc_fileo.createVariable('cloud_mask', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    cloud.long_name = 'cloud mask'
    cloud.add_offset = 0.00
    cloud.grid_mapping = 'crs'
    cloud.scale_factor = Scaling_factor
    cloud.set_auto_maskandscale(False)
    reclin = nc_fileo.createVariable('reconst_lin', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    reclin.long_name = 'reconst brut'
    reclin.add_offset = 0.00
    reclin.grid_mapping = 'crs'
    reclin.scale_factor = Scaling_factor
    reclin.set_auto_maskandscale(False)
    recbrt = nc_fileo.createVariable('reconst_brt', 'f8',
                                              ('time','latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    recbrt.long_name = 'reconst brut'
    recbrt.add_offset = 0.00
    recbrt.grid_mapping = 'crs'
    recbrt.scale_factor = Scaling_factor
    recbrt.set_auto_maskandscale(False)
    rmsebrt = nc_fileo.createVariable('rmse_brt', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmsebrt.long_name = 'rmse brut'
    rmsebrt.add_offset = 0.00
    rmsebrt.grid_mapping = 'crs'
    rmsebrt.scale_factor = Scaling_factor
    rmsebrt.set_auto_maskandscale(False)
    rmselin = nc_fileo.createVariable('rmse_lin', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmselin.long_name = 'rmse lin'
    rmselin.add_offset = 0.00
    rmselin.grid_mapping = 'crs'
    rmselin.scale_factor = Scaling_factor
    rmselin.set_auto_maskandscale(False)
    rmsehants = nc_fileo.createVariable('rmse_hants', 'f8',
                                              ('latitude', 'longitude'),fill_value=-9999,
                                        zlib=True, least_significant_digit=0)
    rmsehants.long_name = 'rmse hants'
    rmsehants.add_offset = 0.00
    rmsehants.grid_mapping = 'crs'
    rmsehants.scale_factor = Scaling_factor
    rmsehants.set_auto_maskandscale(False)

    ################ Apply linear fit and reconstruct data #############################

    val_reconst = val_dmax_resized + values_modis_h_all

    values_modis_o_corr=values_modis_o_all
    values_modis_o_corr= np.ma.masked_where(dT_modis<CT, values_modis_o_all)

    val_dmax_corr=np.ma.masked_where(dT_modis<CT, val_dmax_resized)
    val_dmax_corr=np.ma.masked_where(val_dmax_corr<-50., val_dmax_corr)

    dT_modis_corr=np.ma.masked_where(dT_modis<CT, dT_modis)

    X = val_dmax_corr
    X.mask = np.nan    
    X.mask = False

    Y = dT_modis_corr
    Y.mask = np.nan
    Y.mask = False

    Y[Y<-9000] = np.nan
    X[X<-9000] = np.nan

    m = ((np.sum(np.where(np.isnan(X),0,1),axis = 0) * np.nansum(X * Y, axis = 0)) - (np.nansum(X, axis = 0) * np.nansum(Y, axis = 0)))/((np.sum(np.where(np.isnan(X),0,1),axis = 0)* np.nansum(X * X, axis = 0)) - (np.nansum(X, axis = 0) * np.nansum(X, axis = 0)))
    c = ((np.nansum(X * X, axis = 0) * np.nansum(Y, axis = 0)) - (np.nansum(X, axis = 0) * np.nansum(X * Y, axis = 0)))/((np.sum(np.where(np.isnan(X),0,1),axis = 0) * np.nansum(X * X, axis = 0)) - (np.nansum(X, axis = 0) * np.nansum(X, axis = 0)))

    yy_new = m[None,:,:] * val_dmax_resized + c[None,:,:]

    val_reconst_corr = yy_new + values_modis_h_all

    cloud_end = values_modis_o_all-values_modis_h_all
    reclin_end = np.where(dT_modis<CT, val_reconst_corr,values_modis_o_all) 
    reclin_end[reclin_end<-100] = val_reconst_corr[reclin_end<-100]
    reclin_end[np.isnan(reclin_end)] = val_reconst_corr[np.isnan(reclin_end)]    
    recbrt_end =values_modis_h_all+val_dmax_resized
    rmsebrt_end =np.sqrt(np.nanmean((values_modis_o_corr-val_reconst)**2, axis = 0))
    rmselin_end =np.sqrt(np.nanmean((values_modis_o_corr-val_reconst_corr)**2, axis = 0))
    rmsehants_end =np.sqrt(np.nanmean((values_modis_o_corr-values_modis_h_all)**2, axis = 0))

    ############################# Save the data ###################################

    cloud_end.mask = False
    recbrt_end.mask = False
    rmsebrt_end.mask = False
    rmselin_end.mask = False
    rmsehants_end.mask = False

    cloud_end[np.isnan(cloud_end)] = -9999 * np.float(Scaling_factor)
    cloud_end = np.int_(cloud_end * 1./np.float(Scaling_factor))
    recbrt_end[np.isnan(recbrt_end)] = -9999 * np.float(Scaling_factor)
    recbrt_end = np.int_(recbrt_end * 1./np.float(Scaling_factor))
    reclin_end[np.isnan(reclin_end)] = -9999 * np.float(Scaling_factor)
    reclin_end = np.int_(reclin_end * 1./np.float(Scaling_factor))
    rmsebrt_end[np.isnan(rmsebrt_end)] = -9999 * np.float(Scaling_factor)
    rmsebrt_end = np.int_(rmsebrt_end * 1./np.float(Scaling_factor))
    rmselin_end[np.isnan(rmselin_end)] = -9999 * np.float(Scaling_factor)
    rmselin_end = np.int_(rmselin_end * 1./np.float(Scaling_factor))
    rmsehants_end[np.isnan(rmsehants_end)] = -9999 * np.float(Scaling_factor)
    rmsehants_end = np.int_(rmsehants_end * 1./np.float(Scaling_factor))

    timeo[:] = np.float_(time)

    lat_var[:] = latx
    lon_var[:] = lonx

    for i in range(len(time)):
        cloud[i,:,:] = cloud_end[i,:,:]
        reclin[i,:,:] = reclin_end[i,:,:]
        recbrt[i,:,:] = recbrt_end[i,:,:]
    rmsebrt[:,:] = rmsebrt_end[:,:]
    rmselin[:,:] = rmselin_end[:,:]
    rmsehants[:,:] = rmsehants_end[:,:]

    nc_fileo.close()

########################## Functions ##########################################

def gap_filling(data,NoDataValue, method = 1):
    """
    This function fills the no data gaps in a numpy array

    Keyword arguments:
    dataset -- 'C:/'  path to the source data (dataset that must be filled)
    NoDataValue -- Value that must be filled
    """

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


    EndProduct = data_end

    return (EndProduct)

def Get_epsg(g, extension = 'tiff'):
    """
    This function reads the projection of a GEOGCS file or tiff file

    Keyword arguments:
    g -- string
        Filename to the file that must be read
    extension -- tiff or GEOGCS
        Define the extension of the dataset (default is tiff)
    """
    try:
        if extension == 'tiff':
            # Get info of the dataset that is used for transforming
            g_proj = g.GetProjection()
            Projection=g_proj.split('EPSG","')
        if extension == 'GEOGCS':
            Projection = g
        epsg_to=int((str(Projection[-1]).split(']')[0])[0:-1])
    except:
       epsg_to=4326
       #print 'Was not able to get the projection, so WGS84 is assumed'
    return(epsg_to)

def reproject_dataset_example(dataset, dataset_example, method=1):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The user can define the wanted projection and shape by defining an example dataset.

    Keywords arguments:
    dataset -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    dataset_example -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    method -- 1,2,3,4 default = 1
        1 = Nearest Neighbour, 2 = Bilinear, 3 = lanzcos, 4 = average
    """
    # open dataset that must be transformed
    try:
        if os.path.splitext(dataset)[-1] == '.tif':
            g = gdal.Open(dataset)
        else:
            g = dataset
    except:
            g = dataset
    epsg_from = Get_epsg(g)

    #exceptions
    if epsg_from == 9001:
        epsg_from = 5070

    # open dataset that is used for transforming the dataset
    gland = dataset_example
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
    if method is 1:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    if method is 2:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Bilinear)
    if method is 3:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Lanczos)
    if method is 4:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
    return(dest1)


def Save_as_MEM(data='', geo='', projection=''):
    """
    This function save the array as a memory file

    Keyword arguments:
    data -- [array], dataset of the geotiff
    geo -- [minimum lon, pixelsize, rotation, maximum lat, rotation,
            pixelsize], (geospatial dataset)
    projection -- interger, the EPSG code
    """
    # save as a geotiff
    driver = gdal.GetDriverByName("MEM")
    dst_ds = driver.Create('', int(data.shape[1]), int(data.shape[0]), 1,
                           gdal.GDT_Float32)
    srse = osr.SpatialReference()
    if projection == '':
        srse.SetWellKnownGeogCS("WGS84")
    else:
        srse.SetWellKnownGeogCS(projection)
    dst_ds.SetProjection(srse.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo)
    dst_ds.GetRasterBand(1).WriteArray(data)
    return(dst_ds)


def Spatial_Reference(epsg, return_string=True):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    if return_string:
        return srs.ExportToWkt()
    else:
        return srs
