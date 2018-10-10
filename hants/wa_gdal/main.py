# -*- coding: utf-8 -*-
"""
Authors: Gonzalo E. Espinoza-Dávalos, Wim G.M. Bastiaanssen, Boaz Bett, and
         Xueliang Cai
         IHE Delft 2017
Contact: g.espinoza@un-ihe.org
Repository: https://github.com/gespinoza/hants
Module: hants
"""

from __future__ import division
import netCDF4
import pandas as pd
import numpy as np
import datetime
import math
import os
import osr
import glob
from copy import deepcopy
import matplotlib.pyplot as plt
import warnings
import gdal
from joblib import Parallel, delayed

def run_HANTS(rasters_path_inp, name_format,
              start_date, end_date, latlim, lonlim, cellsize, nc_path,
              nb, nf, HiLo, low, high, fet, dod, delta, Scaling_factor = 0.001,
              epsg=4326, cores=1):
    '''
    This function runs the python implementation of the HANTS algorithm. It
    takes a folder with geotiffs raster data as an input, creates a netcdf
    file, and optionally export the data back to geotiffs.
    '''
    nc_paths = create_netcdf(rasters_path_inp, name_format, start_date, end_date,
                  latlim, lonlim, cellsize, nc_path, Scaling_factor,
                  epsg)
    args = [nb, nf, HiLo, low, high, fet, dod, delta, Scaling_factor]
    print('\tApply HANTS on tiles...')
    results = Parallel(n_jobs=cores)(delayed(HANTS_netcdf)(nc_path, args)
                                         for nc_path in nc_paths)

    if len(nc_paths) > 1:
        Merge_NC_Tiles(nc_paths, nc_path, start_date, end_date, latlim, lonlim, cellsize, epsg, Scaling_factor)

    return nc_path


def create_netcdf(rasters_path, name_format, start_date, end_date,
                  latlim, lonlim, cellsize, nc_path, Scaling_factor,
                  epsg=4326):
    '''
    This function creates a netcdf file from a folder with geotiffs rasters to
    be used to run HANTS.
    '''
    # Latitude and longitude
    lat_ls = pd.np.arange(latlim[0] + 0.5*cellsize, latlim[1],
                          cellsize)
    lat_ls = lat_ls[::-1]  # ArcGIS numpy
    lon_ls = pd.np.arange(lonlim[0] + 0.5*cellsize, lonlim[1],
                          cellsize)
    lat_n = len(lat_ls)
    lon_n = len(lon_ls)
    spa_ref = Spatial_Reference(epsg)
    # ll_corner = [lonlim[0], latlim[0]]

    # Rasters
    dates_dt = pd.date_range(start_date, end_date, freq='D')
    dates_ls = [d.toordinal() for d in dates_dt]
    os.chdir(rasters_path)
    ras_ls = glob.glob('*.tif')

    # Create tile parts
    if (lat_n > 200 or lon_n > 200):

        lat_n_amount = np.maximum(1,int(np.floor(lat_n/100)))
        lon_n_amount = np.maximum(1,int(np.floor(lon_n/100)))
        nc_path_part_names = nc_path.split('.')
        nc_path_tiles = []
        for lat_n_one in range(0, lat_n_amount):
            for lon_n_one in range(0, lon_n_amount):
                nc_path_tile = ''.join(nc_path_part_names[0] + "_h%03d_v%03d.nc" %(lon_n_one, lat_n_one))
                nc_path_tiles = np.append(nc_path_tiles, nc_path_tile)

    else:
        nc_path_tiles = nc_path

    i = 0
    # Loop over the nc_paths
    for nc_path_tile in nc_path_tiles:
        i += 1
        if lat_n_amount > 1:
            lat_part = int(nc_path_tile[-6:-3])
            
            lat_start = lat_part * 100
            if int(lat_part) is not int(lat_n_amount-1):
                lat_end = int((lat_part + 1) * 100)
            else:
                lat_end = int(lat_n)
        else:
            lat_start = int(0)
            lat_end = int(lat_n)
         
        if lon_n_amount > 1:
            lon_part = int(nc_path_tile[-11:-8])

            lon_start = int(lon_part * 100)
            if int(lon_part) is not int(lon_n_amount-1):
                lon_end = int((lon_part + 1) * 100)
            else:
                lon_end = int(lon_n)
        else:
            lon_start = int(0)
            lon_end = int(lon_n)


        # Define space dimention
        lat_range = lat_ls[lat_start:lat_end]
        lon_range = lon_ls[lon_start:lon_end]
        geo_ex = tuple([lon_range[0] - 0.5*cellsize, cellsize, 0, lat_range[0] + cellsize * 0.5, 0, -cellsize])

        # Create netcdf file
        print('Creating netCDF file tile %s out of %s...' %(i,len(nc_path_tiles)))
        nc_file = netCDF4.Dataset(nc_path_tile, 'w', format="NETCDF4_CLASSIC")

        # Create Dimensions
        lat_dim = nc_file.createDimension('latitude', lat_end - lat_start)
        lon_dim = nc_file.createDimension('longitude', lon_end - lon_start)
        time_dim = nc_file.createDimension('time', len(dates_ls))

        # Create Variables
        crso = nc_file.createVariable('crs', 'i4')
        crso.long_name = 'Lon/Lat Coords in WGS84'
        crso.standard_name = 'crs'
        crso.grid_mapping_name = 'latitude_longitude'
        crso.projection = spa_ref
        crso.longitude_of_prime_meridian = 0.0
        crso.semi_major_axis = 6378137.0
        crso.inverse_flattening = 298.257223563
        crso.geo_reference = geo_ex

        lat_var = nc_file.createVariable('latitude', 'f8', ('latitude',))
        lat_var.units = 'degrees_north'
        lat_var.standard_name = 'latitude'

        lon_var = nc_file.createVariable('longitude', 'f8', ('longitude',))
        lon_var.units = 'degrees_east'
        lon_var.standard_name = 'longitude'

        time_var = nc_file.createVariable('time', 'l', ('time',))
        time_var.standard_name = 'time'
        time_var.calendar = 'gregorian'

        original_var = nc_file.createVariable('original_values', 'i',
                                              ('time', 'latitude', 'longitude'),
                                              fill_value=-9999, zlib=True, least_significant_digit=0)
        original_var.long_name = 'original_values'
        original_var.grid_mapping = 'crs'
        original_var.add_offset = 0.00
        original_var.scale_factor = Scaling_factor
        original_var.set_auto_maskandscale(False)
        print('\tVariables created')

        # Fill in time and space dimention
        lat_var[:] = lat_range
        lon_var[:] = lon_range
        time_var[:] = dates_ls

        # Create memory example file
        # empty array
        empty_vec = pd.np.empty((lat_end - lat_start, lon_end - lon_start))
        empty_vec[:] = -9999 * np.float(Scaling_factor)
        dest_ex = Save_as_MEM(empty_vec, geo_ex, str(epsg))

        # Raster loop
        print('\tExtracting data from rasters...')
        for tt in range(len(dates_ls)):

            Date_now = datetime.datetime.fromordinal(dates_ls[tt])
            yyyy = str(Date_now.year)
            mm = '%02d' %int(Date_now.month)            
            dd = '%02d' %int(Date_now.day)
             
            # Raster
            ras = name_format.format(yyyy=yyyy,mm=mm,dd=dd)

            if ras in ras_ls:

                data_in = os.path.join(rasters_path, ras)
                dest = reproject_dataset_example(data_in, dest_ex)
                array_tt = dest.GetRasterBand(1).ReadAsArray()
                array_tt[array_tt<-9999] = -9999 * np.float(Scaling_factor)
                original_var[tt, :, :] = np.int_(array_tt * 1./np.float(Scaling_factor))

            else:
                # Store values
                original_var[tt, :, :] = np.int_(empty_vec * 1./np.float(Scaling_factor))

        # Close file
        nc_file.close()
        print('NetCDF %s file created' %i)

    # Return
    return nc_path_tiles


def HANTS_netcdf(nc_path, args):
    '''
    This function runs the python implementation of the HANTS algorithm. It
    takes the input netcdf file and fills the 'hants_values',
    'combined_values', and 'outliers' variables.
    '''
    nb, nf, HiLo, low, high, fet, dod, delta, Scaling_factor = args

    # Read netcdfs
    nc_file = netCDF4.Dataset(nc_path, 'r+', format="NETCDF4_CLASSIC")
    nc_file.set_fill_on()
    
    time_var = nc_file.variables['time'][:]
    original_values = nc_file.variables['original_values'][:]

    [ztime, rows, cols] = original_values.shape
    size_st = cols*rows

    values_hants = pd.np.empty((ztime, rows, cols))
    outliers_hants = pd.np.empty((ztime, rows, cols))

    values_hants[:] = pd.np.nan
    outliers_hants[:] = pd.np.nan

    # Additional parameters
    ni = len(time_var)
    ts = range(ni)

    # Loop
    counter = 1
    #print('Running HANTS...')
    for m in range(rows):
        for n in range(cols):
            #print('\t{0}/{1}'.format(counter, size_st))

            y = pd.np.array(original_values[:, m, n])

            y[pd.np.isnan(y)] = -9999

            [yr, outliers] = HANTS(ni, nb, nf, y, ts, HiLo,
                                   low, high, fet, dod, delta)

            values_hants[:, m, n] = yr
            outliers_hants[:, m, n] = outliers

            counter = counter + 1
            
    values_hants[values_hants<-9999] = -9999 * np.float(Scaling_factor)

   
    hants_var = nc_file.createVariable('hants_values', 'i',
                                       ('time', 'latitude', 'longitude'),
                                       fill_value=-9999, zlib=True, least_significant_digit=0)
    hants_var.long_name = 'hants_values'
    hants_var.grid_mapping = 'crs'
    hants_var.add_offset = 0.00
    hants_var.scale_factor = Scaling_factor
    hants_var.set_auto_maskandscale(False)

    combined_var = nc_file.createVariable('combined_values', 'i',
                                          ('time', 'latitude', 'longitude'),
                                          fill_value=-9999, zlib=True, least_significant_digit=0)
    combined_var.long_name = 'combined_values'
    combined_var.grid_mapping = 'crs'
    combined_var.add_offset = 0.00
    combined_var.scale_factor = Scaling_factor
    combined_var.set_auto_maskandscale(False)   
    
    outliers_var = nc_file.createVariable('outliers', 'i4',
                                          ('time', 'latitude', 'longitude'),
                                          fill_value=-9999)
    outliers_var.long_name = 'outliers'
    outliers_var.grid_mapping = 'crs'    
    
    hants_var[:,:,:]= np.int_(values_hants * 1./np.float(Scaling_factor))
    outliers_var[:,:,:] = outliers_hants
    combined_var[:,:,:] = pd.np.where(outliers_hants,
                                      np.int_(values_hants * 1./np.float(Scaling_factor)),
                                      np.int_(original_values * 1./np.float(Scaling_factor)))
    # Close netcdf file
    nc_file.close()


def HANTS_singlepoint(nc_path, point, nb, nf, HiLo, low, high, fet, dod,
                      delta):
    '''
    This function runs the python implementation of the HANTS algorithm for a
    single point (lat, lon). It plots the fit and returns a data frame with
    the 'original' and the 'hants' time series.
    '''
    # Location
    lonx = point[0]
    latx = point[1]

    nc_file = netCDF4.Dataset(nc_path, 'r', format="NETCDF4_CLASSIC")

    time = [pd.to_datetime(i, format='%Y%m%d')
            for i in nc_file.variables['time'][:]]

    lat = nc_file.variables['latitude'][:]
    lon = nc_file.variables['longitude'][:]

    # Check that the point falls within the extent of the netcdf file
    lon_max = max(lon)
    lon_min = min(lon)
    lat_max = max(lat)
    lat_min = min(lat)
    if not (lon_min < lonx < lon_max) or not (lat_min < latx < lat_max):
        warnings.warn('The point lies outside the extent of the netcd file. '
                      'The closest cell is plotted.')
        if lonx > lon_max:
            lonx = lon_max
        elif lonx < lon_min:
            lonx = lon_min
        if latx > lat_max:
            latx = lat_max
        elif latx < lat_min:
            latx = lat_min

    # Get lat-lon index in the netcdf file
    lat_closest = lat.flat[pd.np.abs(lat - latx).argmin()]
    lon_closest = lon.flat[pd.np.abs(lon - lonx).argmin()]

    lat_i = pd.np.where(lat == lat_closest)[0][0]
    lon_i = pd.np.where(lon == lon_closest)[0][0]

    # Read values
    original_values = nc_file.variables['original_values'][:, lat_i, lon_i]

    # Additional parameters
    ni = len(time)
    ts = range(ni)

    # HANTS
    y = pd.np.array(original_values)

    y[pd.np.isnan(y)] = -9999

    [hants_values, outliers] = HANTS(ni, nb, nf, y, ts, HiLo, low, high, fet,
                                     dod, delta)
    # Plot
    top = 1.15*max(pd.np.nanmax(original_values),
                   pd.np.nanmax(hants_values))
    bottom = 1.15*min(pd.np.nanmin(original_values),
                      pd.np.nanmin(hants_values))
    ylim = [bottom, top]

    plt.plot(time, hants_values, 'r-', label='HANTS')
    plt.plot(time, original_values, 'b.', label='Original data')

    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc=4)
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: lon {0:.2f}, lat {1:.2f}'.format(lon_closest,
                                                                  lat_closest))
    plt.axes().set_aspect(0.5*(time[-1] - time[0]).days/(ylim[1] - ylim[0]))

    plt.show()

    # Close netcdf file
    nc_file.close()

    # Data frame
    df = pd.DataFrame({'time': time,
                       'original': original_values,
                       'hants': hants_values})

    # Return
    return df


def HANTS(ni, nb, nf, y, ts, HiLo, low, high, fet, dod, delta):
    '''
    This function applies the Harmonic ANalysis of Time Series (HANTS)
    algorithm originally developed by the Netherlands Aerospace Centre (NLR)
    (http://www.nlr.org/space/earth-observation/).

    This python implementation was based on two previous implementations
    available at the following links:
    https://codereview.stackexchange.com/questions/71489/harmonic-analysis-of-time-series-applied-to-arrays
    http://nl.mathworks.com/matlabcentral/fileexchange/38841-matlab-implementation-of-harmonic-analysis-of-time-series--hants-
    '''
    # Arrays
    mat = pd.np.zeros((min(2*nf+1, ni), ni))
    # amp = np.zeros((nf + 1, 1))

    # phi = np.zeros((nf+1, 1))
    yr = pd.np.zeros((ni, 1))
    outliers = pd.np.zeros((1, len(y)))

    # Filter
    sHiLo = 0
    if HiLo == 'Hi':
        sHiLo = -1
    elif HiLo == 'Lo':
        sHiLo = 1

    nr = min(2*nf+1, ni)
    noutmax = ni - nr - dod
    # dg = 180.0/math.pi
    mat[0, :] = 1.0

    ang = 2*math.pi*pd.np.arange(nb)/nb
    cs = pd.np.cos(ang)
    sn = pd.np.sin(ang)

    i = pd.np.arange(1, nf+1)
    for j in pd.np.arange(ni):
        index = pd.np.mod(i*ts[j], nb)
        mat[2 * i-1, j] = cs.take(index)
        mat[2 * i, j] = sn.take(index)

    p = pd.np.ones_like(y)
    bool_out = (y < low) | (y > high)
    p[bool_out] = 0
    outliers[bool_out.reshape(1, y.shape[0])] = 1
    nout = pd.np.sum(p == 0)

    if nout > noutmax:
        if pd.np.isclose(y, -9999).any():
            ready = pd.np.array([True])
            yr = y
            outliers = pd.np.zeros((y.shape[0]), dtype=int)
            outliers[:] = -9999
        else:
            raise Exception('Not enough data points.')
    else:
        ready = pd.np.zeros((y.shape[0]), dtype=bool)

    nloop = 0
    nloopmax = ni

    while ((not ready.all()) & (nloop < nloopmax)):

        nloop += 1
        za = pd.np.matmul(mat, p*y)

        A = pd.np.matmul(pd.np.matmul(mat, pd.np.diag(p)),
                         pd.np.transpose(mat))
        A = A + pd.np.identity(nr)*delta
        A[0, 0] = A[0, 0] - delta

        zr = pd.np.linalg.solve(A, za)

        yr = pd.np.matmul(pd.np.transpose(mat), zr)
        diffVec = sHiLo*(yr-y)
        err = p*diffVec

        err_ls = list(err)
        err_sort = deepcopy(err)
        err_sort.sort()

        rankVec = [err_ls.index(f) for f in err_sort]

        maxerr = diffVec[rankVec[-1]]
        ready = (maxerr <= fet) | (nout == noutmax)

        if (not ready):
            i = ni - 1
            j = rankVec[i]
            while ((p[j]*diffVec[j] > 0.5*maxerr) & (nout < noutmax)):
                p[j] = 0
                outliers[0, j] = 1
                nout += 1
                i -= 1
                if i == 0:
                    j = 0
                else:
                    j = 1

    return [yr, outliers]


def plot_point(nc_path, point, ylim=None):
    '''
    This function plots the original time series and the HANTS time series.
    It can be used to assess the fit.
    '''
    # Location
    lonx = point[0]
    latx = point[1]

    nc_file = netCDF4.Dataset(nc_path, 'r', format="NETCDF4_CLASSIC")

    time = [pd.to_datetime(i, format='%Y%m%d')
            for i in nc_file.variables['time'][:]]

    lat = nc_file.variables['latitude'][:]
    lon = nc_file.variables['longitude'][:]

    # Check that the point falls within the extent of the netcdf file
    lon_max = max(lon)
    lon_min = min(lon)
    lat_max = max(lat)
    lat_min = min(lat)
    if not (lon_min < lonx < lon_max) or not (lat_min < latx < lat_max):
        warnings.warn('The point lies outside the extent of the netcd file. '
                      'The closest cell is plotted.')
        if lonx > lon_max:
            lonx = lon_max
        elif lonx < lon_min:
            lonx = lon_min
        if latx > lat_max:
            latx = lat_max
        elif latx < lat_min:
            latx = lat_min

    # Get lat-lon index in the netcdf file
    lat_closest = lat.flat[pd.np.abs(lat - latx).argmin()]
    lon_closest = lon.flat[pd.np.abs(lon - lonx).argmin()]

    lat_i = pd.np.where(lat == lat_closest)[0][0]
    lon_i = pd.np.where(lon == lon_closest)[0][0]

    # Read values
    values_o = nc_file.variables['original_values'][lat_i, lon_i, :]
    values_h = nc_file.variables['hants_values'][lat_i, lon_i, :]

    if not ylim:
        top = 1.15*max(pd.np.nanmax(values_o),
                       pd.np.nanmax(values_h))
        bottom = 1.15*min(pd.np.nanmin(values_o),
                          pd.np.nanmin(values_h))
        ylim = [bottom, top]

    # Plot
    plt.plot(time, values_h, 'r-', label='HANTS')
    plt.plot(time, values_o, 'b.', label='Original data')

    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc=4)
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: lon {0:.2f}, lat {1:.2f}'.format(lon_closest,
                                                                  lat_closest))
    plt.axes().set_aspect(0.5*(time[-1] - time[0]).days/(ylim[1] - ylim[0]))

    plt.show()

    # Close netcdf file
    nc_file.close()

    # Return
    return True


def Merge_NC_Tiles(nc_paths, nc_path, start_date, end_date, latlim, lonlim, cellsize, epsg, Scaling_factor):


    # Latitude and longitude
    lat_ls = pd.np.arange(latlim[0] + 0.5*cellsize, latlim[1],
                          cellsize)
    lat_ls = lat_ls[::-1]  # ArcGIS numpy
    lon_ls = pd.np.arange(lonlim[0] + 0.5*cellsize, lonlim[1],
                          cellsize)
    lat_n = len(lat_ls)
    lon_n = len(lon_ls)
    spa_ref = Spatial_Reference(epsg)
    geo_ex = tuple([lon_ls[0] - 0.5*cellsize, cellsize, 0, lat_ls[0] - cellsize * 0.5, 0, -cellsize])
    dates_dt = pd.date_range(start_date, end_date, freq='D')
    dates_ls = [d.toordinal() for d in dates_dt]

    # Create netcdf file
    print('Merging netCDF files...')
    nc_file = netCDF4.Dataset(nc_path, 'w', format="NETCDF4_CLASSIC")
    nc_file.set_fill_on()
    
    # Create Dimensions
    lat_dim = nc_file.createDimension('latitude', lat_n)
    lon_dim = nc_file.createDimension('longitude', lon_n)
    time_dim = nc_file.createDimension('time', len(dates_ls))

    # Create Variables
    crso = nc_file.createVariable('crs', 'i4')
    crso.long_name = 'Lon/Lat Coords in WGS84'
    crso.standard_name = 'crs'
    crso.grid_mapping_name = 'latitude_longitude'
    crso.projection = spa_ref
    crso.longitude_of_prime_meridian = 0.0
    crso.semi_major_axis = 6378137.0
    crso.inverse_flattening = 298.257223563
    crso.geo_reference = geo_ex

    lat_var = nc_file.createVariable('latitude', 'f8', ('latitude',))
    lat_var.units = 'degrees_north'
    lat_var.standard_name = 'latitude'
    lat_var.pixel_size = cellsize
    
    lon_var = nc_file.createVariable('longitude', 'f8', ('longitude',))
    lon_var.units = 'degrees_east'
    lon_var.standard_name = 'longitude'
    lon_var.pixel_size = cellsize

    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.standard_name = 'time'
    time_var.calendar = 'gregorian'

    outliers_var = nc_file.createVariable('outliers', 'i4',
                                          ('time', 'latitude', 'longitude'),
                                          fill_value=-9999)
    outliers_var.long_name = 'outliers'
    outliers_var.grid_mapping = 'crs'

    original_var = nc_file.createVariable('original_values', 'i',
                                          ('time', 'latitude', 'longitude'),
                                          fill_value=-9999, zlib=True, least_significant_digit=0)
    original_var.long_name = 'original_values'
    original_var.grid_mapping = 'crs'  
    original_var.scale_factor = Scaling_factor
    original_var.add_offset = 0.00      
    original_var.set_auto_maskandscale(False)

    hants_var = nc_file.createVariable('hants_values', 'i',
                                       ('time', 'latitude', 'longitude'),
                                       fill_value=-9999, zlib=True, least_significant_digit=0)
    hants_var.long_name = 'hants_values'
    hants_var.grid_mapping = 'crs'
    hants_var.scale_factor = Scaling_factor
    hants_var.add_offset = 0.00    
    hants_var.set_auto_maskandscale(False)

    combined_var = nc_file.createVariable('combined_values', 'i',
                                          ('time', 'latitude', 'longitude'),
                                          fill_value=-9999, zlib=True, least_significant_digit=0)
    combined_var.long_name = 'combined_values'
    combined_var.grid_mapping = 'crs'
    combined_var.scale_factor = Scaling_factor
    combined_var.add_offset = 0.00    
    combined_var.set_auto_maskandscale(False)

    print('\tFill in End Variables')

    # Fill in time and space dimention
    lat_var[:] = lat_ls
    lon_var[:] = lon_ls
    time_var[:] = dates_ls

    parameter = 'outliers'
    Array = Get_Array(nc_path, nc_paths, parameter)
    Array[np.isnan(Array)] = -9999
    outliers_var[:,:,:] = np.int_(Array)
    del Array

    parameter = 'original_values'
    Array = Get_Array(nc_path, nc_paths, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array[Array < -9999] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    original_var[:,:,:] = Array
    del Array

    parameter = 'hants_values'
    Array = Get_Array(nc_path, nc_paths, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array[Array < -9999] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    hants_var[:,:,:] = Array
    del Array

    parameter = 'combined_values'
    Array = Get_Array(nc_path, nc_paths, parameter)
    Array[np.isnan(Array)] = -9999 * np.float(Scaling_factor)
    Array[Array < -9999] = -9999 * np.float(Scaling_factor)
    Array = np.int_(Array * 1./np.float(Scaling_factor))
    combined_var[:,:,:] = Array
    del Array

    nc_file.close()
    return()

def Get_Array(nc_path, nc_paths, parameter):

    nc_file = netCDF4.Dataset(nc_path, 'r', format="NETCDF4_CLASSIC")
    latx = nc_file.variables['latitude'][:]
    lonx = nc_file.variables['longitude'][:]
    timeo = nc_file.variables['time'][:]

    Array = np.ones([len(timeo), len(latx), len(lonx)]) * np.nan

    for file_nc in nc_paths:
        nc_file_part = netCDF4.Dataset(file_nc)

        latx_part = nc_file_part.variables['latitude'][:]
        lonx_part = nc_file_part.variables['longitude'][:]
        start_y = np.argwhere(latx == latx_part[0])[0][0]
        start_x = np.argwhere(lonx == lonx_part[0])[0][0]
        Array_part = nc_file_part.variables[parameter][:,:,:]
        Array[:, start_y: start_y + Array_part.shape[1],start_x: start_x + Array_part.shape[2]] = Array_part

    return(Array)

def makediag3d(M):
    '''
    Computing diagonal for each row of a 2d array.
    Reference: http://stackoverflow.com/q/27214027/2459096
    '''
    b = pd.np.zeros((M.shape[0], M.shape[1]*M.shape[1]))
    b[:, ::M.shape[1]+1] = M
    # Return
    return b.reshape(M.shape[0], M.shape[1], M.shape[1])

def Spatial_Reference(epsg, return_string=True):
    """
    Obtain a spatial reference from the EPSG parameter
    """
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    if return_string:
        return srs.ExportToWkt()
    else:
        return srs

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
