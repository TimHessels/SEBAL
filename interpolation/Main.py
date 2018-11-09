# -*- coding: utf-8 -*-
"""
This is the Main module and depends on modules 'Functions' and 'QGIS_'
It gets ETo values and multiplies it with Kc maps to get ETa values
"""

import Functions as fnx
import numpy as np
import copy
import os
import QGIS_ as qgis
#%%set environments

print("Choose folder containing the 'Kc' files")
input_folder = fnx.SelectFolderDialog()
print("Choose the folder to save the output")
output_folder =fnx.SelectFolderDialog()
print("select the 'ET0.csv' file")
ET0_daily = fnx.OpenFileDialog()
#%%
#input_folder = r'D:\OneDrive\Courses\Environmental Programming\CJ\kc'
#output_folder =  r'V:\ETa Program'
Kc_files = fnx.ListFilesInFolder (input_folder,'.tif')                                  #get list of kc files

while len(Kc_files) == 0:
    print("The selected folder does not have 'Kc' files \nPlease choose a folder with 'Kc' files")
    input_folder = fnx.SelectFolderDialog()
    Kc_files = fnx.ListFilesInFolder (input_folder,'.tif')

ET0_daily = fnx.import_txt(ET0_daily) #import the ET0 file

#%%
ET0_array = []                                                      #array for importing ET0 data
time_stamp_int_et0=[]                                               #array for timestamps
for i in ET0_daily:                                                 #looping through elements of the ET0_daily imported data array
    ET0_array.append(i.split(';')[1] +',' + i.split(';')[2])        #Getting array of ET0 data 
    time_stamp_int_et0.append(int(i.split(';')[1]))                 #populating array for timestamps

#%%  
# Time-stapms
time_stamp_int_kc=[]
time_stamp_str_kc=[]

for filename in Kc_files:
    time_step=filename[-7:-4]
    time_stamp_int_kc.append(int(time_step))
    time_stamp_str_kc.append(time_step)
#%%
#get window available for the calculations
minimum_kc = np.min(time_stamp_int_kc)
minimum_et0 = np.min(time_stamp_int_et0)
lower_time_limit=max(minimum_kc,minimum_et0)

maximum_kc = np.max(time_stamp_int_kc)
maximum_et0 = np.max(time_stamp_int_et0)
upper_time_limit=min(maximum_kc,maximum_et0)

#%%
#Make user chose the start date and end date for ETa output
print "Please select the start day for ETa Calculation"
starting_day = str(raw_input("Select between " + str(lower_time_limit) + " and " + str(upper_time_limit) + ": "))
try:
    start_day = int(starting_day)
except:
    start_day = -1

while not(lower_time_limit<=start_day<upper_time_limit):                                                             #checking that the selection is within range of available data
    if start_day == -1:
        print "Enter a number"
    elif start_day<lower_time_limit:
        print ("\n\nThe selected day is below the limit of " + str(lower_time_limit))       
    elif start_day>=upper_time_limit:
        print "\n\nStart day cannot be the end day or a day after"
    starting_day = str(raw_input("Select between " + str(lower_time_limit) + " and " + str(upper_time_limit) + ": "))
    try:
        int(starting_day)
    except:
        start_day = -1
        continue
    start_day = int(starting_day)
        

print "\n\nPlease select the end day for ETa Calculation"
ending_day = str(raw_input("Select between " + str(start_day) + " and " + str(upper_time_limit) + ": "))
try:
    end_day = int(ending_day)
except:
    end_day = -1

while not(start_day<end_day<=upper_time_limit):                                                                     #checking that the selection is within range of available data
    if end_day == -1:
        print "Enter a number"
    elif end_day>start_day:
        print ("\n\nThe selected day is beyond the limit of " + str(upper_time_limit))       
    elif end_day<=start_day:
        print "\n\nEnd day cannot be the start day or a day before"
    ending_day = str(raw_input("Select between " + str(start_day) + " and " + str(upper_time_limit) + ": "))
    try:
        int(ending_day)
    except:
        end_day = -1
        continue
    end_day = int(ending_day)

print ("selected period is from day " + str(start_day) + " to day " + str(end_day))

#%%
#get floor through ceiling rasters
start_index = 1000*100
beginning_raster=copy.copy(start_day)
while start_index == 1000*100:
    try:
        if beginning_raster<100:
            start_index = Kc_files.index(input_folder+ "\\" + "0"+str(beginning_raster)+".tif")
        elif beginning_raster>=100:
            start_index = Kc_files.index(input_folder+ "\\" + str(beginning_raster)+".tif")
    except:
        beginning_raster = beginning_raster-1        
        
end_index = 1000*100
last_raster=copy.copy(end_day)
while end_index == 1000*100:
    try:
        if last_raster<100:
            end_index = Kc_files.index(input_folder+ "\\" + "0"+str(last_raster)+".tif")
        elif last_raster>=100:
            end_index = Kc_files.index(input_folder+ "\\" + str(last_raster)+".tif")
    except:
        last_raster = last_raster+1 
 
#%%
#get geoinfo
driver, NDV, xsize, ysize, GeoT, Projection = qgis.GetGeoInfo(Kc_files[start_index], Subdataset = 0)
# make dictionary with raster of each kc
Ras_dictionary={}
print "Importing and interpolating for missing values in the raster files"
for i in range(start_index,end_index+1):
    d_element = 'timestep_' + str(int(Kc_files[i][-7:-4])) #dictionary element
    tif2np = qgis.OpenAsArray(Kc_files[i], Bandnumber = 1, dtype = 'float32', nan_values = False)   # nan_values = False
    Ras_dictionary[d_element] = tif2np
    
#%% Selecting method of interpolation
interpolation_methods = str(raw_input("Select interpolation method ;  1 = Fixed, 2 = Linear:  "))
try:
    interpolation_method = int(interpolation_methods)
except:
    interpolation_method = int(3)
        
while not ((interpolation_method == int(1) or interpolation_method == int(2))):
    print "Choose either 1 = Fixed or 2 = Linear"
    interpolation_methods = str(raw_input("Select interpolation method ;  1 = Fixed, 2 = Linear:  "))
    try:
        interpolation_method = int(interpolation_methods)
    except:
        interpolation_method = int(3)

#%%        
# Fixed interpolation of kc_values between timesteps
if interpolation_method == int(1):
    fnx.Do_Fixed_Inerpolation(start_day,end_day,start_index,end_index,Kc_files,Ras_dictionary)

#%%
if interpolation_method == int(2):
    fnx.Do_Linear_Interpolation(start_index,end_index,start_day,end_day,Ras_dictionary,Kc_files)
#%%
ETa={}
ETa['seasonalETa']=np.empty_like(Ras_dictionary['timestep_' + str(start_day)])

if start_day!=beginning_raster:                                 #keep base_raster if user includes it in time range
    del Ras_dictionary['timestep_' + str(beginning_raster)]
if end_day!=last_raster:
    del Ras_dictionary['timestep_' + str(last_raster)]          #keep ceilling_raster if user includes it in time range

#%%
#ETa calculation

print '\nGetting ET0 data and Multiplying with Kc Maps...'
for i in ET0_array:
    dict_key = 'timestep_' + str(i.split(',')[0])
    if dict_key in Ras_dictionary:
        ETa[dict_key] = Ras_dictionary[dict_key] * float(i.split(',')[1])       #get daily ETa
        if int(i.split(',')[0])>start_day:
            ETa['seasonalETa'] = ETa['seasonalETa'] + ETa[dict_key]
        else:
            ETa['seasonalETa'] = ETa[dict_key]                                          #get seasonal ETa
#%%
ETa_out = str(raw_input('Name the folder you want to save your ETa results: '))

daily_ETA_subdir = '/' + ETa_out + '/Daily ETa' + '/'
seasonal_ETA_subdir = '/' + ETa_out + '/Seasonal ETa' + '/'

print '\nCreating ETa Tif images...'
for key in ETa:   
    if key=='seasonalETa':
        qgis.CreateGeoTiff(key, ETa[key], driver, NDV, xsize, ysize, GeoT, Projection,output_folder,seasonal_ETA_subdir)        #saving seasonal tiff
    else:
        qgis.CreateGeoTiff(key, ETa[key], driver, NDV, xsize, ysize, GeoT, Projection,output_folder,daily_ETA_subdir)        #saving daily tiff
print '\nDone...'


#%%
print'Biomass interpolation?'
bio_permission =fnx.yes_or_no()
if bio_permission==int(1):
    print "Choose folder containing the 'Wp' files"
    wp_folder = fnx.SelectFolderDialog()
    wp_files = fnx.ListFilesInFolder (wp_folder,'.tif')                                  #get list of wp files
    
    while len(wp_files) == 0:
        print "The selected folder does not have 'Wp' files \nPlease choose a folder with 'Wp' files"
        wp_folder = fnx.SelectFolderDialog()
        wp_files = fnx.ListFilesInFolder (wp_folder,'.tif')

    print'Reading wp files'
    RasWP_dictionary={}
    for i in range(start_index,end_index+1):
        d_element = 'timestep_' + str(int(wp_files[i][-7:-4])) #dictionary element
        tif2np = qgis.OpenAsArray(wp_files[i], Bandnumber = 1, dtype = 'float32', nan_values = False)   # nan_values = False
        RasWP_dictionary[d_element] = tif2np
    
    print'Calculating seasonal average wp of available rasters in selected period'
    average_wp = fnx.Average_WP(RasWP_dictionary,start_index,end_index,wp_files)
    
    Biomass_subdir = '/' + ETa_out + '/Seasonal Biomass/'
    qgis.CreateGeoTiff('Seasonal_Average_WP',average_wp, driver, NDV, xsize, ysize, GeoT, Projection,output_folder,Biomass_subdir)        #saving seasonal tiff

    biomass_total1 = np.multiply(np.array(ETa['seasonalETa']),average_wp)
    biomass_total = np.multiply(biomass_total1,10)
    qgis.CreateGeoTiff('Seasonal_Total_biomass',biomass_total, driver, NDV, xsize, ysize, GeoT, Projection,output_folder,Biomass_subdir)        #saving seasonal tiff
    
#%%
print'\n\n\nMake a plot of ETa over time for selected pixel?'
graph_permission = fnx.yes_or_no()
if graph_permission==int(1):
    fnx.Plot_ETa(xsize,ysize, start_day,end_day,ETa,output_folder, ETa_out)
#%% Open Results Folder
resultfolder = output_folder + '\\' + ETa_out
os.startfile(resultfolder)
#%%
