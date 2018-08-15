'''
This Module contains the functions used in the Main program
'''

import os
import matplotlib.pyplot as plt
import numpy as np
import copy
from tkFileDialog import askdirectory, askopenfilename
from Tkinter import Tk
#%%
#Functions

def OpenFileDialog():
    Tk().withdraw()
    file_name = askopenfilename()
    return file_name

def SelectFolderDialog():
    Tk().withdraw()
    folder = askdirectory()
    return folder
    
def yes_or_no():
    permissions = str(raw_input("Yes = 1, No = 0 : "))
    try:
        permission = int(permissions)
    except:
        permission = 3
        print "Chose Yes or No"
        
    while not((permission==int(1)) or (permission==int(0))):
        permissions = str(raw_input("Yes = 1, No = 0 : "))
        try:
            permission = int(permissions)
        except:
            permission = 3
            
    return permission

def ListFilesInFolder(folder,extension):
    list_of_files = [os.path.join(folder,fn) for fn in next(os.walk(folder))[2] if fn[-4:] == extension]
    return list_of_files

def import_txt (filename):
    f=open(filename, 'r')
    text=f.readlines()
    f.close()
    return text

#%%   
def Plot_ETa(xsize,ysize, start_day,end_day,ETa_Dictionary,output_folder, ETa_out):
    while True:
        print ("\n\n\nSelect the coordinates (X,Y) for the point you want to plot.\nX from 1 to " + str(xsize) + " and Y from 1 to " + str(ysize))
        coordinates = str(raw_input("X,Y (e.g, 4,6) : "))
        try:
            coordinates.split(",")[1]
        except:
            coordinates = '0,0'
            print 'Please use correct format'
            
        while not(( 0 < int(coordinates.split(",")[1]) <= ysize ) and (0 < int(coordinates.split(",")[0]) <= xsize)): #checking correct format and ranges
            print "Selected pixel does not exist on grid, check the limits" + "\nX from 1 to " + str(xsize) + " and Y from 1 to " + str(ysize)
            coordinates = str(raw_input("X,Y (e.g, 4,6) : "))
            try:
                coordinates.split(",")[1]
            except:
                coordinates = '0,0'
                print 'Please use correct format'
                continue
            
        time_axis=[]
        ETa_value=[]
        for d in range (start_day,end_day+1):
            time_axis.append(d)
            ETa_value.append(ETa_Dictionary[('timestep_' + str(d))][int(coordinates.split(",")[1])-1][int(coordinates.split(",")[0])-1])
            
            
        if not os.path.isdir(output_folder + '/' + ETa_out + "/" + "Figures"):
            os.makedirs(output_folder + '/' + ETa_out + "/" + "Figures")
        
        plt.clf()   # Clear figure
        plt.plot(time_axis,ETa_value)
        plt.xlabel('Time (day)')
        plt.ylabel('ETa (mm)')
        plt.title('Daily ETa for pixel (' + coordinates + ')')
        plt.xlim(start_day,end_day)
        plt.grid(True)
        plt.savefig(output_folder + '/' + ETa_out + "/" + "Figures" "/" + 'ETa_timeseries_for_' + coordinates + '.jpg')
        
        print 'Make another plot?'
        answer=yes_or_no()
        if answer == 1:
            continue
        else:
            break
#%%
def Do_Linear_Interpolation(start_index,end_index,start_day,end_day,Raster_Dictionary,List_of_files):
    print '\nUsing Linear Interpolation...'
    
    lower_step = 1000**1000
    #lower_raster = ''
    print '\nInterpolating Kc Maps...'
    for k in range (start_index,end_index+1):
        upper_step = int(List_of_files[k][-7:-4])
        
        if upper_step < lower_step:
            pass
        else:
            day_count = upper_step -lower_step
            weights = np.linspace(0, 1, day_count+1)
            base_raster = 'timestep_' + str(lower_step) #used inside if statement in following for loop
            ceilling_raster = 'timestep_' + str(upper_step) #used inside if statement in following for loop
            for m in range (1,day_count):
                interpolated_raster = 'timestep_' + str(lower_step + m)
                x = weights[m]
                if (start_day <= lower_step + m <= end_day):
                    Raster_Dictionary[interpolated_raster] = (1-x)*(Raster_Dictionary[base_raster])  +  x*(Raster_Dictionary[ceilling_raster])
                
        lower_step = copy.copy(upper_step)
        #lower_raster = List_of_files[k]
#%%
def Do_Fixed_Inerpolation(start_day,end_day,start_index,end_index,List_of_Files,Raster_Dictionary):
        #block type
    lower_step = 1000**1000
    #lower_raster = ''
    
    for k in range (start_index,end_index+1):
        upper_step = int(List_of_Files[k][-7:-4])
        
        if upper_step < lower_step:
            pass
        else:
            day_count = upper_step -lower_step
            base_raster = 'timestep_' + str(lower_step) #used inside if statement in following for loop
            ceilling_raster = 'timestep_' + str(upper_step) #used inside if statement in following for loop
            for j in range (1,day_count):
                interpolated_raster = 'timestep_' + str(lower_step + j)
                
                if not ((int(lower_step + j) <= end_day) and (int(lower_step + j)>=start_day)):
                    continue
                if j < (int(day_count/2)):
                    Raster_Dictionary[interpolated_raster] = copy.deepcopy(Raster_Dictionary[base_raster])
                elif j >=  (int(day_count/2)):
                    Raster_Dictionary[interpolated_raster] = copy.deepcopy(Raster_Dictionary[ceilling_raster])
                
        lower_step = copy.copy(upper_step)
        #lower_raster = List_of_Files[k]

#%%
def Average_WP(Raster_Dictionary,start_index,end_index,List_of_Files):
    for k in range(start_index,end_index+1):
        if k==start_index:
            name='timestep_' + str(int(List_of_Files[k][-7:-4]))
            raster=np.array(Raster_Dictionary[name])
            sum_wp=raster
        else:
            name='timestep_' + str(int(List_of_Files[k][-7:-4]))
            raster=np.array(Raster_Dictionary[name])
            stack_wp=np.stack((sum_wp,raster))
            sum_wp=np.sum(stack_wp,axis=0)
    mean_wp=np.divide(sum_wp,(end_index-start_index))
    
    return mean_wp
    