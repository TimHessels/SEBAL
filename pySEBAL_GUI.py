# -*- coding: utf-8 -*-
"""
Created on Thu Jan 04 15:04:39 2018

@author: tih
"""

from PyQt5 import QtCore, QtWidgets
import SEBAL.GUI.pySEBAL_pyGUI as pySEBAL_pyGUI
import SEBAL
import time
import os
import glob
import re
import datetime as dt
from osgeo import gdal, gdal_array
import osr
import numpy as np
import pandas as pd

class Interface_pySEBAL (QtWidgets.QTabWidget, pySEBAL_pyGUI.Ui_pySEBALGUI):

    def __init__(self):
        super(self.__class__,self).__init__()
        self.setupUi(self)
        
        # Fixed Size of Application
        self.setFixedSize(self.width(),self.height())
        
# =============================================================================
#         #Signals and Slots
# =============================================================================
        '''
        # Actual Date
        now_time=time.strftime("%Y%m%d")
            
        # SEBAL
        self.dateEdit_from.setDate(QtCore.QDate(int(now_time[0:4]),int(now_time[4:6]),int(now_time[6:8])))
        self.dateEdit_to.setDate(QtCore.QDate(int(now_time[0:4]),int(now_time[4:6]),int(now_time[6:8])))
        
        # =============================================================================
        #         SEBAL
        # =============================================================================
        
        self.dateEdit_from.dateChanged.connect(self.setdate_sebal)
        self.dateEdit_to.dateChanged.connect(self.setdate_sebal)
        self.checkBox_SebalLot.stateChanged.connect(self.verify)
        self.radioButton_Landsat.clicked.connect(self.Hot_Cold_pixels)
        self.radioButton_VIIRS_PROBAV_100.clicked.connect(self.Hot_Cold_pixels)
        self.radioButton_VIIRS_375.clicked.connect(self.Hot_Cold_pixels)
        self.BrowseButton_DEMSebal.clicked.connect(self.Set_folders)
        self.BrowseButton_InSebal.clicked.connect(self.Set_folders)
        self.BrowseButton_OutSebal.clicked.connect(self.Set_folders)
        self.BrowseButton_MeteoSebal.clicked.connect(self.Set_folders)
        self.BrowseButton_SoilSebal.clicked.connect(self.Set_folders)
        self.buttonBox_OKCancel.accepted.connect(self.ApplySEBAL)
        self.ButtonBox.accepted.connect(self.Apply_Test)
        '''
        
        self.pushButton.clicked.connect(self.browse_folder)

        
# =============================================================================
# #Functions or Methodes call
# =============================================================================
  
    def reject(self):
        self.pySEBAL_pyGUI.close
        
    def browse_folder(self):
            sender = self.sender()
            objectName = sender.objectName()
            
            if objectName == 'pushButton':
                a=self.lineEdit.text()
                filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Select Excel File', 'C:\\', 'Excel Files (*.xlsx)')
                if filename=='':
                    self.lineEdit.setText(a)
                else:
                    self.lineEdit.setText(filename[0])
    
    def accept(self):
        self.ExcelFile_sebal=str(self.lineEdit.text())
        self.Startvalue=int(self.spinBox.value())
        self.Endvalue=int(self.spinBox_2.value())
                
        if (self.ExcelFile_sebal==''):
            QtWidgets.QMessageBox.warning(self, "Warning",
                                       "Select the SEBAL excel file",
                                       QtWidgets.QMessageBox.Ok)
        else:
            self.progdialog = QtWidgets.QProgressDialog( "SEBAL in Progress. Please Wait...", "Cancel", 0 , 100, self)
            self.progdialog.setWindowTitle("SEBAL Progression")
            self.progdialog.setWindowModality(QtCore.Qt.ApplicationModal)
            self.progdialog.setFixedSize(self.progdialog.size())
            self.progdialog.show()  
            
            self.TT_sebal = Sebal_Thread(self)
            self.TT_sebal.updateProgress.connect(self.Sebal_Progressed)
            self.TT_sebal.exception.connect(self.Sebal_exception)
            self.TT_sebal.start()
            
    def Sebal_Progressed(self,value_pg):
        self.progdialog.setLabelText("SEBAL Accomplishment Percentage")
        self.progdialog.setValue(value_pg)
        
    def Sebal_exception(self, message):
        QtWidgets.QMessageBox.critical(self, "Error", message, QtWidgets.QMessageBox.Ok)
        
# =============================================================================
#         #Threads and Dialog Class
# =============================================================================
    
class Sebal_Thread(QtCore.QThread):

    updateProgress = QtCore.pyqtSignal(int)
    exception = QtCore.pyqtSignal(str)
    
    def __init__(self,Interface_pySEBAL):
        QtCore.QThread.__init__(self)
        
        self.Interface=Interface_pySEBAL
        
    def run(self):

        ExcelFile_sebal=self.Interface.ExcelFile_sebal
        Startvalue=int(self.Interface.Startvalue)
        Endvalue=int(self.Interface.Endvalue)
        i = 0
        
        for number in range(Startvalue,Endvalue+1):
            
            try:
                SEBAL.pySEBAL.pySEBAL_code.main(number,ExcelFile_sebal)
            except:
                message="SEBAL Model Processing on Landsat Images are not well executed"
                self.exception.emit(message)
            self.updateProgress.emit(int(round((i+1)*(100./len(range(Startvalue,Endvalue+1))))))
            i+=1

if __name__== '__main__':

    import sys
    app=QtWidgets.QApplication(sys.argv)
    form=Interface_pySEBAL()
    form.show()
    app.exec_()

