# -*- coding: utf-8 -*-
"""
Created on Tue May 03 13:12:18 2016

@author: tih
"""

import SEBAL

inputExcel =  r"D:\Project_WAKAS\Input_Data\SEBAL\Excel_SEBAL_v3_4_4.xlsx"

for number in range(3,13):
    SEBAL.pySEBAL.pySEBAL_code.main(number,inputExcel)


