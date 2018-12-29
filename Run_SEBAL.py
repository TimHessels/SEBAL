# -*- coding: utf-8 -*-
"""
Created on Tue May 03 13:12:18 2016

@author: tih
"""

import SEBAL

inputExcel =  r"G:\SEBAL_Tadla\Excel_SEBAL_v3_4_0.xlsx"

for number in range(2,3):
    SEBAL.pySEBAL.pySEBAL_code_new.main(number,inputExcel)



