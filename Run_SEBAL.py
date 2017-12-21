# -*- coding: utf-8 -*-
"""
Created on Tue May 03 13:12:18 2016

@author: tih
"""

import SEBAL

inputExcel = r"$HOME\SEBAL\Excel_SEBAL_v3_3_8.xlsx"

for number in range(2,4):
    SEBAL.pySEBAL.pySEBAL_code.main(number,inputExcel)

        

