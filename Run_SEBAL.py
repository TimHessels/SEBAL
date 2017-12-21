# -*- coding: utf-8 -*-
"""
Created on Tue May 03 13:12:18 2016

@author: tih
"""

import SEBAL

inputExcel = r"$HOME\SEBAL_Codes\InputEXCEL_v3_3_7.xlsx"

for number in range(2,4):
    try:
        SEBAL.SEBALcode(number,inputExcel)
    except:
        print 'SEBAL did not run line %d fully' % number
        

