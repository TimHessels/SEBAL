# -*- coding: utf-8 -*-
"""
Created on Tue May 03 13:12:18 2016

@author: tih
"""

import SEBAL

inputExcel = r"$HOME\InputEXCEL_v3_3_6.xlsx"

for number in range(2,3):
    try:
        SEBAL.SEBALcode(number,inputExcel)
    except:
        print 'SEBAL did not run line %d fully' % number
        

