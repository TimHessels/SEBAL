# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
         IHE Delft 2017
Contact: t.hessels@un-ihe.org
Repository: https://github.com/wateraccounting/SEBAL
Module: SEBAL/SEBAL


Description:
This module contains a compilation of scripts and functions to run pySEBAL
(http://www.wateraccounting.org/)
"""


from .EESEBAL_code import main as EESEBAL

__all__ = ['EESEBAL']

__version__ = '0.1'
