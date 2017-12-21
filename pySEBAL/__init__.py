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


from .pySEBAL_code import main as pySEBAL

__all__ = ['pySEBAL']

__version__ = '0.1'
