# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
         IHE Delft 2017
Contact: t.hessels@un-ihe.org
Repository: https://github.com/wateraccounting/SEBAL
Module: SEBAL/preSEBAL


Description:
This module contains a compilation of scripts and functions to run pySEBAL
(http://www.wateraccounting.org/)
"""


from .preSEBAL import main as preSEBAL
from .preSEBAL_light import main as preSEBAL_light

__all__ = ['preSEBAL', 'preSEBAL_light']

__version__ = '0.1'
