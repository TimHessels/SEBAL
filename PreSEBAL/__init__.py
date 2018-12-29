# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
         IHE Delft 2017
Contact: t.hessels@un-ihe.org
Repository: https://github.com/TimHessels/SEBAL
Module: SEBAL/preSEBAL


Description:
This module contains a compilation of scripts and functions to run pySEBAL
"""


from .preSEBAL_0 import main
from .preSEBAL_1 import main
from .preSEBAL_1b import main
from .preSEBAL_2 import main
from .preSEBAL_2b import main
from .pyKARIMET import main
#from .preSEBAL_3 import main
#from .preSEBAL_4 import main
#from .preSEBAL_4b import main

__all__ = ['preSEBAL_0', 'preSEBAL_1', 'preSEBAL_1b', 'preSEBAL_2', 'preSEBAL_2b', 'pyKARIMET']

__version__ = '0.1'
