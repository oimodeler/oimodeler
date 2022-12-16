# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""
#from oimOptions import *
from oimOptions import oimOptions
from oimFTBackends import *
from oimParam import *
from oimParam import _standardParameters
from oimModel import *
from oimComponent import *
from oimBasicFourierComponents import *
import oimCustomComponents as oimCustomComponents
from oimData import *
from oimData import _oimDataType, _oimDataTypeErr, _oimDataTypeArr
from oimFitter import *
from oimSimulator import *
from oimUtils import *
from oimDataFilter import *
from oimPlots import *


import numpy as np
np.seterr(invalid='ignore')

from os.path import join, dirname, split
import inspect
import matplotlib.projections as proj

proj.register_projection(oimAxes)

__pkg_dir__ = dirname(inspect.getfile(inspect.currentframe()))


if split(__pkg_dir__)[-1] == "":
    __git_dir__ = dirname(split(__pkg_dir__)[0])
else:
    __git_dir__ = split(__pkg_dir__)[0]

