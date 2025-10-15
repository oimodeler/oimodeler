# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""
import inspect
from os.path import split
from pathlib import Path

import matplotlib.projections as proj
import numpy as np

from .oimBasicFourierComponents import *
from .oimComponent import *
from .oimCustomComponents import *
from .oimData import *
from .oimDataFilter import *
from .oimFTBackends import *
from .oimFitter import *
from .oimFluxData import *
from .oimModel import *
from .oimOptions import oimOptions
from .oimParam import *
from .oimParam import _interpolators, _standardParameters
from .oimPlots import *
from .oimSimulator import *
from .oimUtils import *
from .oimUtils import _oimDataType, _oimDataTypeArr, _oimDataTypeErr

np.seterr(invalid="ignore")
proj.register_projection(oimAxes)

__version__ = "0.9.0"
__pkg_dir__ = Path(inspect.getfile(inspect.currentframe())).parent

if split(__pkg_dir__)[-1] == "":
    __git_dir__ = str(Path(split(__pkg_dir__)[0]).parent)
else:
    __git_dir__ = str(split(__pkg_dir__)[0])
__pkg_dir__ = str(__pkg_dir__)
