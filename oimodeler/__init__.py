# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""
import inspect
from pathlib import Path
from os.path import split

import matplotlib.projections as proj
import numpy as np

from .oimBasicFourierComponents import *
from .oimCustomComponents import *
from .oimComponent import *
from .oimData import oimDataGetWl, oimDataType, oimGetDataValErrAndTypeFlag, oimData, \
    oimDataCheckData, oimDataGetVectCoord, loadOifitsData
from .oimFluxData import *
from .oimDataFilter import *
from .oimFTBackends import numpyFFTBackend, FFTWBackend
from .oimFitter import *
from .oimModel import *
from .oimOptions import oimOptions
from .oimParam import *
from .oimParam import _standardParameters, _interpolators
from .oimPlots import *
from .oimSimulator import *
from .oimUtils import *
from .oimUtils import  _oimDataType, _oimDataTypeArr, _oimDataTypeErr

np.seterr(invalid='ignore')
proj.register_projection(oimAxes)

__version__ = "0.8.1"

# TODO: After pathlib change of all `oimodeler` modules, remove str here
__pkg_dir__ = Path(inspect.getfile(inspect.currentframe())).parent

if split(__pkg_dir__)[-1] == "":
    __git_dir__ = str(Path(split(__pkg_dir__)[0]).parent)
else:
    __git_dir__ = str(split(__pkg_dir__)[0])
__pkg_dir__ = str(__pkg_dir__)
