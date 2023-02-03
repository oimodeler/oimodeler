# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""
import inspect
from os.path import join, dirname, split

import matplotlib.projections as proj
import numpy as np

from .oimBasicFourierComponents import *
from .oimCustomComponents import *
from .oimData import oimData
from .oimFitter import oimFitterEmcee
from .oimModel import oimModel
from .oimPlots import oimAxes, oimPlot
from .oimSimulator import oimSimulator


np.seterr(invalid='ignore')

proj.register_projection(oimAxes)

__pkg_dir__ = dirname(inspect.getfile(inspect.currentframe()))


if split(__pkg_dir__)[-1] == "":
    __git_dir__ = dirname(split(__pkg_dir__)[0])
else:
    __git_dir__ = split(__pkg_dir__)[0]
