# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""

from .oimModel import *
from .oimModel import _standardParameters
from .oimData import *
from .oimFitter import *
from .oimSimulator import *
from .oimUtils import *
from .oimDataFilter import *

from os.path import join, dirname, split
import inspect

__pkg_dir__ = dirname(inspect.getfile(inspect.currentframe()))


if split(__pkg_dir__)[-1] == "":
    __git_dir__ = dirname(split(__pkg_dir__)[0])
else:
    __git_dir__ = split(__pkg_dir__)[0]

