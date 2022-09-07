# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""

import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import astropy.units as u
import os

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

data=oim.OImData(files)

print(data)
data.vectorizeData()
