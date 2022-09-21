# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""

import oimodeler as oim
import os

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

data=oim.oimData(files)

print(data)
