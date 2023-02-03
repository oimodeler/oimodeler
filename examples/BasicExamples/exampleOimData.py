# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
import os

import oimodeler as oim

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData)]

data=oim.oimData(files)

print(data.data)

data.prepareData()
print(data.vect_u)
print(data.vect_v)   
print(data.vect_wl)  
print(data.vect_u.shape)
