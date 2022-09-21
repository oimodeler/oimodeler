# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:59:15 2022

@author: Ame
"""

import oimodeler as oim
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import os
from datetime import datetime
from astropy.io import fits
import corner.corner

path = os.path.dirname(oim.__file__)

# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

#Building a oimodeler model with the same parameters
ud=oim.oimUD(d=3,f=0.5)
pt=oim.oimPt(f=1)
model=oim.oimModel([ud,pt])

#setting limits of the parameter space but also setting x,y of the UD as a 
#free parameters and f of the pt as fixed (to 1)
ud.params['d'].set(min=0.01,max=20)
ud.params['x'].set(min=-50,max=50,free=True)
ud.params['y'].set(min=-50,max=50,free=True)
ud.params['f'].set(min=0.,max=10.)
pt.params['f'].free=False

#Create a new fitter wwith 32 walkers and the list of oifits files and the model
fit=oim.oimFitterEmcee(files,model,nwalkers=32)
#print(fit._logProbability([0,10,1,5]))


# Prepare the fitter.Here we ste the intial positions of all walkers to 
#the current parameters values of our model. 
fit.prepare(init="random")

#Printing the initial values of the walkers
print(fit.initialParams)

#run a 1000 steps fit with fixed starting inital and 1000 steps 
fit.run(nsteps=1000,progress=True)

#%%

#Saving corner and walkers plots
filename_corner=os.path.join(path,os.pardir,"images","SimpleFitCorner.png")
filename_steps=os.path.join(path,os.pardir,"images","SimpleFitWalkers.png")
fit.cornerPlot(discard=400,savefig=filename_corner)
fit.walkersPlot(cmap="plasma_r",savefig=filename_steps)
