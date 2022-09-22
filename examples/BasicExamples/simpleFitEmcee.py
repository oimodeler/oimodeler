# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:59:15 2022

@author: Ame
"""

import oimodeler as oim
import os

path = os.path.dirname(oim.__file__)

# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

#Building a oimodeler model with the same parameters
ud=oim.oimUD(d=3,f=0.5)
pt=oim.oimPt(f=1)
model=oim.oimModel([ud,pt])

#%%

#setting limits of the parameter space but also setting x,y of the UD as a 
#free parameters and f of the pt as fixed (to 1)
ud.params['d'].set(min=0.01,max=20)
ud.params['x'].set(min=-50,max=50,free=True)
ud.params['y'].set(min=-50,max=50,free=True)
ud.params['f'].set(min=0.,max=10.)
pt.params['f'].free=False


#%%
#Create a new fitter with 32 walkers and the list of oifits files and the model
fit=oim.oimFitterEmcee(files,model,nwalkers=32)
#print(fit._logProbability([0,10,1,5]))


# Prepare the fitter.Here we ste the intial positions of all walkers to 
#the current parameters values of our model. 
fit.prepare(init="random")

#Printing the initial values of the walkers
print("Initial values of the free parameters for the {} walkers".format(fit.params["nwalkers"].value))
print(fit.initialParams)

#run a 1000 steps fit with fixed starting inital and 1000 steps 
fit.run(nsteps=10000,progress=True)


#%% Getting results from the mcmc run
# Result method return the models values and 1 sigma uncertainties 
#(low, up and avg) computed from the distrib qnuatiles (0.16,0.84)
#It also compute simulated Data with the best/mean/median model

# Results can be computed with various options : 
    
#Best return the  model with the lowest chi2
best,err_l,err_u,err=fit.getResults(mode='best')

#Mean return the mean value of the parameters excluding first "discard" steps
# and models with a chi2 >chi2min*chi2limifact 
mean,err_l,err_u,err=fit.getResults(mode='mean',discard=3000,chi2limfact=20)

#Median return the median model
median,err_l,err_u,err=fit.getResults(mode='median',discard=2000,chi2limfact=20)

#%% plotting corner and walkers plots and saving them
figCorner,axeCorner=fit.cornerPlot(discard=2000,savefig=
                os.path.join(path,os.pardir,"images","SimpleFitCorner.png"))

figWalkers,axeWalkers=fit.walkersPlot(cmap="plasma_r",savefig=
                os.path.join(path,os.pardir,"images","SimpleFitWalkers.png"))


#%%
#Plotting the median model (last computed) over the data
fit.simulator.compute(computeChi2=True,computeSimulatedData=True)
figSim,axSim=fit.simulator.plot(["VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI"])

median,err_l,err_u,err=fit.getResults(mode='median',discard=1000,chi2limfact=20)

