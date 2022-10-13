# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:59:15 2022

@author: Ame
"""

import oimodeler as oim

#%%
import os
path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData)]
#%%
ud=oim.oimUD(d=3,f=0.5,x=5,y=-5)
pt=oim.oimPt(f=1)

#%%
print(ud)
print(ud.params['d'])
print(ud.params['x'])

#%%
ud.params['d'].set(min=0.01,max=20)
ud.params['x'].set(min=-50,max=50,free=True)
ud.params['y'].set(min=-50,max=50,free=True)
ud.params['f'].set(min=0.,max=10.)
pt.params['f'].free=False

#%%
model=oim.oimModel([ud,pt])

#%%

print(model.getParameters())
print(model.getFreeParameters())

#%%
sim=oim.oimSimulator(data=files,model=model)
sim.compute(computeChi2=True,computeSimulatedData=True)

#%%
print("Chi2r = {}".format(sim.chi2r))


fig0,ax0= sim.plot(["VIS2DATA","T3PHI"],
    savefig=os.path.join(path,os.pardir,"images","gettingStarted_model0.png"))

#%%

fit=oim.oimFitterEmcee(files,model,nwalkers=10)
fit.prepare(init="random")
print(fit.initialParams)
fit.run(nsteps=2000,progress=True)

#%%

figWalkers,axeWalkers=fit.walkersPlot(
    savefig=os.path.join(path,os.pardir,"images","gettingStarted_Walkers.png"))

figCorner,axeCorner=fit.cornerPlot(discard=1000, 
     savefig=os.path.join(path,os.pardir,"images","gettingStarted_corner.png"))

#%%
median,err_l,err_u,err=fit.getResults(mode='median',discard=1000)

figSim,axSim=fit.simulator.plot(["VIS2DATA","T3PHI"],
    savefig=os.path.join(path,os.pardir,"images","gettingStarted_modelFinal.png"))
print("Chi2r = {}".format(fit.simulator.chi2r))

figImg,axImg,im=model.showModel(512,0.1,normPow=0.1,figsize=(6,4.8),
    savefig=os.path.join(path,os.pardir,"images","gettingStarted_modelImage.png"))

