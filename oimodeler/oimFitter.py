# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 14:32:59 2022

@author: Ame
"""

import numpy as np
from oimodeler.oimModel import oimParam
import oimodeler as oim
import emcee
import matplotlib.pyplot as plt
import corner

class oimFitter(object):
    params={}
    def __init__(self,*args,**kwargs):
        nargs=len(args)
        if nargs==2:
            self.simulator=oim.oimSimulator(args[0],args[1])
        elif nargs==1:
            self.simulator=args[0]
        else:
            raise TypeError("Wrong number of arguments")
            
        self.data=self.simulator.data
        self.model=self.simulator.model 
        
        self.isPrepared=False
        
        self._eval(**kwargs)

    def _eval(self,**kwargs):
        
        for key in self.params.keys():
            if key in kwargs.keys(): 
                self.params[key].value=kwargs.pop(key)
        return kwargs
                    
                    
    
    def prepare(self,**kwargs):
        self.freeParams=self.model.getFreeParameters()
        self.nfree=len(self.freeParams)
        
        self.limits={}
        for key in self.freeParams:
            self.limits[key] = (self.freeParams[key].min, 
                                self.freeParams[key].max)
            
        kwargs=self._eval(**kwargs)
        
        self.simulator.prepareData()
        kwargs=self._prepare(**kwargs)  
        self.isPrepared=True
        return kwargs
        
    
    def run(self,**kwargs):
        if self.isPrepared==False:
            raise TypeError("Fitter not initialized")
        self._run(**kwargs)
        return kwargs
        
    def getResults(**kwargs):
        return kwargs
    
    def _prepare(self,**kwargs):
        return kwargs
        
    def _run(self,**kwargs):
        return kwargs    



###############################################################################
  
class oimFitterEmcee(oimFitter):
    
    def __init__(self,*args,**kwargs):
        
        self.params["nwalkers"]=oimParam(name="nwalkers",value=16,mini=1,
                    description="Number of walkers")  

        
        super().__init__(*args,**kwargs)
        

    def _prepare(self,**kwargs):
        
        if not('init' in kwargs):
            init='random'
        else:
            init=kwargs.pop('init')
        if init=="random":
            self.initialParams=self._initRandom()
        elif init=="fixed":
            self.initialParams=self._initFixed()
            
        print(self.params["nwalkers"].value)
        self.sampler = emcee.EnsembleSampler(self.params["nwalkers"].value, 
                        self.nfree,self._logProbability,
                        moves=[(emcee.moves.DEMove(), 0.8), 
                               (emcee.moves.DESnookerMove(), 0.2)],**kwargs)
        return kwargs   
            
    def _initRandom(self):
        
        nw=self.params['nwalkers'].value
        initialParams=np.ndarray([nw,self.nfree])
        
        for iparam,parami in enumerate(self.freeParams.values()):
            initialParams[:,iparam]=np.random.random(
                self.params['nwalkers'].value)*(parami.max-parami.min)+parami.min
        
        return initialParams
              
    def _initFixed(self):
        
        nw=self.params['nwalkers'].value
        
        initialParams=np.ndarray([nw,self.nfree])
        
        for iparam,parami in enumerate(self.freeParams.values()):
            initialParams[:,iparam]=np.ones(nw)*parami.value
                    
        return initialParams
        
        
        
    def _run(self,**kwargs):
        self.sampler.run_mcmc(self.initialParams,**kwargs )
        return kwargs  
    
    
    def _logProbability(self, theta):

        for iparam,parami in enumerate(self.freeParams.values()):
            parami.value=theta[iparam]

        
        for i,key in enumerate(self.freeParams):
            val = theta[i]
            low,up = self.limits[key]
            if not (low < val < up):
                return -np.inf
        
        self.simulator.compute(computeChi2=True)
        return -0.5 * self.simulator.chi2r  
        
    
    def cornerPlot(self,discard=0,chi2lim=20,savefig=None,**kwargs):
        pnames=list(self.freeParams.keys())
        punits=[p.unit for p in list(self.freeParams.values())]
        
        labels=[]
        for namei,uniti in zip(pnames,punits):
            txt=namei
            if uniti.to_string()!="":
                txt+=" ("+uniti.to_string()+")"  
            labels.append(txt)
        
        c=self.sampler.get_chain(discard=discard,flat=True)
        chi2=-2*self.sampler.get_log_prob(discard=discard,flat=True)

        idx=np.where(chi2<chi2lim)[0]
        c2=c[idx,:]

        corner.corner(c2,labels=labels,quantiles=[0.16, 0.5, 0.84],show_titles=True,bins=50
                      ,smooth=2,smooth1d=2,fontsize=8,title_kwargs={'fontsize':8},use_math_text=True)

        if savefig!=None:
            plt.savefig(savefig)

    def walkersPlot(self,savefig=None,**kwargs):
        fig, ax = plt.subplots(4, figsize=(10, 7), sharex=True)
        samples=self.sampler.get_chain()
        chi2=-2*self.sampler.get_log_prob()
        pnames=list(self.freeParams.keys())
        punits=[p.unit for p in list(self.freeParams.values())]
        x=np.arange(chi2.shape[0])

        samples=self.sampler.get_chain()


        xm=np.outer(x,np.ones(chi2.shape[1]))
        chi2f=chi2.flatten()
        xf=xm.flatten()
        idx=np.argsort(-1*chi2f)

        chi2f=chi2f[idx]
        xf=xf[idx]
        samples=samples.reshape([samples.shape[0]*samples.shape[1],samples.shape[2]])[idx,:]
        for i in range(self.nfree):
            
            scale=ax[i].scatter(xf,samples[:,i],c=chi2f,marker=".",
                           vmin=chi2.min(), vmax=50*chi2.min(),s=0.1,**kwargs)
            ax[i].set_xlim(0, np.max(xf))
            
            txt=pnames[i]
            if punits[i].to_string()!="":
                txt+=" ("+punits[i].to_string()+")"
            
            
            ax[i].set_ylabel(txt)
            ax[i].yaxis.set_label_coords(-0.1, 0.5)

        fig.colorbar(scale, ax=ax.ravel().tolist(),label="$\\chi^2_r$ ")

        ax[-1].set_xlabel("step number");

        if savefig!=None:
            plt.savefig(savefig)
    
    