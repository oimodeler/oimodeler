# -*- coding: utf-8 -*-
"""model fitting"""
from multiprocessing import Pool

import corner
import emcee
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import minimize

from .oimParam import oimParam
from .oimSimulator import oimSimulator


class oimFitter(object):
    params = {}

    def __init__(self, *args, **kwargs):
        nargs = len(args)
        if nargs == 2:
            self.simulator = oimSimulator(args[0], args[1])
        elif nargs == 1:
            self.simulator = args[0]
        else:
            raise TypeError("Wrong number of arguments")


        try:
            self.dataTypes = kwargs.pop("dataTypes")
        except:
            self.dataTypes = None

        self.data = self.simulator.data
        self.model = self.simulator.model
        self.pool = None
        self.isPrepared = False

        self._eval(**kwargs)

    def _eval(self, **kwargs):

        for key in self.params.keys():
            if key in kwargs.keys():
                self.params[key].value = kwargs.pop(key)
        return kwargs

    def prepare(self, **kwargs):
        self.freeParams = self.model.getFreeParameters()
        self.nfree = len(self.freeParams)

        self.limits = {}
        for key in self.freeParams:
            self.limits[key] = (self.freeParams[key].min,
                                self.freeParams[key].max)

        kwargs = self._eval(**kwargs)

        self.simulator.prepareData()
        kwargs = self._prepare(**kwargs)
        self.isPrepared = True
        return kwargs

    def run(self, **kwargs):
        if self.isPrepared == False:
            raise TypeError("Fitter not initialized")
        self._run(**kwargs)
        return kwargs

    def getResults(self,**kwargs):
        return 0

    def _prepare(self, **kwargs):
        return kwargs

    def _run(self, **kwargs):
        return kwargs


class oimFitterEmcee(oimFitter):
    def __init__(self, *args, **kwargs):

        self.params["nwalkers"] = oimParam(name="nwalkers", value=16, mini=1,
                                           description="Number of walkers")

        super().__init__(*args, **kwargs)

    def _prepare(self, **kwargs):

        if not ('init' in kwargs):
            init = 'random'
        else:
            init = kwargs.pop('init')
        if init == "random":
            self.initialParams = self._initRandom()
        elif init == "fixed":
            self.initialParams = self._initFixed()
        elif init == "gaussian":
            self.initialParams = self._initGaussian()

        if "moves" in kwargs:
            moves = kwargs.pop['moves']
        else:
            moves = [(emcee.moves.DEMove(), 0.8),
                     (emcee.moves.DESnookerMove(), 0.2)]

        if not ('samplerFile' in kwargs):
            samplerFile = None
        else :
            samplerFile = kwargs.pop('samplerFile')

        if  samplerFile == None:
            self.sampler = emcee.EnsembleSampler(self.params["nwalkers"].value,
                                 self.nfree, self._logProbability, moves=moves,
                                      **kwargs)
        else :
            backend = emcee.backends.HDFBackend(samplerFile)
            self.sampler = emcee.EnsembleSampler(self.params["nwalkers"].value,
                     self.nfree, self._logProbability, moves=moves,
                            backend=backend, **kwargs)
        return kwargs

    def _initGaussian(self):
        nw = self.params['nwalkers'].value
        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = \
                np.random.normal(parami.value,parami.error,
                                 self.params['nwalkers'].value)
        return initialParams

    def _initRandom(self):
        nw = self.params['nwalkers'].value
        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.random.random(
                self.params['nwalkers'].value)*(parami.max-parami.min)+parami.min

        return initialParams

    def _initFixed(self):
        nw = self.params['nwalkers'].value

        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.ones(nw)*parami.value

        return initialParams

    def _run(self, **kwargs):
        self.sampler.run_mcmc(self.initialParams, **kwargs)
        self.getResults()

        return kwargs

    def _logProbability(self, theta):
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]

        for i, key in enumerate(self.freeParams):
            val = theta[i]
            low, up = self.limits[key]
            if not low < val < up:
                return -np.inf

        self.simulator.compute(computeChi2=True,dataTypes = self.dataTypes)
        return -0.5 * self.simulator.chi2r

    def getResults(self, mode='best', discard=0, chi2limfact=20, **kwargs):
        chi2 = -2*self.sampler.get_log_prob(discard=discard, flat=True)
        chain = self.sampler.get_chain(discard=discard, flat=True)

        idx = np.where(chi2 <= chi2limfact*chi2.min())[0]
        chain2 = chain[idx, :]

        if mode == 'best':
            idx = np.argmin(chi2)
            res = chain[idx]
        elif mode == 'mean':
            res = np.mean(chain2, axis=0)
        elif mode == 'median':
            res = np.median(chain2, axis=0)
        else:
            raise NameError("'mode' should be either 'best', 'mean' or 'median'")

        nparam = chain.shape[1]
        err_m = np.ndarray([nparam])
        err_p = np.ndarray([nparam])
        for iparam in range(nparam):
            err_m[iparam] = np.abs(np.quantile(
                chain2[:, iparam], 0.16)-res[iparam])
            err_p[iparam] = np.abs(np.quantile(
                chain2[:, iparam], 0.84)-res[iparam])

        err = 0.5*(err_m+err_p)

        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = res[iparam]
            parami.error = err[iparam]

        self.simulator.compute(computeChi2=True, computeSimulatedData=True
                               ,dataTypes=self.dataTypes)

        return res, err, err_m, err_p

    def cornerPlot(self, discard=0, chi2limfact=20, savefig=None, **kwargs):
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]


        kwargs0 = dict(quantiles=[ 0.16, 0.5, 0.84], show_titles=True, bins=50,
            smooth=2, smooth1d=2, fontsize=8, title_kwargs={'fontsize': 8},
            use_math_text=True)

        kwargs = {**kwargs0,**kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += " ("+uniti.to_string()+")"
            labels.append(txt)

        c = self.sampler.get_chain(discard=discard, flat=True)
        chi2 = -2*self.sampler.get_log_prob(discard=discard, flat=True)

        idx = np.where(chi2 < chi2limfact*chi2.min())[0]
        c2 = c[idx, :]

        fig = corner.corner(c2, labels=labels,**kwargs )

        if savefig != None:
            plt.savefig(savefig)

        return fig, fig.axes

    def walkersPlot(self, savefig=None, chi2limfact=20, labelsize=10,
                       ncolors=128, **kwargs):
        fig, ax = plt.subplots(self.nfree, figsize=(10, 7), sharex=True)
        if self.nfree == 1:
            ax = np.array([ax])

        samples = self.sampler.get_chain()
        chi2 = -2*self.sampler.get_log_prob()
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]
        x = np.arange(chi2.shape[0])

        samples = self.sampler.get_chain()

        xm = np.outer(x, np.ones(chi2.shape[1]))
        chi2f = chi2.flatten()
        xf = xm.flatten()
        idx = np.argsort(-1*chi2f)

        chi2f = chi2f[idx]
        xf = xf[idx]
        samples = samples.reshape(
            [samples.shape[0]*samples.shape[1], samples.shape[2]])[idx, :]

        chi2min=chi2f.min()
        chi2max=chi2limfact*chi2min
        chi2bins=np.linspace(chi2max,chi2min,ncolors)
        if 'cmap' in kwargs:
            cmap=cm.get_cmap(kwargs.pop('cmap'), ncolors)
        else:
            cmap=cm.get_cmap(mpl.rcParams['image.cmap'],ncolors)

        for i in range(self.nfree):
            for icol in range(ncolors):
                if icol!=0:
                    idxi=np.where((chi2f>chi2bins[icol]) & 
                                  (chi2f<=chi2bins[icol-1]))[0]
                else:
                    idxi=np.where(chi2f>chi2bins[icol])[0]

                ax[i].scatter(xf[idxi], samples[idxi, i],
                              color=cmap(ncolors-icol-1),
                              marker=".", s=0.1, **kwargs)

            ax[i].set_xlim(0, np.max(xf))

            txt = pnames[i]

            unit_txt = ""
            if punits[i].to_string() != "":
                unit_txt += " ("+punits[i].to_string()+")"

            ax[i].set_ylabel(unit_txt)

            c = (np.max(samples[:, i])+np.min(samples[:, i]))/2
            ax[i].text(0.02*np.max(xf), c, txt, va="center", ha="left", fontsize=labelsize,
                       backgroundcolor=(1, 1, 1, 0.5))
            ax[i].yaxis.set_label_coords(-0.1, 0.5)

        norm = mpl.colors.Normalize(vmin=chi2min, vmax=chi2max)
        sm =cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm,ax=ax.ravel().tolist(), label="$\\chi^2_r$ ")
        #fig.colorbar(scale, ax=ax.ravel().tolist(), label="$\\chi^2_r$ ")

        ax[-1].set_xlabel("step number")

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax


    def printResults(self,format=".5f",**kwargs):
    
        res=self.getResults(**kwargs)
        chi2r=self.simulator.chi2r
        pm=u'\xb1'
        for iparam,parami in enumerate(self.freeParams):
            print(f"{parami} = {res[0][iparam]:{format}} {pm} {res[3][iparam]:{format}} {self.freeParams[parami].unit}")
        print(f"chi2r = {chi2r:{format}}")
        
        

class oimFitterMinimize(oimFitter):
    def __init__(self, *args, **kwargs):

        self.params["method"] = oimParam(name="method", value=None, mini=1,
                                           description="Number of walkers")
        super().__init__(*args, **kwargs)
    
    def _prepare(self, **kwargs):
        self.initialParams = []
        if not ('initialParams' in kwargs):
            for iparam, parami in enumerate(self.freeParams.values()):
                self.initialParams.append(parami.value) 
        else:
            self.initialParams = kwargs['initialParams']
        return kwargs    
        
    def _getChi2r(self, theta):
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]
        self.simulator.compute(computeChi2=True,dataTypes = self.dataTypes)
        return self.simulator.chi2r
            
    def _run(self, **kwargs):
        self.res=minimize(self._getChi2r,self.initialParams)
        self.getResults()

        return kwargs
    
    def getResults(self,**kwargs):
        values=self.res.x
        errors=np.diag(self.res.hess_inv)**0.5
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = values[iparam]
            parami.error = errors[iparam]

        self.simulator.compute(computeChi2=True, computeSimulatedData=True
                               ,dataTypes=self.dataTypes)
        
        return values, errors
    
    def printResults(self,format=".5f",**kwargs):
    
        res=self.getResults(**kwargs)
        chi2r=self.simulator.chi2r
        pm=u'\xb1'
        nfree = len(self.freeParams)
        print(f"Model : {self.model.shortname}")
        print(f"{nfree} free parameters fitted:")
        for iparam,parami in enumerate(self.freeParams):
            print(f"{parami} = {res[0][iparam]:{format}} {pm} {res[1][iparam]:{format}} {self.freeParams[parami].unit}")
        print(f"chi2r = {chi2r:{format}}")
        
        
        
