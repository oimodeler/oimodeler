# -*- coding: utf-8 -*-
"""model fitting"""
from typing import Dict, Optional, Tuple
from pathlib import Path

import corner
import emcee
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.optimize import minimize
from dynesty import DynamicNestedSampler, NestedSampler
from dynesty import plotting as dyplot

from .oimSimulator import oimSimulator


class oimFitter:
    """A class that fits the model to the data.

    Attributes
    ----------
    simulator : .oimSimulator.oimSimulator
        The simulator object that contains the data and the model.
    dataTypes : list
        The list of data types to be used in the computation of the chi2
    data : .oimData.oimData
        The data object that contains the data.
    model : .oimModel.oimModel
        The model object that will be fitted to the data.
    isPrepared : bool
        A boolean that indicates if the fitter has been prepared.
    params : dict
        The parameters of the fitter.
    freeParams : dict
        A dictionary of the free parameters of the model.
    nfree : int
        The number of free parameters of the model.
    limits : dict
        The limits of the free parameters.
    """
    params = {}

    def __init__(self, *args, **kwargs):
        nargs = len(args)
        if nargs == 2:
            self.simulator = oimSimulator(args[0], args[1])
        elif nargs == 1:
            self.simulator = args[0]
        else:
            raise TypeError("Wrong number of arguments")

        self.dataTypes = kwargs.pop("dataTypes", None)
        self.data = self.simulator.data
        self.model = self.simulator.model

        self.isPrepared = False

        self._eval(**kwargs)

    def _eval(self, **kwargs) -> Dict:
        """Evaluate the parameters of the fitter."""
        for key in self.params.keys():
            if key in kwargs.keys():
                self.params[key].value = kwargs.pop(key)

        return kwargs

    def prepare(self, **kwargs) -> Dict:
        """Prepare the fitter for the model-fitting."""
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

    def run(self, **kwargs) -> Dict:
        """Run the model-fitting."""
        if not self.isPrepared:
            raise TypeError("Fitter not initialized")
        self._run(**kwargs)
        return kwargs

    def getResults(self, **kwargs):
        """Get the results of the model-fitting and closes the pool if it is not None.
        """
        return 0

    def printResults(self, format: Optional[str] = ".5f", **kwargs) -> None:
        """Print the results of the model-fitting.

        Parameters
        ----------
        format : str
            The format of the output.
        """
        res = self.getResults(**kwargs)
        chi2r = self.simulator.chi2r
        pm = u'\xb1'
        for iparam,parami in enumerate(self.freeParams):
            print(f"{parami} = {res[0][iparam]:{format}} {pm} {res[3][iparam]:{format}} {self.freeParams[parami].unit}")
        print(f"chi2r = {chi2r:{format}}")

    def _prepare(self, **kwargs):
        """Prepare the fitter for the model-fitting."""
        return kwargs

    def _run(self, **kwargs):
        """Run the model-fitting."""
        return kwargs


class oimFitterEmcee(oimFitter):
    """The emcee fitter. It uses the emcee package to fit the model to the data.

    Attributes
    ----------
    simulator : .oimSimulator.oimSimulator
        The simulator object that contains the data and the model.
    dataTypes : list
        The list of data types to be used in the computation of the chi2
    data : .oimData.oimData
        The data object that contains the data.
    model : .oimModel.oimModel
        The model object that will be fitted to the data.
    isPrepared : bool
        A boolean that indicates if the fitter has been prepared.
    params : dict
        The parameters of the fitter.
    freeParams : dict
        A dictionary of the free parameters of the model.
    nfree : int
        The number of free parameters of the model.
    nwalkers : int
        The number of walkers.
    limits : dict
        The limits of the free parameters.
    initialParams : np.ndarray
        The initial parameters of the model.
    sampler : emcee.EnsembleSampler
        The emcee sampler.
    """
    def _prepare(self, **kwargs) -> Dict:
        """Prepare the emcee fitter.

        Parameters
        ----------
        nwalkers : int, optional
            The number of walkers. Default is 16.
        pool : multiprocessing.Pool, optional
            The pool of workers. Default is None.
        init : str, optional
            The initialisation method. It can be either 'random', 'fixed' or 'gaussian'.
            Default is 'random'.
        """
        self.nwalkers = kwargs.pop("nwalkers", 16)
        pool = kwargs.pop("pool", None)

        init = kwargs.pop("init", "random")
        if init == "random":
            self.initialParams = self._initRandom()
        elif init == "fixed":
            self.initialParams = self._initFixed()
        elif init == "gaussian":
            self.initialParams = self._initGaussian()

        moves = kwargs.pop("moves",
                           [(emcee.moves.DEMove(), 0.8),
                            (emcee.moves.DESnookerMove(), 0.2)])

        samplerFile = kwargs.pop("samplerFile", None)
        if samplerFile is None:
            self.sampler = emcee.EnsembleSampler(
                    self.nwalkers, self.nfree,
                    self._logProbability,
                    pool=pool, moves=moves, **kwargs)
        else:
            backend = emcee.backends.HDFBackend(samplerFile)
            self.sampler = emcee.EnsembleSampler(
                    self.nwalkers, self.nfree,
                    self._logProbability,
                    pool=pool, moves=moves,
                    backend=backend, **kwargs)
        return kwargs

    def _initGaussian(self) -> np.ndarray:
        """Initialise the parameters with a Gaussian distribution."""
        initialParams = np.ndarray([self.nwalkers, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = \
                np.random.normal(parami.value, parami.error, self.nwalkers)
        return initialParams

    def _initRandom(self) -> np.ndarray:
        """Initialise the parameters with a random distribution."""
        initialParams = np.ndarray([self.nwalkers, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.random.random(
                self.nwalkers)*(parami.max-parami.min)+parami.min

        return initialParams

    def _initFixed(self) -> np.ndarray:
        """Initialise the parameters with a fixed distribution."""
        initialParams = np.ndarray([self.nwalkers, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.ones(self.nwalkers)*parami.value

        return initialParams

    def _run(self, **kwargs) -> Dict:
        """Run the model-fitting."""
        self.sampler.run_mcmc(self.initialParams, **kwargs)
        return kwargs

    # TODO: Maybe make it possible for end-user to input their own
    # parametrisation
    def _logProbability(self, theta: np.ndarray):
        """The log probability."""
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]

        for i, key in enumerate(self.freeParams):
            val = theta[i]
            low, up = self.limits[key]
            if not low < val < up:
                return -np.inf

        self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
        return -0.5 * self.simulator.chi2r

    def getResults(self, mode: Optional[str] = "best",
                   discard: Optional[int] = 0,
                   chi2limfact: Optional[float] = 20.,
                   **kwargs) -> Tuple[np.ndarray]:
        """Get the results of the model-fitting.

        Parameters
        ----------
        mode : str, optional
            The mode of the results. It can be either 'best', 'mean' or 'median'.
            Default is 'best'.
        discard : int, optional
            The number of steps to discard. Default is 0.
        chi2limfact : float, optional
            The factor of the minimum chi2r to consider. Default is 20.

        Returns
        -------
        res : np.ndarray
            The results of the model-fitting.
        err : np.ndarray
            The errors of the results.
        err_m : np.ndarray
            The negative errors of the results.
        err_p : np.ndarray
            The positive errors of the results.
        """
        chi2 = -2*self.sampler.get_log_prob(discard=discard, flat=True)
        chain = self.sampler.get_chain(discard=discard, flat=True)

        idx = np.where(chi2 <= chi2limfact*chi2.min())[0]
        chain2 = chain[idx, :]

        if mode == "best":
            idx = np.argmin(chi2)
            res = chain[idx]
        elif mode == "mean":
            res = np.mean(chain2, axis=0)
        elif mode == "median":
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

        self.simulator.compute(computeChi2=True, computeSimulatedData=True,
                               dataTypes=self.dataTypes)

        return res, err, err_m, err_p

    def cornerPlot(self, discard: Optional[int] = 0,
                   chi2limfact: Optional[float] = 20.,
                   savefig: Optional[Path] = None, **kwargs) -> Tuple[Figure, Axes]:
        """Make a corner plot of the results.

        Parameters
        ----------
        discard : int, optional
            The number of steps to discard. Default is 0.
        chi2limfact : float, optional
            The factor of the minimum chi2r to consider. Default is 20.
        savefig : str, optional
            The name of the file to save the figure. Default is None.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure of the corner plot.
        ax : matplotlib.axes.Axes
            The axes of the corner plot.
        """
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(quantiles=[0.16, 0.5, 0.84], show_titles=True,
                       bins=50, smooth=2, smooth1d=2, fontsize=8,
                       title_kwargs={'fontsize': 8}, use_math_text=True)
        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += f" ({uniti.to_string()})"
            labels.append(txt)

        c = self.sampler.get_chain(discard=discard, flat=True)
        chi2 = -2*self.sampler.get_log_prob(discard=discard, flat=True)

        idx = np.where(chi2 < chi2limfact*chi2.min())[0]
        c2 = c[idx, :]

        fig = corner.corner(c2, labels=labels, **kwargs)

        if savefig != None:
            plt.savefig(savefig)

        return fig, fig.axes

    def walkersPlot(self, savefig: Optional[Path] = None,
                    chi2limfact: Optional[float] = 20.,
                    labelsize: Optional[int] = 10,
                    ncolors: Optional[int] = 128, **kwargs) -> Tuple[Figure, Axes]:
        """Make a walkers plot of the results.

        Parameters
        ----------
        savefig : str, optional
            The name of the file to save the figure. Default is None.
        chi2limfact : float, optional
            The factor of the minimum chi2r to consider. Default is 20.
        labelsize : int, optional
            The size of the labels. Default is 10.
        ncolors : int, optional
            The number of colors. Default is 128.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure of the walkers plot.
        ax : matplotlib.axes.Axes
            The axes of the walkers plot.
        """
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

        chi2min = chi2f.min()
        chi2max = chi2limfact*chi2min
        chi2bins = np.linspace(chi2max, chi2min, ncolors)
        if "cmap" in kwargs:
            cmap = cm.get_cmap(kwargs.pop("cmap"), ncolors)
        else:
            cmap = cm.get_cmap(mpl.rcParams["image.cmap"], ncolors)

        for i in range(self.nfree):
            for icol in range(ncolors):      
                if icol != 0:
                    idxi = np.where((chi2f > chi2bins[icol])
                                    & (chi2f <= chi2bins[icol-1]))[0]

                else:
                    idxi = np.where(chi2f > chi2bins[icol])[0]

                ax[i].scatter(xf[idxi], samples[idxi, i],
                              color=cmap(ncolors-icol-1),
                              marker=".", s=0.1, **kwargs)

            ax[i].set_xlim(0, np.max(xf))

            txt, unit_txt = pnames[i], ""
            if punits[i].to_string() != "":
                unit_txt += f" ({punits[i].to_string()})"

            ax[i].set_ylabel(unit_txt)

            c = (np.max(samples[:, i])+np.min(samples[:, i]))/2
            ax[i].text(0.02*np.max(xf), c, txt, va="center", ha="left",
                       fontsize=labelsize, backgroundcolor=(1, 1, 1, 0.5))
            ax[i].yaxis.set_label_coords(-0.1, 0.5)

        norm = mpl.colors.Normalize(vmin=chi2min, vmax=chi2max)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax.ravel().tolist(), label=r"$\chi^2_r$ ")
        # fig.colorbar(scale, ax=ax.ravel().tolist(), label="$\\chi^2_r$ ")

        ax[-1].set_xlabel("step number")

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax
        

class oimFitterDynesty(oimFitter):
    """The dynesty fitter. It uses the dynesty package to fit the model to the data.

    Attributes
    ----------
    simulator : .oimSimulator.oimSimulator
        The simulator object that contains the data and the model.
    dataTypes : list
        The list of data types to be used in the computation of the chi2
    data : .oimData.oimData
        The data object that contains the data.
    model : .oimModel.oimModel
        The model object that will be fitted to the data.
    isPrepared : bool
        A boolean that indicates if the fitter has been prepared.
    params : dict
        The parameters of the fitter.
    freeParams : dict
        A dictionary of the free parameters of the model.
    nfree : int
        The number of free parameters of the model.
    limits : dict
        The limits of the free parameters.
    method : str
        The method used for sampling, either "dynamic" or "static".
    sampler : dynesty.NestedSampler or dynesty.DynamicNestedSampler
        The dynesty sampler.
    """
    def _prepare(self, **kwargs) -> Dict:
        """Prepares the dynesty fitter.

        Parameters
        ----------
        method : str, optional
            The method used for sampling, either "dynamic" or "static".
            Default is "static".
        ncores : int, optional
            The number of cores. Default is 1.
        pool : multiprocessing.Pool, optional
            The pool of workers. Default is None.
        nlive : int, optional
            The number of live points. Default is 1000.
        sample : str, optional
            The sample method. Default is "rwalk".
        bound : str, optional
            The bound method. Default is "multi".
        samplerFile : str, optional
            The name of the sampler file. Default is None.
        """
        self.method = kwargs.pop("method", "static")
        samplers = {"dynamic": DynamicNestedSampler, "static": NestedSampler}
        self.sampler = samplers[self.method]

        ncores = kwargs.pop("ncores", 1)
        pool = kwargs.pop("pool", None)
        nlive = kwargs.pop("nlive", 1000)
        sample = kwargs.pop("sample", "rwalk")
        bound = kwargs.pop("bound", "multi")

        # TODO: Implement the loading of the sampler
        samplerFile = kwargs.pop("samplerFile", None)
        if samplerFile is None:
            self.sampler = self.sampler(
                    self._logProbability, self._ptform, self.nfree,
                    nlive=nlive, sample=sample, bound=bound,
                    queue_size=ncores, pool=pool,
                    update_interval=self.nfree, **kwargs)
        else:
            ...

        return kwargs

    def _run(self, **kwargs):
        print_progress = kwargs.pop("progress", False)
        dlogz = kwargs.pop("dlogz", 0.010)

        sampler_kwargs = {"dlogz": dlogz, "print_progress": print_progress}
        if self.method == "dynamic":
            del sampler_kwargs["dlogz"]

        # TODO: Implement checkpoint file here
        self.sampler.run_nested(**sampler_kwargs, **kwargs)
        return kwargs

    # TODO: Maybe make it possible for end-user to input their own
    # parametrisation
    def _ptform(self, uniform_samples: np.ndarray) -> np.ndarray:
        """The transformation for uniform sampled values to the
        uniform parameter space."""
        priors = np.array([(param.min, param.max) for param in self.freeParams.values()])
        return priors[:, 0] + (priors[:, 1] - priors[:, 0])*uniform_samples

    def _logProbability(self, theta: np.ndarray) -> float:
        """The log probability."""
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]

        self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
        return -0.5 * self.simulator.chi2r

    def getResults(self, mode: Optional[str] = "median", **kwargs) -> Tuple[np.ndarray]:
        """Get the results of the model-fitting.

        Parameters
        ----------
        mode : str, optional
            The mode of the results. It can be either 'best', 'mean' or 'median'.
            Default is 'median'.

        Returns
        -------
        res : np.ndarray
            The results of the model-fitting.
        err : np.ndarray
            The errors of the results.
        err_m : np.ndarray
            The negative errors of the results.
        err_p : np.ndarray
            The positive errors of the results.
        """
        if mode == "median":
            samples = self.sampler.results.samples
            quantiles = np.percentile(samples, [10, 50, 84], axis=0)
            res = quantiles[1]
            err_m, err_p = np.diff(quantiles, axis=0)
        else:
            raise NameError("'mode' should be either 'best', 'mean' or 'median'")

        err = 0.5*(err_m+err_p)

        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = res[iparam]
            parami.error = err[iparam]

        self.simulator.compute(computeChi2=True, computeSimulatedData=True,
                               dataTypes=self.dataTypes)

        return res, err, err_m, err_p

    def cornerPlot(self, savefig: Optional[Path] = None, **kwargs) -> Tuple[Figure, Axes]:
        """Make a corner plot of the results.

        Parameters
        ----------
        savefig : str, optional
            The name of the file to save the figure. Default is None.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure of the corner plot.
        ax : matplotlib.axes.Axes
            The axes of the corner plot.
        """
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(quantiles=[0.16, 0.5, 0.84], show_titles=True,
                       title_quantiles=[0.16, 0.5, 0.84],
                       use_math_text=True, max_n_ticks=3)

        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += f" ({uniti.to_string()})"
            labels.append(txt)

        results = self.sampler.results
        fig, _ = dyplot.cornerplot(results, color="blue",
                                   truths=np.zeros(len(pnames)),
                                   labels=labels, truth_color="black", **kwargs)

        if savefig is not None:
            plt.savefig(savefig)

        return fig, fig.axes

    def walkersPlot(self, savefig: Optional[Path] = None, **kwargs) -> Tuple[Figure, Axes]:
        """Make a walkers plot of the results.

        Parameters
        ----------
        savefig : str, optional
            The name of the file to save the figure. Default is None.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure of the walkers plot.
        ax : matplotlib.axes.Axes
            The axes of the walkers plot.
        """
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, use_math_text=True)
        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += f" ({uniti.to_string()})"
            labels.append(txt)

        results = self.sampler.results
        fig, ax = dyplot.traceplot(results, labels=labels,
                                   truths=np.zeros(len(labels)),
                                   truth_color="black", connect=True,
                                   trace_cmap="viridis",
                                   connect_highlight=range(5), **kwargs)
        if savefig is not None:
            plt.savefig(savefig)

        return fig, ax
        

class oimFitterMinimize(oimFitter):
    """The minimize fitter. It uses the scipy.optimize.minimize function
    to fit the model to the data.


    Attributes
    ----------
    simulator : .oimSimulator.oimSimulator
        The simulator object that contains the data and the model.
    dataTypes : list
        The list of data types to be used in the computation of the chi2
    data : .oimData.oimData
        The data object that contains the data.
    model : .oimModel.oimModel
        The model object that will be fitted to the data.
    isPrepared : bool
        A boolean that indicates if the fitter has been prepared.
    params : dict
        The parameters of the fitter.
    freeParams : dict
        A dictionary of the free parameters of the model.
    nfree : int
        The number of free parameters of the model.
    limits : dict
        The limits of the free parameters.
    method : str
        The method used for minimization. Default is None.
    initialParams : np.ndarray
        The initial parameters of the model.
    """
    def _prepare(self, **kwargs) -> Dict:
        """Prepare the minimize fitter.

        Parameters
        ----------
        method : str, optional
            The method used for minimization. Default is None.
        pool : multiprocessing.Pool, optional
            The pool of workers. Default is None.
        initialParams : list, optional
            The initial parameters of the model. Default is [].
        """
        self.method = kwargs.pop("method", None)
        pool = kwargs.pop("pool", None)

        self.initialParams = kwargs.pop("initialParams", [])
        if self.initialParams:
            for _, parami in enumerate(self.freeParams.values()):
                self.initialParams.append(parami.value) 
        return kwargs    
        
    def _getChi2r(self, theta: np.ndarray) -> float:
        """The chi2r function."""
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]
        self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
        return self.simulator.chi2r
            
    def _run(self, **kwargs):
        """Run the model-fitting."""
        self.res = minimize(self._getChi2r, self.initialParams)
        self.getResults()
        return kwargs
    
    def getResults(self, **kwargs) -> Tuple[np.ndarray]:
        """Get the results of the model-fitting.

        Returns
        -------
        values : np.ndarray
            The values of the model's parameters.
        errors : np.ndarray
            The errors of the model's parameters.
        """
        values, errors = self.res.x, np.diag(self.res.hess_inv)**0.5
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = values[iparam]
            parami.error = errors[iparam]

        self.simulator.compute(computeChi2=True, computeSimulatedData=True,
                               dataTypes=self.dataTypes)
        
        return values, errors
