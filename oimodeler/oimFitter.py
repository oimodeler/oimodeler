# -*- coding: utf-8 -*-
"""model fitting"""
# from multiprocessing import Pool

import astropy.units as unit
import corner
import emcee
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from dynesty import DynamicNestedSampler, NestedSampler
from dynesty import plotting as dyplot
from matplotlib import cm
from scipy.optimize import least_squares, minimize
from tqdm import tqdm

from .oimParam import oimParam
from .oimSimulator import oimSimulator


class oimFitter:
    params = {}
    description = "Abstract class for model-fitting"

    def __init__(self, *args, **kwargs):

        self.cprior = kwargs.get("cprior", None)
        nargs = len(args)
        if nargs == 2:
            self.simulator = oimSimulator(args[0], args[1], cprior=self.cprior)
        elif nargs == 1:
            self.simulator = args[0]
        else:
            raise TypeError("Wrong number of arguments")

        self.dataTypes = kwargs.pop("dataTypes", None)
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
            self.limits[key] = (
                self.freeParams[key].min,
                self.freeParams[key].max,
            )

        kwargs = self._eval(**kwargs)

        self.simulator.prepareData()
        kwargs = self._prepare(**kwargs)
        self.isPrepared = True
        return kwargs

    def run(self, **kwargs):
        if not self.isPrepared:
            raise TypeError("Fitter not initialized")
        self._run(**kwargs)
        return kwargs

    def getResults(self, **kwargs):
        return 0

    def printResults(self, format=".5f", **kwargs):
        res = self.getResults(**kwargs)
        chi2r = self.simulator.chi2r
        pm = "\xb1"
        for iparam, parami in enumerate(self.freeParams):
            print(
                f"{parami} = {res[0][iparam]:{format}} {pm} {res[1][iparam]:{format}} {self.freeParams[parami].unit}"
            )
        print(f"chi2r = {chi2r:{format}}")

    def _prepare(self, **kwargs):
        return kwargs

    def _run(self, **kwargs):
        return kwargs


class oimFitterEmcee(oimFitter):
    description = "MCMC sampler based on the emcee python module"

    def __init__(self, *args, **kwargs):
        self.params["nwalkers"] = oimParam(
            name="nwalkers", value=16, mini=1, description="Number of walkers"
        )

        super().__init__(*args, **kwargs)

    def _prepare(self, **kwargs):
        init = kwargs.pop("init", "random")
        if init == "random":
            self.initialParams = self._initRandom()
        elif init == "fixed":
            self.initialParams = self._initFixed()
        elif init == "gaussian":
            self.initialParams = self._initGaussian()

        moves = kwargs.pop(
            "moves",
            [
                (emcee.moves.DEMove(), 0.8),
                (emcee.moves.DESnookerMove(), 0.2),
            ],
        )

        samplerFile = kwargs.pop("samplerFile", None)
        if samplerFile is None:
            self.sampler = emcee.EnsembleSampler(
                self.params["nwalkers"].value,
                self.nfree,
                self._logProbability,
                moves=moves,
                **kwargs,
            )
        else:
            backend = emcee.backends.HDFBackend(samplerFile)
            self.sampler = emcee.EnsembleSampler(
                self.params["nwalkers"].value,
                self.nfree,
                self._logProbability,
                moves=moves,
                backend=backend,
                **kwargs,
            )
        return kwargs

    def _initGaussian(self):
        nw = self.params["nwalkers"].value
        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.random.normal(
                parami.value, parami.error, self.params["nwalkers"].value
            )
        return initialParams

    def _initRandom(self):
        nw = self.params["nwalkers"].value
        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = (
                np.random.random(self.params["nwalkers"].value)
                * (parami.max - parami.min)
                + parami.min
            )

        return initialParams

    def _initFixed(self):
        nw = self.params["nwalkers"].value
        initialParams = np.ndarray([nw, self.nfree])

        for iparam, parami in enumerate(self.freeParams.values()):
            initialParams[:, iparam] = np.ones(nw) * parami.value

        return initialParams

    def _run(self, **kwargs):
        if kwargs.get("reset", False):
            state = self.initialParams
        else:
            if self.sampler.iteration == 0:
                state = self.initialParams
            else:
                state = None
        self.sampler.run_mcmc(state, **kwargs)
        self.getResults()

        return kwargs

    # TODO: Maybe make it possible for end-user to input their own
    # parametrisation
    def _logProbability(self, theta):
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]

        for i, key in enumerate(self.freeParams):
            val = theta[i]
            lower, upper = self.limits[key]
            if not lower < val < upper:
                return -np.inf

        self.simulator.compute(
            computeChi2=True, dataTypes=self.dataTypes, cprior=self.cprior
        )
        return -0.5 * self.simulator.chi2

    def getResults(self, mode="best", discard=0, chi2limfact=20, **kwargs):
        chi2 = -2 * self.sampler.get_log_prob(discard=discard, flat=True)
        chain = self.sampler.get_chain(discard=discard, flat=True)

        idx = np.where(chi2 <= chi2limfact * chi2.min())[0]
        chain2 = chain[idx, :]

        if mode == "best":
            idx = np.argmin(chi2)
            res = chain[idx]
        elif mode == "mean":
            res = np.mean(chain2, axis=0)
        elif mode == "median":
            res = np.median(chain2, axis=0)
        else:
            raise NameError(
                "'mode' should be either 'best', 'mean' or 'median'"
            )

        nparam = chain.shape[1]
        err_m = np.ndarray([nparam])
        err_p = np.ndarray([nparam])
        for iparam in range(nparam):
            err_m[iparam] = np.abs(
                np.quantile(chain2[:, iparam], 0.16) - res[iparam]
            )
            err_p[iparam] = np.abs(
                np.quantile(chain2[:, iparam], 0.84) - res[iparam]
            )

        err = 0.5 * (err_m + err_p)
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = res[iparam]
            parami.error = err[iparam]

        self.simulator.compute(
            computeChi2=True,
            computeSimulatedData=True,
            dataTypes=self.dataTypes,
            cprior=self.cprior,
        )

        return res, err, err_m, err_p

    # TODO: Change the chi2limfact implementation so it doesn't break if
    # over-fitting takes place
    def cornerPlot(self, discard=0, chi2limfact=20, savefig=None, **kwargs):
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            bins=50,
            smooth=2,
            smooth1d=2,
            fontsize=8,
            title_kwargs={"fontsize": 8},
            use_math_text=True,
        )
        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += f" ({uniti.to_string()})"
            labels.append(txt)

        c = self.sampler.get_chain(discard=discard, flat=True)
        chi2 = -2 * self.sampler.get_log_prob(discard=discard, flat=True)

        idx = np.where(chi2 < chi2limfact * chi2.min())[0]
        c2 = c[idx, :]
        if c2.size == 0:
            raise ValueError(
                "The emcee chain does not contain enough valid samples for the corner plot. "
                "Potentially the `chi2limfact` parameter is too stringent or the min and max "
                "values of your priors are switched."
            )

        fig = corner.corner(c2, labels=labels, **kwargs)
        if savefig is not None:
            plt.savefig(savefig)

        return fig, fig.axes

    def walkersPlot(
        self, savefig=None, chi2limfact=20, labelsize=10, ncolors=128, **kwargs
    ):
        fig, ax = plt.subplots(self.nfree, figsize=(10, 7), sharex=True)
        if self.nfree == 1:
            ax = np.array([ax])

        samples = self.sampler.get_chain()
        chi2 = -2 * self.sampler.get_log_prob()
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]
        x = np.arange(chi2.shape[0])

        samples = self.sampler.get_chain()

        xm = np.outer(x, np.ones(chi2.shape[1]))
        chi2f = chi2.flatten() / self.simulator.nelChi2
        xf = xm.flatten()
        idx = np.argsort(-1 * chi2f)

        chi2f = chi2f[idx]
        xf = xf[idx]
        samples = samples.reshape(
            [samples.shape[0] * samples.shape[1], samples.shape[2]]
        )[idx, :]

        chi2min = chi2f.min()
        chi2max = chi2limfact * chi2min
        chi2bins = np.linspace(chi2max, chi2min, ncolors)
        if "cmap" in kwargs:
            cmap = cm.get_cmap(kwargs.pop("cmap"), ncolors)
        else:
            cmap = cm.get_cmap(mpl.rcParams["image.cmap"], ncolors)

        for i in range(self.nfree):
            for icol in range(ncolors):
                if icol != 0:
                    idxi = np.where(
                        (chi2f > chi2bins[icol])
                        & (chi2f <= chi2bins[icol - 1])
                    )[0]

                else:
                    idxi = np.where(chi2f > chi2bins[icol])[0]

                ax[i].scatter(
                    xf[idxi],
                    samples[idxi, i],
                    color=cmap(ncolors - icol - 1),
                    marker=".",
                    s=0.1,
                    **kwargs,
                )

            ax[i].set_xlim(0, np.max(xf))

            txt, unit_txt = pnames[i], ""
            if punits[i].to_string() != "":
                unit_txt += f" ({punits[i].to_string()})"

            ax[i].set_ylabel(unit_txt)

            c = (np.max(samples[:, i]) + np.min(samples[:, i])) / 2
            ax[i].text(
                0.02 * np.max(xf),
                c,
                txt,
                va="center",
                ha="left",
                fontsize=labelsize,
                backgroundcolor=(1, 1, 1, 0.5),
            )
            ax[i].yaxis.set_label_coords(-0.1, 0.5)

        norm = mpl.colors.Normalize(vmin=chi2min, vmax=chi2max)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax.ravel().tolist(), label=r"$\chi^2_r$ ")
        # fig.colorbar(scale, ax=ax.ravel().tolist(), label="$\\chi^2_r$ ")

        ax[-1].set_xlabel("step number")

        if savefig is not None:
            plt.savefig(savefig)

        return fig, ax


class oimFitterDynesty(oimFitter):
    """A multinested fitter that has generally a better coverage of the
    global (than MCMC) parameter space."""

    description = "a dynamic nested sampler based on the dynesty python module"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        samplers = {"dynamic": DynamicNestedSampler, "static": NestedSampler}
        self.method = kwargs.pop("method", "dynamic")
        self.sampler = samplers[self.method]

    def _prepare(self, **kwargs):
        """Prepares the dynesty fitter."""
        samplerFile = kwargs.pop("samplerFile", None)
        sampler_kwargs = {
            "sample": kwargs.pop("sample", "rwalk"),
            "bound": kwargs.pop("bound", "multi"),
            "periodic": kwargs.pop("periodic", None),
            "reflective": kwargs.pop("reflective", None),
        }

        if self.method != "dynamic":
            sampler_kwargs["nlive"] = kwargs.pop("nlive", 1000)

        # TODO: Implement the loading of the sampler
        if samplerFile is None:
            self.sampler = self.sampler(
                self._logProbability,
                self._ptform,
                self.nfree,
                update_interval=self.nfree,
                **sampler_kwargs,
                **kwargs,
            )
        else:
            ...

        return kwargs

    def _run(self, **kwargs):
        print_progress = kwargs.pop("progress", False)
        if self.method == "dynamic":
            run_kwargs = {
                "nlive_batch": kwargs.pop("nlive_batch", 1000),
                "maxbatch": kwargs.pop("maxbatch", 100),
                "dlogz_init": kwargs.pop("dlogz_init", 0.01),
                "nlive_init": kwargs.pop("nlive_init", 1000),
            }
        else:
            run_kwargs = {"dlogz": kwargs.pop("dlogz", 0.01)}

        # TODO: Implement checkpoint file here
        self.sampler.run_nested(
            print_progress=print_progress, **run_kwargs, **kwargs
        )
        self.getResults()
        return kwargs

    # TODO: Maybe make it possible for end-user to input their own
    # parametrisation
    def _ptform(self, uniform_samples: np.ndarray) -> np.ndarray:
        """The transformation for uniform sampled values to the
        uniform parameter space."""
        priors = np.array(
            [(param.min, param.max) for param in self.freeParams.values()]
        )
        return priors[:, 0] + (priors[:, 1] - priors[:, 0]) * uniform_samples

    def _logProbability(self, theta: np.ndarray) -> float:
        """The log probability."""
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]

        self.simulator.compute(
            computeChi2=True, dataTypes=self.dataTypes, cprior=self.cprior
        )
        return -0.5 * self.simulator.chi2r

    def getResults(self, mode="median", **kwargs):
        if mode == "median":
            samples = self.sampler.results.samples
            quantiles = np.percentile(samples, [10, 50, 84], axis=0)
            res = quantiles[1]
            err_m, err_p = np.diff(quantiles, axis=0)
        else:
            raise NameError(
                "'mode' should be either 'best', 'mean' or 'median'"
            )

        err = 0.5 * (err_m + err_p)

        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = res[iparam]
            parami.error = err[iparam]

        self.simulator.compute(
            computeChi2=True,
            computeSimulatedData=True,
            dataTypes=self.dataTypes,
            cprior=self.cprior,
        )

        return res, err, err_m, err_p

    def cornerPlot(self, savefig=None, **kwargs):
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_quantiles=[0.16, 0.5, 0.84],
            use_math_text=True,
            max_n_ticks=3,
        )

        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += " (" + uniti.to_string() + ")"
            labels.append(txt)

        results = self.sampler.results
        fig, _ = dyplot.cornerplot(
            results,
            color="blue",
            truths=np.zeros(len(pnames)),
            labels=labels,
            truth_color="black",
            **kwargs,
        )

        if savefig is not None:
            plt.savefig(savefig)

        return fig, fig.axes

    def walkersPlot(self, savefig=None, **kwargs):
        pnames = list(self.freeParams.keys())
        punits = [p.unit for p in list(self.freeParams.values())]

        kwargs0 = dict(
            quantiles=[0.16, 0.5, 0.84], show_titles=True, use_math_text=True
        )
        kwargs = {**kwargs0, **kwargs}

        labels = []
        for namei, uniti in zip(pnames, punits):
            txt = namei
            if uniti.to_string() != "":
                txt += f" ({uniti.to_string()})"
            labels.append(txt)

        results = self.sampler.results
        fig, ax = dyplot.traceplot(
            results,
            labels=labels,
            truths=np.zeros(len(labels)),
            truth_color="black",
            connect=True,
            trace_cmap="viridis",
            connect_highlight=range(5),
            **kwargs,
        )
        if savefig is not None:
            plt.savefig(savefig)

        return fig, ax


class oimFitterMinimize(oimFitter):
    description = (
        r"a simple :math:`\chi^2` minimizer using the numpy Minimize function"
    )

    def __init__(self, *args, **kwargs):

        self.params["method"] = oimParam(
            name="method",
            value="trf",
            mini=1,
            description="minimization method",
        )
        super().__init__(*args, **kwargs)

    def _prepare(self, **kwargs):
        self.initialParams = []
        if "initialParams" not in kwargs:
            for _, parami in enumerate(self.freeParams.values()):
                self.initialParams.append(parami.value)
        else:
            self.initialParams = kwargs["initialParams"]
        return kwargs

    def _getChi2r(self, theta):
        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = theta[iparam]
            if theta[iparam] < parami.min or theta[iparam] > parami.max:
                return np.inf
        self.simulator.compute(
            computeChi2=True, dataTypes=self.dataTypes, cprior=self.cprior
        )

        return self.simulator.chi2

    def _run(self, **kwargs):

        self.res = minimize(self._getChi2r, self.initialParams, **kwargs)
        # self.res = least_squares(self._getChi2r, self.initialParams,method=
        #                         self.params["method"].value)
        self.getResults()
        return kwargs

    def getResults(self, **kwargs):
        values = self.res.x

        try:
            errors = np.diag(self.res.hess_inv) ** 0.5
        except:
            cov = np.linalg.inv(self.res.jac.T.dot(self.res.jac))
            errors = np.sqrt(np.diagonal(cov))

        for iparam, parami in enumerate(self.freeParams.values()):
            parami.value = values[iparam]
            parami.error = errors[iparam]

        self.simulator.compute(
            computeChi2=True,
            computeSimulatedData=True,
            dataTypes=self.dataTypes,
            cprior=self.cprior,
        )

        return values, errors

    def printResults(self, format=".5f", **kwargs):
        res = self.getResults(**kwargs)
        chi2r = self.simulator.chi2r
        pm = "\xb1"
        for iparam, parami in enumerate(self.freeParams):
            print(
                f"{parami} = {res[0][iparam]:{format}} {pm} {res[1][iparam]:{format}} {self.freeParams[parami].unit}"
            )
        print(f"chi2r = {chi2r:{format}}")


class oimFitterRegularGrid(oimFitter):
    description = r"regular grid with :math:`\chi^2` explorer"

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def _prepare(self, **kwargs):

        if "params" in kwargs:
            self.gridParams = kwargs["params"]
        else:
            self.gridParams = list(self.model.getFreeParameters().values())

        if "max" in kwargs:
            self.gridMax = kwargs["max"]
        else:
            self.gridMax = []
            for parami in self.gridParams:
                self.gridMax.append(parami.max)
        if "min" in kwargs:
            self.gridMin = kwargs["min"]
        else:
            self.gridMin = []
            for parami in self.gridParams:
                self.gridMin.append(parami.min)
        if "steps" in kwargs:
            self.gridSteps = kwargs["steps"]
        else:
            raise TypeError("steps argument needed in oimFitterGrid")

        self.grid = []
        self.gridSize = []
        if (
            len(self.gridParams) == len(self.gridMax)
            and len(self.gridParams) == len(self.gridMin)
            and len(self.gridParams) == len(self.gridSteps)
        ):
            for iparam in range(len(self.gridMax)):
                gi = np.arange(
                    self.gridMin[iparam],
                    self.gridMax[iparam] + self.gridSteps[iparam],
                    self.gridSteps[iparam],
                )

                self.gridSize.append(gi.size)
                self.grid.append(gi)
        else:

            raise TypeError(
                "Steps, min, max and gridParams arrays should have the same sizes"
            )

        return kwargs

    def _run(self, **kwargs):
        self.chi2rMap = np.zeros(self.gridSize)

        n = self.chi2rMap.size

        for i in tqdm(range(n)):
            igrid = np.unravel_index(i, self.gridSize)
            for iparam in range(len(igrid)):
                self.gridParams[iparam].value = self.grid[iparam][
                    igrid[iparam]
                ]
            self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
            self.chi2rMap[igrid] = self.simulator.chi2r
        return kwargs

    def getResults(self, **kwargs):

        idx_best = np.unravel_index(np.argmin(self.chi2rMap), self.gridSize)

        best = []
        for iparam in range(len(idx_best)):
            self.gridParams[iparam].value = self.grid[iparam][idx_best[iparam]]
            best.append(self.gridParams[iparam].value)

        self.simulator.compute(
            computeChi2=True,
            computeSimulatedData=True,
            dataTypes=self.dataTypes,
            cprior=self.cprior,
        )

        return best

    def printResults(self, format=".5f", **kwargs):
        res = self.getResults(**kwargs)
        chi2r = self.simulator.chi2r
        for iparam, parami in enumerate(self.gridParams):
            print(
                f"{parami.name} = {res[iparam]:{format}} {self.gridParams[iparam].unit}"
            )
        print(f"chi2r = {chi2r:{format}}")

    def plotMap(
        self,
        params=None,
        fixedValues="best",
        plotContour=False,
        plotMinLines=False,
        plotMin=True,
        minLines_kwargs={},
        contour_kwargs={},
        clabel_kwargs={},
        min_kwargs={},
        **kwargs,
    ):

        min_kwargs0 = dict(marker="o", color="r")
        min_kwargs = min_kwargs0 | min_kwargs

        minLines_kwargs0 = dict(ls="--", color="r")
        minLines_kwargs = minLines_kwargs0 | minLines_kwargs

        contour_kwargs0 = dict(colors="r", levels=[2])
        contour_kwargs = contour_kwargs0 | contour_kwargs

        clabel_kwargs0 = dict(inline=True, fmt="%.1f", fontsize=10)
        clabel_kwargs = clabel_kwargs0 | clabel_kwargs

        if not ("origin" in kwargs):
            kwargs["origin"] = "lower"
        if not ("aspect" in kwargs):
            kwargs["aspect"] = "auto"

        ndims = len(self.gridSize)

        min_idx = np.argmin(self.chi2rMap)
        chi2rmin = np.min(self.chi2rMap)
        chi2rmax = np.max(self.chi2rMap)
        if ndims == 1:
            fig, ax = plt.subplots()
            xmin = self.grid[0][min_idx]

            ax.plot(self.grid[0], self.chi2rMap, color="r")
            if plotMinLines:
                label = (
                    rf"min $\chi^2_r$ = {chi2rmin:.1f}"
                    f" at {self.gridParams[0].name}={xmin}"
                    f" {self.gridParams[0].unit.to_string(format='latex')} "
                )
                ax.plot(
                    [xmin, xmin],
                    [chi2rmin, chi2rmax],
                    label=label,
                    **min_kwargs,
                )
                ax.legend()
            txt = self.gridParams[0].name
            if self.gridParams[0].unit != unit.one:
                txt += f"({self.gridParams[0].unit.to_string(format='latex')})"
            ax.set_xlabel(txt)
            ax.set_ylabel("$\\chi^2_r$")

        else:
            if len(self.gridSize) != 2 and params == None:
                raise TypeError(
                    "The 2 parameters to plot should be specified "
                    "when the grid dimension is higher than 2"
                )
            elif ndims == 2 and params == None:

                im = self.chi2rMap

                im_params = [self.gridParams[0], self.gridParams[1]]
                dim1 = 0
                dim2 = 1

            elif ndims > 2 and params != None:
                dim1 = np.where(np.array(self.gridParams) == params[0])[0][0]
                dim2 = np.where(np.array(self.gridParams) == params[1])[0][0]

                mins = np.unravel_index(
                    np.argmin(self.chi2rMap), self.gridSize
                )

                # That is very ugly but I don't know how to programatically select axes
                txt = "self.chi2rMap["
                for idim in range(ndims):
                    if idim == dim1:
                        txt += ":,"
                    elif idim == dim2:
                        txt += ":,"
                    else:
                        txt += str(mins[idim])
                        txt += ","
                txt = txt[:-1]
                txt += "]"
                im = eval(txt)

                if dim1 > dim2:
                    im = np.transpose(im)

                im_params = params

            min_idx = np.unravel_index(min_idx, self.chi2rMap.shape)

            fig, ax = plt.subplots()

            contour_kwargs["levels"] = (
                np.array(contour_kwargs["levels"]) * chi2rmin
            )

            sm = ax.imshow(
                np.transpose(im),
                extent=[
                    self.grid[dim1][0],
                    self.grid[dim1][-1],
                    self.grid[dim2][0],
                    self.grid[dim2][-1],
                ],
                **kwargs,
            )
            if plotContour:
                cs = ax.contour(
                    self.grid[dim1],
                    self.grid[dim2],
                    np.transpose(im),
                    **contour_kwargs,
                )
                ax.clabel(cs, cs.levels, **clabel_kwargs)

            xmin = self.grid[dim1][min_idx[0]]
            ymin = self.grid[dim2][min_idx[1]]

            if plotMinLines:
                ax.plot(
                    [self.grid[dim1][0], self.grid[dim1][-1]],
                    [ymin, ymin],
                    **minLines_kwargs,
                )
                ax.plot(
                    [xmin, xmin],
                    [self.grid[dim2][0], self.grid[dim2][-1]],
                    **minLines_kwargs,
                )

            if plotMin:

                label = (
                    rf"min $\chi^2_r$ = {chi2rmin:.1f}"
                    f" at {im_params[0].name}={xmin:.1f}"
                    f" {im_params[0].unit.to_string(format='latex')}"
                    f" and {im_params[1].name}={ymin:.1f}"
                    f" {im_params[1].unit.to_string(format='latex')}"
                )
                ax.scatter(xmin, ymin, label=label, **min_kwargs)
                ax.legend()

            xlabel = im_params[0].name
            if im_params[0].unit != unit.one:
                xlabel += " (" + im_params[0].unit.to_string() + ")"

            ylabel = im_params[1].name
            if im_params[1].unit != unit.one:
                ylabel += " (" + im_params[1].unit.to_string() + ")"

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            fig.colorbar(sm, ax=ax, label="$\\chi^2_r$")

        return fig, ax
