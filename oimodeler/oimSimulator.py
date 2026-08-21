# -*- coding: utf-8 -*-
"""Data/model simulation"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .oimData import oimData, oimDataType
from .oimPlots import (
    _errorplot,
    oimPlotParamArr,
    oimPlotParamError,
    oimPlotParamLabel,
    oimPlotParamLabelShort,
    oimPlotParamName,
    oimWlTemplatePlots,
)
from .oimUtils import hdulistDeepCopy

oimDataArrDict = {
    "OI_FLUX": {"data": ["FLUXDATA"], "err": ["FLUXERR"]},
    "OI_T3": {
        "data": ["T3AMP", "T3PHI"],
        "err": ["T3AMPERR", "T3PHIERR"],
    },
    "OI_VIS": {
        "data": ["VISAMP", "VISPHI"],
        "err": ["VISAMPERR", "VISPHIERR"],
    },
    "OI_VIS2": {"data": ["VIS2DATA"], "err": ["VIS2ERR"]},
}


def corrFlux2Vis2(vcompl):
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB - 1), vcompl[0, :])
    return np.abs(vcompl[1:, :] / norm) ** 2


def corrFlux2VisAmpAbs(vcompl: np.ndarray) -> np.ndarray:
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB - 1), vcompl[0, :])
    return np.abs(vcompl[1:, :] / norm)


# FIXME : Not real formula for differential visibilities
def corrFlux2VisAmpDif(vcompl: np.ndarray) -> np.ndarray:
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    return np.abs(vcompl[1:, :] / norm)


def corrFlux2VisAmpCor(vcompl: np.ndarray) -> np.ndarray:
    return np.abs(vcompl[1:, :])


def corrFlux2VisPhiAbs(vcompl: np.ndarray) -> np.ndarray:
    return np.angle(vcompl[1:, :], deg=True)


# FIXME : Not real formula for differential phases
def corrFlux2VisPhiDif(vcompl: np.ndarray) -> np.ndarray:
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    return np.angle(vcompl[1:, :] * np.conjugate(norm), deg=True)


# TODO: Special function doing T3Amp and T3Phi simultaneously
def corrFlux2T3Amp(vcompl: np.ndarray) -> np.ndarray:
    nB = vcompl.shape[0]
    nCP = (nB - 1) // 3
    norm = np.outer(np.ones(nCP), vcompl[0, :])
    BS = (
        vcompl[1 : nCP + 1, :]
        * vcompl[nCP + 1 : 2 * nCP + 1, :]
        * np.conjugate(vcompl[2 * nCP + 1 :, :])
        / norm**3
    )
    return np.abs(BS)


def corrFlux2T3Phi(vcompl: np.ndarray) -> np.ndarray:
    nB = vcompl.shape[0]
    nCP = (nB - 1) // 3
    norm = np.outer(np.ones(nCP), vcompl[0, :])
    BS = (
        vcompl[1 : nCP + 1, :]
        * vcompl[nCP + 1 : 2 * nCP + 1, :]
        * np.conjugate(vcompl[2 * nCP + 1 :, :])
        / norm**3
    )
    return np.angle(BS, deg=True)


def corrFlux2Flux(vcompl: np.ndarray) -> np.ndarray:
    return np.abs(vcompl)


class oimSimulator:
    """Contains"""

    def __init__(
        self,
        data=None,
        model=None,
        fitter=None,
        cprior=None,
        priorWeight=1,
        **kwargs,
    ):

        self.data = oimData()
        self.simulatedData = None
        self.model = None
        self.cprior = cprior
        self.priorWeight = priorWeight

        if data is not None:
            if isinstance(data, oimData):
                self.data = data
            else:
                self.addData(data)

        if model is not None:
            self.setModel(model)

        if model is not None and data is not None:
            self.compute(
                computeChi2=True, computeSimulatedData=True, cprior=self.cprior
            )

    def setModel(self, model):
        self.model = model

    def addData(self, data):
        self.data.addData(data)

    def prepareData(self):
        self.data.prepareData()
        self.simulatedData = oimData()
        for datai in self.data.data:
            self.simulatedData.addData(hdulistDeepCopy(datai))

    def prepareBootstrap(self):
        self.data.prepareData()
        self.bootstrapData = oimData()
        for datai in self.simulatedData.data:
            dataic = hdulistDeepCopy(datai)

            for dataij in dataic:
                if dataij.name in oimDataArrDict:
                    for dataName, errName in zip(
                        oimDataArrDict[dataij.name]["data"],
                        oimDataArrDict[dataij.name]["err"],
                    ):

                        shape = dataij.data[dataName].shape
                        err = dataij.data[errName]

                        dataij.data[dataName] += np.random.randn(*shape) * err

            self.bootstrapData.addData(dataic)

    def compute(
        self,
        computeChi2: bool = False,
        computeSimulatedData: bool = False,
        checkSimulatedData: bool = True,
        dataTypes=None,
        cprior=None,
    ):
        if dataTypes is None:
            dataTypes = [
                "VIS2DATA",
                "VISAMP",
                "VISPHI",
                "T3AMP",
                "T3PHI",
                "FLUXDATA",
            ]

        self.vcompl = self.model.getComplexCoherentFlux(
            self.data.vect_u,
            self.data.vect_v,
            self.data.vect_wl,
            self.data.vect_mjd,
        )

        nelChi2 = 0
        chi2 = 0
        chi2List = []
        residuals = []

        if computeSimulatedData == True and (
            checkSimulatedData == True or self.simulatedData == None
        ):
            self.simulatedData = oimData()
            for datai in self.data.data:
                self.simulatedData.addData(hdulistDeepCopy(datai))

        data = self.data

        if (computeChi2 == True) | (computeSimulatedData == True):
            idx = 0
            nfiles = len(data.struct_u)
            for ifile in range(nfiles):
                narr = len(data.struct_arrType[ifile])
                for iarr in range(narr):
                    arrNum = data.struct_arrNum[ifile][iarr]
                    arrType = data.struct_arrType[ifile][iarr]
                    dataType = data.struct_dataType[ifile][iarr]
                    nB = data.struct_nB[ifile][iarr]
                    nwl = data.struct_nwl[ifile][iarr]
                    vcompli = self.vcompl[idx : idx + nB * nwl]
                    vcompli = np.reshape(vcompli, [nB, nwl])

                    dataVal = data.struct_val[ifile][iarr]
                    dataErr = data.struct_err[ifile][iarr]
                    flag = data.struct_flag[ifile][iarr]

                    idx += nB * nwl
                    quantities = []
                    val = []

                    # NOTE: Computing all observables from complex Coherent Flux
                    if arrType == "OI_VIS2":
                        val.append(corrFlux2Vis2(vcompli))
                        quantities.append("VIS2DATA")

                    elif arrType == "OI_VIS":
                        if dataType & oimDataType.VISAMP_ABS:
                            val.append(corrFlux2VisAmpAbs(vcompli))
                            quantities.append("VISAMP")
                        elif dataType & oimDataType.VISAMP_DIF:
                            val.append(corrFlux2VisAmpDif(vcompli))
                            quantities.append("VISAMP")
                        elif dataType & oimDataType.VISAMP_COR:
                            val.append(corrFlux2VisAmpCor(vcompli))
                            quantities.append("VISAMP")

                        if dataType & oimDataType.VISPHI_ABS:
                            val.append(corrFlux2VisPhiAbs(vcompli))
                            quantities.append("VISPHI")
                        elif dataType & oimDataType.VISPHI_DIF:
                            val.append(corrFlux2VisPhiDif(vcompli))
                            quantities.append("VISPHI")

                    elif arrType == "OI_T3":
                        if dataType & oimDataType.T3AMP:
                            val.append(corrFlux2T3Amp(vcompli))
                            quantities.append("T3AMP")
                        if dataType & oimDataType.T3PHI:
                            val.append(corrFlux2T3Phi(vcompli))
                            quantities.append("T3PHI")

                    elif arrType == "OI_FLUX":
                        val.append(corrFlux2Flux(vcompli))
                        # Fucking GRAVITY patch!
                        if "FLUXDATA" in [
                            c.name
                            for c in data.data[ifile][arrNum].data.columns
                        ]:
                            quantities.append("FLUXDATA")
                        else:
                            quantities.append("FLUX")

                    # NOTE: Filling the simulatedData astropy array with the computed values
                    if computeSimulatedData:
                        for ival in range(len(val)):
                            try:
                                self.simulatedData.data[ifile][arrNum].data[
                                    quantities[ival]
                                ] = val[ival]
                            except:
                                self.simulatedData.data[ifile][arrNum].data[
                                    quantities[ival]
                                ] = np.squeeze(val[ival])

                    # NOTE: Computing the chi2
                    if computeChi2 == True:
                        for ival in range(len(val)):

                            if nwl == 1 and len(dataVal[ival].shape) == 1:
                                dataVal[ival] = dataVal[ival][:, None]
                                dataErr[ival] = dataErr[ival][:, None]
                                flag[ival] = flag[ival][:, None]

                            # NOTE: For phase quantities go to the complex plane
                            if quantities[ival] in dataTypes:
                                if quantities[ival] in ["VISPHI", "T3PHI"]:
                                    dphi = np.rad2deg(
                                        np.angle(
                                            np.exp(
                                                1j * np.deg2rad(dataVal[ival])
                                            )
                                            * np.exp(
                                                -1j * np.deg2rad(val[ival])
                                            )
                                        )
                                    )
                                    resi = (
                                        dphi
                                        * np.logical_not(flag[ival])
                                        / dataErr[ival]
                                    )
                                    chi2i = resi**2

                                else:
                                    resi = (
                                        (dataVal[ival] - val[ival])
                                        * np.logical_not(flag[ival])
                                        / dataErr[ival]
                                    )
                                    chi2i = resi**2

                                nelChi2 += np.sum(
                                    (dataErr[ival] != 0)
                                    * np.logical_not(flag[ival])
                                )
                                chi2 += np.sum(np.nan_to_num(chi2i, nan=0))

                                chi2List.append(chi2i)
                                residuals.append(resi)
        if computeChi2:
            self.chi2List = chi2List
            self.chi2_0 = chi2
            self.nelChi2_0 = nelChi2
            self.residuals = residuals
            self.degFree = nelChi2 - len(self.model.getFreeParameters())

            if nelChi2 != 0:
                self.chi2r_0 = self.chi2_0 / self.degFree
            else:
                self.chi2r_0 = 0

            if self.cprior:
                self.chi2Prior = self.cprior()

                if nelChi2 != 0:
                    self.chi2 = (
                        self.chi2_0
                        + self.priorWeight * self.chi2Prior * self.nelChi2_0
                    )
                    self.degFree = self.nelChi2_0 * (
                        1 + self.priorWeight
                    ) - len(self.model.getFreeParameters())
                else:
                    self.chi2 = self.chi2Prior * self.priorWeight
                    self.degFree = self.priorWeight - len(
                        self.model.getFreeParameters()
                    )
            else:
                self.chi2Prior = 0
                self.chi2 = self.chi2_0

            self.chi2r = self.chi2 / self.degFree
            self.nelChi2 = self.degFree

    def computeAll(self, checkSimulatedData=True, dataTypes=None, cprior=None):
        self.compute(
            computeChi2=True,
            computeSimulatedData=True,
            checkSimulatedData=checkSimulatedData,
            dataTypes=dataTypes,
            cprior=cprior,
        )

    def plotWlTemplate(
        self,
        shape,
        simulated: bool = True,
        savefig=None,
        xunit: str = "m",
        plotFuntionData=_errorplot,
        plotFunctionSimulatedData=plt.Axes.plot,
        kwargsData: dict = {},
        kwargsSimulatedData: dict = {},
        **kwargs,
    ):
        kwargsData0 = {"color": "tab:red", "alpha": 0.5}
        kwargsSimulatedData0 = {"color": "tab:blue", "lw": 2, "alpha": 0.7}
        kwargsData = {**kwargsData0, **kwargsData}
        kwargsSimulatedData = {**kwargsSimulatedData0, **kwargsSimulatedData}

        fig = plt.figure(FigureClass=oimWlTemplatePlots, **kwargs)
        fig.autoShape(self.data.data, shape=shape)
        fig.set_xunit(xunit)
        fig.plot(
            self.data.data,
            plotFunction=plotFuntionData,
            plotFunctionkwarg=kwargsData,
        )
        fig.plot(
            self.simulatedData.data,
            plotFunction=plotFunctionSimulatedData,
            plotFunctionkwarg=kwargsSimulatedData,
        )
        return fig

    def plot(
        self,
        arr: str | list[str],
        simulated: bool = True,
        savefig: str | Path | None = None,
        visLog: bool = False,
        xaxis: str = "SPAFREQ",
        xunit: str = "cycle/rad",
        cname: str = "EFF_WAVE",
        cunit: str = "micron",
        cmap: str = "plasma",
        colorbar: bool = True,
        kwargsData: dict = {},
        kwargsSimulatedData: dict = {},
        fig: Figure | None = None,
        axe: Axes | None = None,
    ) -> None:
        """Plots data vs. simulated data.

        Parameters
        ----------
        arr : str or list of str
            The name(s) of the OIFITS column(s) to be plotted.
        simulated : bool, optional
            If `True`, simulated data are plotted. Default is `True`.
        savefig : str or pathlib.Path, optional
            Saves the plot. Default is `None`.
        visLog : bool, optional
            If `True`, sets the y-scale to logarithmic. Default is `False`.
        xaxis : str, optional
            OIFITS information plotted on the x-axis. Default is `"SPAFREQ"`.
        xunit : str, optional
            Unit of the x-axis. Default is `"cycle/rad"`.
        cname : str, optional
            OIFITS information plotted on the colorbar. Default is `"EFF_WAVE"`.
        cunit : str, optional
            Unit of the colorbar. Default is `"micron"`.
        cmap : str, optional
            Name of the colormap. Default is `"plasma"`.
        colorbar : bool, optional
            If `True`, plots a colorbar. Default is `True`.
        kwargsData : dict, optional
            Data keyword arguments passed to the `oiplot` method. Default is `{}`.
        kwargsSimulatedData : dict, optional
            Simulated data keyword arguments passed to `oiplot` method. Default is `{}`.
        fig : matplotlib.figure.Figure, optional
            Figure used for the plotting. Default is `None`.
        axe : matplotlib.axes.Axes, optional
            Axes used for the plotting. Default is `None`.
        """
        kwargsData0 = {
            "cname": cname,
            "cunit": cunit,
            "lw": 2,
            "cmap": cmap,
            "errorbar": True,
            "label": "data",
        }
        kwargsSimulatedData0 = {
            "color": "k",
            "ls": ":",
            "lw": 1,
            "label": "model",
        }
        kwargsData = {**kwargsData0, **kwargsData}
        kwargsData["cunit"] = u.Unit(kwargsData["cunit"])

        if "color" in kwargsData:
            kwargsData.pop("cmap")
            kwargsData.pop("cname")
            kwargsData.pop("cunit")

        kwargsSimulatedData = {**kwargsSimulatedData0, **kwargsSimulatedData}

        if type(arr) != type([]):
            arr = [arr]

        if fig is None or axe is None:
            fig, axe = plt.subplots(
                len(arr),
                1,
                sharex=True,
                figsize=(8, 6),
                subplot_kw={"projection": "oimAxes"},
            )

        if len(arr) == 1:
            axe = np.array([axe])

        plt.subplots_adjust(left=0.09, top=0.98, right=0.98, hspace=0.14)

        for iax, axi in enumerate(axe):
            scale = axi.oiplot(
                self.data.data,
                xaxis,
                arr[iax],
                xunit=xunit,
                showColorbar=False,
                **kwargsData,
            )

            if simulated:
                axi.oiplot(
                    self.simulatedData.data,
                    xaxis,
                    arr[iax],
                    xunit=xunit,
                    **kwargsSimulatedData,
                )

            if axi != axe[-1]:
                axi.get_xaxis().set_visible(False)
            if axi == axe[0]:
                axi.legend()

            if arr[iax] in ["VIS2DATA", "VISAMP"] and visLog:
                axi.set_yscale("log")

            axi.autolim()
            axi.margins(x=0)

        if colorbar:
            idxC = np.where(oimPlotParamName == kwargsData["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargsData['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=axe.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig is not None:
            plt.savefig(savefig)

        return fig, axe

    # TODO: Make it possible to select if the residuals should
    # be divided by the error or not? Perhaps with another kwarg?
    def plotResiduals(
        self,
        arr,
        xaxis: str = "SPAFREQ",
        xunit: str = "cycle/rad",
        savefig: str | Path | None = None,
        visLog: bool = False,
        cname: str = "EFF_WAVE",
        cunit: str = "micron",
        cmap: str = "plasma",
        colorbar: bool = True,
        marker: str = ".",
        levels: list[int] | None = [1, 2, 3],
        fig: Figure | None = None,
        axe: Axes | None = None,
        **kwargs,
    ) -> tuple[Figure, Axes]:
        """Plots the residuals computed from the data subtracted by the simulatedData
        divided by the errors on the data.

        Parameters
        ----------
        arr : str or list of str
            The name(s) of the OIFITS column(s) to be plotted.
        xaxis : str, optional
            OIFITS information plotted on the x-axis. Default is `"SPAFREQ"`.
        xunit : str, optional
            Unit of the x-axis. Default is `"cycle/rad"`.
        savefig : str or pathlib.Path, optional
            Saves the plot. Default is `None`.
        visLog : bool, optional
            If `True`, sets the y-scale to logarithmic. Default is `False`.
        cname : str, optional
            OIFITS information plotted on the colorbar. Default is `"EFF_WAVE"`.
        cunit : str, optional
            Unit of the colorbar. Default is `"micron"`.
        cmap : str, optional
            Name of the colormap. Default is `"plasma"`.
        colorbar : bool, optional
            If `True`, plots a colorbar. Default is `True`.
        marker : str, optional
            The marker of the residuals. Default is `"."`.
        levels : list of int, optional
            Marks residual levels with horizontal lines. Default is `[1, 2, 3]`.
        fig : matplotlib.figure.Figure, optional
            Figure used for the plotting. Default is `None`.
        axe : matplotlib.axes.Axes, optional
            Axes used for the plotting. Default is `None`.
        """
        kwargs = {
            "cname": cname,
            "cunit": cunit,
            "cmap": cmap,
            **kwargs,
        }
        kwargs["cunit"] = u.Unit(kwargs["cunit"])

        if "color" in kwargs:
            kwargs.pop("cmap")
            kwargs.pop("cname")
            kwargs.pop("cunit")

        if type(arr) != type([]):
            arr = [arr]

        residual_data = oimData()
        for i, (dat, fit_dat) in enumerate(
            zip(self.data.data, self.simulatedData.data)
        ):
            residual_data.addData(hdulistDeepCopy(fit_dat))
            for param in arr:
                idx_p = np.where(oimPlotParamName == param)[0][0]
                p_arr = oimPlotParamArr[idx_p]
                p_err = oimPlotParamError[idx_p]

                if p_arr not in dat:
                    continue

                if param in ["T3PHI", "VISPHI"]:
                    res_ph = (
                        np.rad2deg(
                            np.angle(
                                np.exp(
                                    1j * np.deg2rad(dat[p_arr].data[param])
                                    - 1j
                                    * np.deg2rad(fit_dat[p_arr].data[param])
                                )
                            )
                        )
                    ) / dat[p_arr].data[p_err]
                    residual_data.data[i][p_arr].data[param] = res_ph
                else:
                    res_vis = (
                        dat[p_arr].data[param] - fit_dat[p_arr].data[param]
                    ) / fit_dat[p_arr].data[p_err]
                    residual_data.data[i][p_arr].data[param] = res_vis

        if fig is None or axe is None:
            fig, axe = plt.subplots(
                len(arr),
                1,
                subplot_kw={"projection": "oimAxes"},
                figsize=(14, 8),
                constrained_layout=True,
            )

        if len(arr) == 1:
            axe = np.array([axe])

        for axi, arri in zip(axe, arr):
            idx_p = np.where(oimPlotParamName == arri)[0][0]
            label_p = oimPlotParamLabel[idx_p]

            scale = axi.oiplot(
                residual_data,
                xaxis,
                arri,
                xunit=xunit,
                marker=marker,
                showColorbar=False,
                **kwargs,
            )

            if levels is not None:
                alpha = np.linspace(1, 0.2, num=len(levels))
                axi.axhline(0, ls="-", color="grey")
                for j, leveli in enumerate(levels):
                    for k in range(2):
                        axi.axhline(
                            (2 * k - 1) * leveli,
                            ls="--",
                            color="grey",
                            alpha=alpha[j],
                        )

            axi.set_ylabel(f"{label_p} Residuals")
            axi.margins(x=0)

            if axi != axe[-1]:
                axi.get_xaxis().set_visible(False)

        if colorbar:
            idxC = np.where(oimPlotParamName == kwargs["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargs['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=axe.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig is not None:
            plt.savefig(savefig)

        return fig, axe

    # TODO: Make this a combination of the "plot" and "plotResiduals" methods?
    def plotWithResiduals(
        self,
        arr,
        simulated: bool = True,
        savefig: str | Path | None = None,
        visLog: bool = False,
        xaxis: str = "SPAFREQ",
        xunit: str = "cycle/rad",
        cname: str = "EFF_WAVE",
        cunit: str = "micron",
        cmap: str = "plasma",
        colorbar: bool = True,
        kwargsData: dict = {},
        kwargsSimulatedData: dict = {},
        kwargsResiduals: dict = {},
        levels: list[int] = [1, 2, 3],
        fig: Figure | None = None,
        axe: Axes | None = None,
    ) -> tuple[Figure, Axes]:
        """Plots the data and in a seperate plot underneath, the residuals
        computed from the data subtracted by the simulatedData.

        Parameters
        ----------
        arr : str or list of str
            The name(s) of the OIFITS column(s) to be plotted.
        xaxis : str, optional
            OIFITS information plotted on the x-axis. Default is `"SPAFREQ"`.
        xunit : str, optional
            Unit of the x-axis. Default is `"cycle/rad"`.
        savefig : str or pathlib.Path, optional
            Saves the plot. Default is `None`.
        visLog : bool, optional
            If `True`, sets the y-scale to logarithmic. Default is `False`.
        cname : str, optional
            OIFITS information plotted on the colorbar. Default is `"EFF_WAVE"`.
        cunit : str, optional
            Unit of the colorbar. Default is `"micron"`.
        cmap : str, optional
            Name of the colormap. Default is `"plasma"`.
        colorbar : bool, optional
            If `True`, plots a colorbar. Default is `True`.
        marker : str, optional
            The marker of the residuals. Default is `"."`.
        levels : list of int, optional
            Marks residual levels with horizontal lines. Default is `[1, 2, 3]`.
        fig : matplotlib.figure.Figure, optional
            Figure used for the plotting. Default is `None`.
        axe : matplotlib.axes.Axes, optional
            Axes used for the plotting. Default is `None`.
        """

        # NOTE: Plotting  data and simulated data
        kwargsData0 = {
            "cname": cname,
            "cunit": cunit,
            "lw": 2,
            "cmap": cmap,
            "errorbar": True,
            "label": "data",
        }
        kwargsData = {**kwargsData0, **kwargsData}
        kwargsData["cunit"] = u.Unit(kwargsData["cunit"])

        kwargsSimulatedData0 = {
            "color": "k",
            "ls": ":",
            "lw": 1,
            "label": "model",
        }
        kwargsSimulatedData = {**kwargsSimulatedData0, **kwargsSimulatedData}
        kwargsResiduals0 = {"cname": cname, "cunit": cunit, "cmap": cmap}

        kwargsResiduals = {**kwargsResiduals0, **kwargsResiduals}
        kwargsResiduals["cunit"] = u.Unit(cunit)
        if not ("ls" in kwargsResiduals) and not (
            "linestyle" in kwargsResiduals
        ):
            kwargsResiduals["ls"] = ""
            kwargsResiduals["marker"] = "."

        if "color" in kwargsData:
            kwargsData.pop("cmap")
            kwargsData.pop("cname")
            kwargsData.pop("cunit")
            kwargsResiduals.pop("cmap")
            kwargsResiduals.pop("cname")
            kwargsResiduals.pop("cunit")
            kwargsResiduals["color"] = kwargsData["color"]

        if isinstance(arr, str) or not isinstance(arr, Iterable):
            arr = [arr]

        residuals_data = oimData()
        for i, (dat, fit_dat) in enumerate(
            zip(self.data.data, self.simulatedData.data)
        ):
            residuals_data.addData(hdulistDeepCopy(fit_dat))
            for param in arr:

                idx_p = np.where(oimPlotParamName == param)[0][0]
                p_arr = oimPlotParamArr[idx_p]
                p_err = oimPlotParamError[idx_p]

                if p_arr not in dat:
                    continue

                if param in ["T3PHI", "VISPHI"]:
                    res_ph = (
                        np.rad2deg(
                            np.angle(
                                np.exp(
                                    1j * np.deg2rad(dat[p_arr].data[param])
                                    - 1j
                                    * np.deg2rad(fit_dat[p_arr].data[param])
                                )
                            )
                        )
                    ) / dat[p_arr].data[p_err]
                    residuals_data.data[i][p_arr].data[param] = res_ph
                else:
                    res_vis = (
                        dat[p_arr].data[param] - fit_dat[p_arr].data[param]
                    ) / fit_dat[p_arr].data[p_err]
                    residuals_data.data[i][p_arr].data[param] = res_vis

        # NOTE: Set the projection to oimAxes for all subplots to use oimodeler
        # custom plots
        nplots = len(arr)
        height_ratios = np.arange(1, nplots * 2 + 1) % 2 * 2 + 1

        if fig is None or axe is None:
            fig, axe = plt.subplots(
                len(arr) * 2,
                1,
                sharex=True,
                figsize=(8, 6),
                subplot_kw={"projection": "oimAxes"},
                height_ratios=height_ratios,
            )
        else:
            axe = np.array(axe)

        plt.subplots_adjust(left=0.09, top=0.98, right=0.98, hspace=0.14)

        # NOTE: Plotting loop: Plotting data and simulated data for each data type in arr
        for i in range(nplots):
            # NOTE: Plotting the data with wavelength colorscale + errorbars vs
            # spatial frequencies
            scale = axe[2 * i].oiplot(
                self.data.data,
                xaxis,
                arr[i],
                xunit=xunit,
                showColorbar=False,
                **kwargsData,
            )

            # NOTE: Over-plotting the simulated data as a dotted line vs spatial
            # frequencies
            axe[2 * i].oiplot(
                self.simulatedData.data,
                xaxis,
                arr[i],
                xunit=xunit,
                **kwargsSimulatedData,
            )

            if axe[2 * i] != axe[-1]:
                axe[2 * i].get_xaxis().set_visible(False)
            if axe[2 * i] == axe[0]:
                axe[2 * i].legend()

            # NOTE: Automatic ylim => 0-1 for visibilties, -180,180 for phases
            if arr[i] in ["VIS2DATA", "VISAMP"] and visLog == True:
                axe[2 * i].set_yscale("log")

            axe[2 * i].autolim()
            axe[2 * i].margins(x=0)

        # NOTE: Plotting loop: Plotting residuals
        for i in range(nplots):

            scale = axe[2 * i + 1].oiplot(
                residuals_data,
                xaxis,
                arr[i],
                xunit=xunit,
                showColorbar=False,
                **kwargsResiduals,
            )

            if levels is not None:
                alpha = np.linspace(1, 0.2, num=len(levels))
                axe[2 * i + 1].axhline(0, ls="-", color="grey")
                for j, leveli in enumerate(levels):
                    for k in range(2):
                        axe[2 * i + 1].axhline(
                            (2 * k - 1) * leveli,
                            ls="--",
                            color="grey",
                            alpha=alpha[j],
                        )

            axe[2 * i + 1].set_ylabel(r"$(\sigma)$")
            axe[2 * i + 1].margins(x=0)

            ymax = np.max(np.abs(axe[2 * i + 1].get_ylim()))
            axe[2 * i + 1].set_ylim(-ymax, ymax)

            if axe[i] != axe[-1]:
                axe[i].get_xaxis().set_visible(False)

        # NOTE: Create a colorbar for the data plotted with wavelength colorscale option
       
        if colorbar and cname in kwargsData:
            idxC = np.where(oimPlotParamName == kwargsData["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargsData['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=axe.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig != None:
            plt.savefig(savefig)

        return fig, axe
