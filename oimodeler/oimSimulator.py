# -*- coding: utf-8 -*-
"""Data/model simulation"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

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



oimDataArrDict=dict()

oimDataArrDict["OI_VIS2"]=dict(data=["VIS2DATA"],err=["VIS2ERR"])
oimDataArrDict["OI_VIS"]=dict(data=["VISAMP","VISPHI"],err=["VISAMPERR","VISPHIERR"])
oimDataArrDict["OI_T3"]=dict(data=["T3AMP","T3PHI"],err=["T3AMPERR","T3PHIERR"])
oimDataArrDict["OI_VIS2"]=dict(data=["VIS2DATA"],err=["VIS2ERR"])
oimDataArrDict["OI_FLUX"]=dict(data=["FLUXDATA"],err=["FLUXERR"])





def corrFlux2Vis2(vcompl):
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB - 1), vcompl[0, :])
    return np.abs(vcompl[1:, :] / norm) ** 2


def corrFlux2VisAmpAbs(vcompl):
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB - 1), vcompl[0, :])
    return np.abs(vcompl[1:, :] / norm)


# FIXME : Not real formula for differential visibilities
def corrFlux2VisAmpDif(vcompl):
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    return np.abs(vcompl[1:, :] / norm)


def corrFlux2VisAmpCor(vcompl):
    return np.abs(vcompl[1:, :])


def corrFlux2VisPhiAbs(vcompl):
    return np.angle(vcompl[1:, :], deg=True)


# FIXME : Not real formula for differential phases
def corrFlux2VisPhiDif(vcompl):
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    return np.angle(vcompl[1:, :] * np.conjugate(norm), deg=True)


# TODO: Special function doing T3Amp and T3Phi simultaneously
def corrFlux2T3Amp(vcompl):
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


def corrFlux2T3Phi(vcompl):
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


def corrFlux2Flux(vcompl):
    return np.abs(vcompl)


class oimSimulator:
    """Contains"""

    def __init__(
        self, data=None, model=None, fitter=None, cprior=None, **kwargs
    ):
        self.data = oimData()
        self.simulatedData = None
        self.model = None
        self.cprior = cprior

        if data != None:
            if isinstance(data, oimData):
                self.data = data
            else:
                self.addData(data)

        if model != None:
            self.setModel(model)

        if model != None and not (data is None):
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
                    for dataName,errName in \
                            zip(oimDataArrDict[dataij.name]["data"], 
                                oimDataArrDict[dataij.name]["err"]):
                            
                            shape = dataij.data[dataName].shape
                            err = dataij.data[errName]
                            
                            
                            dataij.data[dataName]+=np.random.randn(*shape)*err
                
            self.bootstrapData.addData(dataic)   

    def compute(
        self,
        computeChi2=False,
        computeSimulatedData=False,
        checkSimulatedData=True,
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

        if (computeChi2 == True) | (computeSimulatedData == True) :
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
                                                -1j * np.deg2rad((val[ival]))
                                            )
                                        )
                                    )
                                    resi  = (
                                        dphi
                                        * np.logical_not(flag[ival])
                                        / dataErr[ival]
                                    )
                                    chi2i = resi**2

                                else:
                                    resi  = (
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

        if computeChi2 and self.cprior is None:
            self.chi2 = chi2
            self.chi2r = chi2 / (nelChi2 - len(self.model.getFreeParameters()))
            self.chi2List = chi2List
            self.nelChi2 = nelChi2
            self.residuals = residuals
        elif computeChi2 and nelChi2==0:
            chi2_prior=self.cprior(self.model.getParameters())

            self.chi2 = chi2_prior
            self.chi2r = chi2_prior
            self.chi2List = None
            self.nelChi2 = 1
            self.residuals = residuals
        elif computeChi2:
            chi2_prior = (
                chi2 + self.cprior(self.model.getParameters()) * nelChi2
            )
            self.chi2 = chi2_prior
            self.chi2r = chi2_prior / (
                nelChi2 - len(self.model.getFreeParameters())
            )

            self.chi2_np = chi2
            self.chi2r_np = chi2 / (
                nelChi2 - len(self.model.getFreeParameters())
            )

            self.chi2List = chi2List
            self.nelChi2 = nelChi2
            self.residuals = residuals

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
        simulated=True,
        savefig=None,
        xunit="m",
        plotFuntionData=_errorplot,
        plotFunctionSimulatedData=plt.Axes.plot,
        kwargsData={},
        kwargsSimulatedData={},
        **kwargs,
    ):
        kwargsData0 = dict(color="tab:red", alpha=0.5)
        kwargsSimulatedData0 = dict(color="tab:blue", lw=2, alpha=0.7)
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

        # fig.set_legends("$LENGTH$m $PA$$^o$","VIS2DATA",fontsize=8)
        # fig.set_legends(0.1,0.8,"$BASELINE$",["VIS2DATA","VISPHI","T3PHI"],fontweight =1000)

        return fig

    def plot(
        self,
        arr,
        simulated=True,
        savefig=None,
        visLog=False,
        xaxis="SPAFREQ",
        xunit="cycle/rad",
        cname="EFF_WAVE",
        cunit="micron",
        cmap="plasma",
        kwargsData={},
        kwargsSimulatedData={},
    ):
        # NOTE: Plotting  data and simulated data
        kwargsData0 = dict(
            cname=cname,
            cunit=cunit,
            lw=2,
            cmap=cmap,
            errorbar=True,
            label="data",
        )

        kwargsSimulatedData0 = dict(color="k", ls=":", lw=1, label="model")

        kwargsData = {**kwargsData0, **kwargsData}
        kwargsData["cunit"] = u.Unit(kwargsData["cunit"])

        if "color" in kwargsData:
            kwargsData.pop("cmap")
            kwargsData.pop("cname")
            kwargsData.pop("cunit")
        kwargsSimulatedData = {**kwargsSimulatedData0, **kwargsSimulatedData}

        if type(arr) != type([]):
            arr = [arr]

        # NOTE: Set the projection to oimAxes for all subplots to use oimodeler
        # custom plots
        fig, ax = plt.subplots(
            len(arr),
            1,
            sharex=True,
            figsize=(8, 6),
            subplot_kw=dict(projection="oimAxes"),
        )

        if len(arr) == 1:
            ax = np.array([ax])

        plt.subplots_adjust(left=0.09, top=0.98, right=0.98, hspace=0.14)

        # NOTE: Plotting loop: Plotting data and simulated data for each data type in arr
        for iax, axi in enumerate(ax):
            # NOTE: Plotting the data with wavelength colorscale + errorbars vs
            # spatial frequencies
            scale = axi.oiplot(
                self.data.data,
                xaxis,
                arr[iax],
                xunit=xunit,
                showColorbar=False,
                **kwargsData,
            )

            # NOTE: Over-plotting the simulated data as a dotted line vs spatial
            # frequencies
            axi.oiplot(
                self.simulatedData.data,
                xaxis,
                arr[iax],
                xunit=xunit,
                **kwargsSimulatedData,
            )

            if axi != ax[-1]:
                axi.get_xaxis().set_visible(False)
            if axi == ax[0]:
                axi.legend()

            # NOTE: Automatic ylim => 0-1 for visibilties, -180,180 for phases
            if arr[iax] in ["VIS2DATA", "VISAMP"] and visLog == True:
                axi.set_yscale("log")

            axi.autolim()
            axi.margins(x=0)

        xmin = 1e99
        xmax = -1e99
        for axi in ax:
            for li in axi.get_lines():
                x = li.get_xdata()
                xmini = np.min(x)
                xmaxi = np.max(x)
                if xmini < xmin:
                    xmin = xmini
                if xmaxi > xmax:
                    xmax = xmaxi

        ax[0].set_xlim(xmin, xmax)

        # NOTE: Create a colorbar for the data plotted with wavelength colorscale option
        if "cname" in kwargsData:
            idxC = np.where(oimPlotParamName == kwargsData["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargsData['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=ax.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax

    def plotResiduals(
        self,
        arr,
        xaxis="SPAFREQ",
        xunit="cycle/rad",
        savefig=None,
        visLog=False,
        cname="EFF_WAVE",
        cunit="micron",
        cmap="plasma",
        marker=".",
        levels=[1, 2, 3],
        **kwargs,
    ):

        kwargs0 = dict(cname=cname, cunit=cunit, lw=2, cmap=cmap)

        kwargs = {**kwargs0, **kwargs}
        kwargs["cunit"] = u.Unit(kwargs["cunit"])

        if "color" in kwargs:
            kwargs.pop("cmap")
            kwargs.pop("cname")
            kwargs.pop("cunit")

        if not ("ls" in kwargs) and not ("linestyle" in kwargs):
            kwargs["ls"] = ""

        residuals_data = oimData()
        for i, (dat, fit_dat) in enumerate(
            zip(self.data.data, self.simulatedData.data)
        ):
            residuals_data.addData(hdulistDeepCopy(fit_dat))
            for param in arr:

                idx_p = np.where(oimPlotParamName == param)[0][0]
                p_arr = oimPlotParamArr[idx_p]
                p_err = oimPlotParamError[idx_p]

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

        # kwargs = dict(cname = cname, cunit=cunit)
        idx_xaxis = np.where(oimPlotParamName == xaxis)[0][0]
        label_xaxis = oimPlotParamLabel[idx_xaxis]

        fig, ax = plt.subplots(
            len(arr),
            1,
            subplot_kw=dict(projection="oimAxes"),
            figsize=(14, 8),
            constrained_layout=True,
        )

        for i in range(len(arr)):
            idx_p = np.where(oimPlotParamName == arr[i])[0][0]
            label_p = oimPlotParamLabel[idx_p]

            scale = ax[i].oiplot(
                residuals_data,
                xaxis,
                arr[i],
                xunit=xunit,
                marker=marker,
                showColorbar=False,
                **kwargs,
            )

            xlim = ax[0].get_xlim()

            if type(levels) != type(None):
                alpha = np.linspace(1, 0.2, num=len(levels))
                ax[i].plot(xlim, [0, 0], ls="-", color="grey")
                for j, leveli in enumerate(levels):
                    for k in range(2):
                        y = np.array([1, 1]) * (2 * k - 1) * (leveli)
                        ax[i].plot(
                            xlim, y, ls="--", color="grey", alpha=alpha[j]
                        )
            # ax[i].set_title(label_p)
            ax[i].set_xlabel(label_xaxis + " (" + xunit + ")")
            ax[i].set_ylabel(f"{label_p} Residuals")
            ax[i].margins(x=0)

            if ax[i] != ax[-1]:
                ax[i].get_xaxis().set_visible(False)
        # NOTE: Create a colorbar for the data plotted with wavelength colorscale option

        if "cname" in kwargs:
            idxC = np.where(oimPlotParamName == kwargs["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargs['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=ax.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax

    def plotWithResiduals(
        self,
        arr,
        simulated=True,
        savefig=None,
        visLog=False,
        xaxis="SPAFREQ",
        xunit="cycle/rad",
        cname="EFF_WAVE",
        cunit="micron",
        cmap="plasma",
        kwargsData={},
        kwargsSimulatedData={},
        kwargsResiduals={},
        levels=[1, 2, 3],
    ):

        # NOTE: Plotting  data and simulated data
        kwargsData0 = dict(
            cname=cname,
            cunit=cunit,
            lw=2,
            cmap=cmap,
            errorbar=True,
            label="data",
        )
        kwargsData = {**kwargsData0, **kwargsData}
        kwargsData["cunit"] = u.Unit(kwargsData["cunit"])

        kwargsSimulatedData0 = dict(color="k", ls=":", lw=1, label="model")
        kwargsSimulatedData = {**kwargsSimulatedData0, **kwargsSimulatedData}

        kwargsResiduals0 = dict(cname=cname, cunit=cunit, cmap=cmap)

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

        if type(arr) != type([]):
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
        fig, ax = plt.subplots(
            len(arr) * 2,
            1,
            sharex=True,
            figsize=(8, 6),
            subplot_kw=dict(projection="oimAxes"),
            height_ratios=height_ratios,
        )

        plt.subplots_adjust(left=0.09, top=0.98, right=0.98, hspace=0.14)

        # NOTE: Plotting loop: Plotting data and simulated data for each data type in arr
        for i in range(nplots):
            # NOTE: Plotting the data with wavelength colorscale + errorbars vs
            # spatial frequencies
            scale = ax[2 * i].oiplot(
                self.data.data,
                xaxis,
                arr[i],
                xunit=xunit,
                showColorbar=False,
                **kwargsData,
            )

            # NOTE: Over-plotting the simulated data as a dotted line vs spatial
            # frequencies
            ax[2 * i].oiplot(
                self.simulatedData.data,
                xaxis,
                arr[i],
                xunit=xunit,
                **kwargsSimulatedData,
            )

            if ax[2 * i] != ax[-1]:
                ax[2 * i].get_xaxis().set_visible(False)
            if ax[2 * i] == ax[0]:
                ax[2 * i].legend()

            # NOTE: Automatic ylim => 0-1 for visibilties, -180,180 for phases
            if arr[i] in ["VIS2DATA", "VISAMP"] and visLog == True:
                ax[2 * i].set_yscale("log")

            ax[2 * i].autolim()
            ax[2 * i].margins(x=0)

        # NOTE: Plotting loop: Plotting residuals
        for i in range(nplots):

            scale = ax[2 * i + 1].oiplot(
                residuals_data,
                xaxis,
                arr[i],
                xunit=xunit,
                showColorbar=False,
                **kwargsResiduals,
            )

            xlim = ax[0].get_xlim()

            if type(levels) != type(None):
                alpha = np.linspace(1, 0.2, num=len(levels))
                ax[2 * i + 1].plot(xlim, [0, 0], ls="-", color="grey")
                for j, leveli in enumerate(levels):
                    for k in range(2):
                        y = np.array([1, 1]) * (2 * k - 1) * (leveli)
                        ax[2 * i + 1].plot(
                            xlim, y, ls="--", color="grey", alpha=alpha[j]
                        )
            # ax[i].set_title(label_p)
            # ax[2*i+1].set_xlabel(label_xaxis+' ('+xunit+')')
            ax[2 * i + 1].set_ylabel("($\\sigma$)")
            ax[2 * i + 1].margins(x=0)

            ymax = np.max(np.abs(ax[2 * i + 1].get_ylim()))
            ax[2 * i + 1].set_ylim(-ymax, ymax)

            if ax[i] != ax[-1]:
                ax[i].get_xaxis().set_visible(False)

        xmin = 1e99
        xmax = -1e99
        for axi in ax:
            for li in axi.get_lines():
                x = li.get_xdata()
                xmini = np.min(x)
                xmaxi = np.max(x)
                if xmini < xmin:
                    xmin = xmini
                if xmaxi > xmax:
                    xmax = xmaxi

        ax[0].set_xlim(xmin, xmax)

        # NOTE: Create a colorbar for the data plotted with wavelength colorscale option
        if "cname" in kwargsData:
            idxC = np.where(oimPlotParamName == kwargsData["cname"])[0][0]
            xlabel = oimPlotParamLabelShort[idxC]
            cunittext = f"{kwargsData['cunit']:latex_inline}"
            fig.colorbar(
                scale, ax=ax.ravel().tolist(), label=f"{xlabel} ({cunittext})"
            )

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax
