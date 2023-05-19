# -*- coding: utf-8 -*-
"""Data/model simulation"""
import matplotlib.pyplot as plt
import numpy as np

from .oimData import oimData, oimDataType
from .oimUtils import hdulistDeepCopy


def corrFlux2Vis2(vcompl):
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB-1), vcompl[0, :])
    return np.abs(vcompl[1:, :]/norm)**2


def corrFlux2VisAmpAbs(vcompl):
    nB = vcompl.shape[0]
    norm = np.outer(np.ones(nB-1), vcompl[0, :])
    return np.abs(vcompl[1:, :]/norm)


# FIXME : Not real formula for differential visibilities
def corrFlux2VisAmpDif(vcompl):
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    return np.abs(vcompl[1:, :]/norm)


def corrFlux2VisAmpCor(vcompl):
    return np.abs(vcompl[1:, :])


def corrFlux2VisPhiAbs(vcompl):
    return np.rad2deg(np.angle(vcompl[1:, :]))


# FIXME : Not real formula for differential phases
def corrFlux2VisPhiDif(vcompl):
    nlam = vcompl.shape[1]
    norm = np.outer(np.mean(vcompl[1:, :], axis=1), np.ones(nlam))
    phi = np.rad2deg(np.angle(vcompl[1:, :]*np.conjugate(norm)))
    return phi


# TODO: Special function doing T3Amp and T3Phi simultaneously
def corrFlux2T3Amp(vcompl):
    nB = vcompl.shape[0]
    nCP = (nB-1)//3
    norm = np.outer(np.ones(nCP), vcompl[0, :])
    BS = vcompl[1:nCP+1, :]*vcompl[nCP+1:2*nCP+1, :] * \
        np.conjugate(vcompl[2*nCP+1:, :])/norm**3
    return np.abs(BS)


def corrFlux2T3Phi(vcompl):
    nB = vcompl.shape[0]
    nCP = (nB-1)//3
    norm = np.outer(np.ones(nCP), vcompl[0, :])
    BS = vcompl[1:nCP+1, :]*vcompl[nCP+1:2*nCP+1, :] * \
        np.conjugate(vcompl[2*nCP+1:, :])/norm**3
    return np.rad2deg(np.angle(BS))


def corrFlux2Flux(vcompl):
    return np.abs(vcompl)


class oimSimulator(object):
    """Contains"""

    def __init__(self, data=None, model=None, fitter=None, **kwargs):
        self.data = oimData()
        self.simulatedData = None
        self.model = None

        if data != None:
            if isinstance(data, oimData):
                self.data = data
            else:
                self.addData(data)

        if model != None:
            self.setModel(model)

        if model != None and not(data is None):
            self.compute(computeChi2=True, computeSimulatedData=True)

    def setModel(self, model):
        self.model = model

    def addData(self, data):
        self.data.addData(data)

    def prepareData(self):
        self.data.prepareData()
        self.simulatedData = oimData()
        for datai in self.data.data:
            self.simulatedData.addData(hdulistDeepCopy(datai))

    def compute(self, computeChi2=False, computeSimulatedData=False, checkSimulatedData=True):
        self.vcompl = self.model.getComplexCoherentFlux(self.data.vect_u, self.data.vect_v, self.data.vect_wl, self.data.vect_mjd)

        nelChi2 = 0
        chi2 = 0
        chi2List = []

        if computeSimulatedData == True and (checkSimulatedData == True or self.simulatedData == None):
            self.simulatedData = oimData()
            for datai in self.data.data:
                self.simulatedData.addData(hdulistDeepCopy(datai))

        data = self.data

        if (computeChi2 == True) | (computeSimulatedData == True):

            idx = 0
            nfiles = len(data.struct_u)
            for ifile in range(nfiles):
                # print("Data {}".format(ifile))
                narr = len(data.struct_arrType[ifile])
                for iarr in range(narr):
                    arrNum = data.struct_arrNum[ifile][iarr]
                    arrType = data.struct_arrType[ifile][iarr]
                    dataType = data.struct_dataType[ifile][iarr]
                    nB = data.struct_nB[ifile][iarr]
                    nwl = data.struct_nwl[ifile][iarr]
                    vcompli = self.vcompl[idx:idx+nB*nwl]
                    vcompli = np.reshape(vcompli, [nB, nwl])

                    dataVal = data.struct_val[ifile][iarr]
                    dataErr = data.struct_err[ifile][iarr]
                    flag = data.struct_flag[ifile][iarr]

                    idx += nB*nwl
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
                        quantities.append("FLUXDATA")

                    # NOTE: Filling the simulatedData astropy array with the computed values
                    if computeSimulatedData == True:
                        for ival in range(len(val)):
                            try:
                                self.simulatedData.data[ifile][arrNum].data[quantities[ival]] = val[ival]
                            except:
                                self.simulatedData.data[ifile][arrNum].data[quantities[ival]] = np.squeeze(
                                    val[ival])

                    # NOTE: Computing the chi2
                    if computeChi2 == True:
                        for ival in range(len(val)):
                            # NOTE: For phase quantities go to the complex plane
                            if quantities[ival] in ["VISPHI", "T3PHI"]:
                                dphi = np.rad2deg(np.angle(np.exp(1j*np.deg2rad(dataVal[ival]))
                                                           * np.exp(-1j*np.deg2rad((val[ival])))))
                                chi2i = (dphi*np.logical_not(flag[ival])/dataErr[ival])**2

                            else:
                                chi2i = ((dataVal[ival]-val[ival])*np.logical_not(flag[ival])/dataErr[ival])**2

                            nelChi2 += np.sum((dataErr[ival] != 0)
                                              * np.logical_not(flag[ival]))
                            chi2 += np.sum(np.nan_to_num(chi2i, nan=0))
                            chi2List.append(chi2i)
    
       
        if computeChi2==True: 
             self.chi2=chi2
             self.chi2r=chi2/(nelChi2-len(self.model.getFreeParameters()))
             self.chi2List=chi2List
             self.nelChi2=nelChi2
                
           
    def plot(self,arr,simulated=True,savefig=None,visLog=False,**kwargs):
        # plotting  data and simulatiedData

        if type(arr)!=type([]):
           
            arr=[arr]
        
        #Set the projection to oimAxes for all subplots to use oimodeler custom plots
        fig,ax=plt.subplots(len(arr),1,sharex=True,figsize=(8,6),
                            subplot_kw=dict(projection='oimAxes'))

        if len(arr)==1:
            ax=np.array([ax])

        plt.subplots_adjust(left=0.09,top=0.98,right=0.98,hspace=0.14)

        # Ploting loop :  plotting data and simulated data for each data type in arr
        for iax,axi in enumerate(ax):
            
            #plotting the data with wavelength colorscale + errorbars vs spatial frequencies
            scale=axi.oiplot(self.data.data,"SPAFREQ",arr[iax] ,xunit="cycles/mas",
                    cname="EFF_WAVE",cunitmultiplier=1e6,lw=2,cmap="coolwarm",
                    errorbar=True,label="data")
            
            #over-plotting the simulated data as a dotted line  vs spatial frequencies
            axi.oiplot(self.simulatedData.data,"SPAFREQ",arr[iax] ,xunit="cycles/mas",
                    color="k",ls=":",lw=1,label="model")
            
            if axi!=ax[-1]: axi.get_xaxis().set_visible(False)
            if axi==ax[0]:axi.legend()
                    
            #automatic ylim => 0-1 for visibilties, -180,180 for phases

            if arr[iax] in ["VIS2DATA","VISAMP"] and visLog==True:
                axi.set_yscale("log")

            axi.autolim()

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
        fig.colorbar(scale, ax=ax.ravel().tolist(),
                     label="$\\lambda$ ($\mu$m)")

        if savefig != None:
            plt.savefig(savefig)

        return fig, ax
