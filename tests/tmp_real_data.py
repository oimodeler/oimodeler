# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:53:04 2022

@author: Ame
"""
from datetime import datetime
from pathlib import Path
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


matplotlib.use('Agg')
name = "realData"

path = Path(__file__).parent.parent.parent

fc1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
fc2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
fcwl = oim.oimWavelengthRangeFilter(targets="all", wlRange=[3.1e-6, 3.9e-6])
filt1 = oim.oimDataFilter([fc1, fc2, fcwl])


# %%
class baseTest:
    name = "basetest"
    truth = []
    chi2r0 = 1
    nwalker = 10
    nstep = 1000
    pathData = None
    fdata = None
    filt = filt1

    def __init__(self):
        self.setModel()
        self.data = oim.oimData((self.pathData / self.fdata).resolve())
        self.data.setFilter(self.filt)

    def setModel(self):
        self.model = None

    def compute(self, progress=False):
        t0 = datetime.now()

        self.fit = oim.oimFitterEmcee(
            self.data, self.model, nwalkers=self.nwalker)
        self.fit.prepare(init="random")
        self.fit.run(nsteps=self.nstep, progress=progress)

        discard = np.max([self.nstep-1000, self.nstep//2])
        self.best, _, _, self.err = self.fit.getResults(
            mode='best', discard=discard)
        self.chi2r = self.fit.simulator.chi2r

        resTest1 = (self.best-self.truth)/self.err < 1
        resTest2 = (self.best-self.truth)/self.best < 0.01
        resTest = np.all(resTest1 | resTest2)
        chi2Test = ((self.chi2r-self.chi2r0)/self.chi2r0 < 0.1)

        dt = (datetime.now()-t0).total_seconds()

        return resTest, chi2Test, dt

    def makePlots(self, prefix="", ext="pdf"):

        fname = prefix+self.name.replace(" ", "_")

        figWalkers, _ = self.fit.walkersPlot(savefig=path / "examples" / "testSuite" / "RESULTS" / f"{fname}_walkers{ext}")
        plt.close(figWalkers)

        figCorner, _ = self.fit.cornerPlot(discard=self.nstep//2,
                                                   savefig=path / "examples" / "testSuite" / "RESULTS" / f"{fname}_corner{ext}")
        plt.close(figCorner)

        figSim, _ = self.fit.simulator.plot(["VIS2DATA", "T3PHI"],
                                                savefig=path / "examples" / "testSuite" / "RESULTS" / f"{fname}_v2CP{ext}")
        plt.close(figSim)


# %%
class matisse75Vir(baseTest):
    name = "MATISSE Binary 75 Vir "
    # PIONIER data generated with ASPRO have somme bias
    truth = np.array([0.45, 4.5, -33.9, -79.3])
    chi2r0 = 16.  # and Chi2r is not close to 1
    nwalker = 50
    nstep = 30000
    pathData = pathData = path / "examples" / "data" / "RealData" / "MATISSE" / "binary75Vir"
    fdata = "2019-05-23T025507_75Vir_A0G1J2J3_IR-LM_LOW_noChop_cal_oifits_0.fits"
    filt = filt1

    def setModel(self):
        ud1 = oim.oimUD()
        pt2 = oim.oimPt()
        pt2.params["x"].set(min=-200, max=200, free=True)
        pt2.params["y"].set(min=-200, max=200, free=True)
        ud1.params["f"].set(free=True)
        ud1.params["d"].set(min=0, max=10)
        pt2.params["f"] = oim.oimParamNorm(ud1.params["f"])
        self.model = oim.oimModel(ud1, pt2)


# %%
tests = [matisse75Vir()]
ntest = len(tests)
for itest, testi in enumerate(tests):
    if itest >= 0:
        res, chi2, dt = testi.compute(progress=True)
        pprint("{}/{}) {} => res={} chi2={} (dt={:.1f}s)".format(itest + 1, ntest, testi.name, res, chi2, dt))
        testi.makePlots(prefix=name+"_"+str(itest+1)+"_")
