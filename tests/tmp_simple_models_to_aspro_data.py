# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:53:04 2022

@author: Ame
"""
# TODO: Finish this test in pytest format
from datetime import datetime
from pathlib import Path
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


matplotlib.use('Agg')
name = "simpleASPRO"

path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "SIMPLE_TESTS"

fc1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
fc2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
filt1 = oim.oimDataFilter([fc1, fc2])


# %%
class baseTest:
    name = "basetest"
    truth = []
    chi2r0 = 1
    nwalker = 10
    nstep = 1000
    fdata = None
    filt = filt1

    def __init__(self):
        self.setModel()
        self.data = oim.oimData((data_dir / self.fdata).resolve())
        self.data.setFilter(self.filt)

    def setModel(self):
        self.model = None

    def compute(self, progress=False):
        t0 = datetime.now()

        self.fit = oim.oimFitterEmcee(
            self.data, self.model, nwalkers=self.nwalker)
        self.fit.prepare(init="random")
        self.fit.run(nsteps=self.nstep, progress=progress)

        self.best, _, _, self.err = self.fit.getResults(
            mode='median', discard=self.nstep//2)
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
class pionierUD(baseTest):
    name = "PIONIER UD"
    # PIONIER data generated with ASPRO have somme bias
    truth = np.array([10.014])
    chi2r0 = 2.934  # and Chi2r is not close to 1
    nwalker = 10
    nstep = 1000
    fdata = "ASPRO_PIONIER_UD10.fits"
    filt = filt1

    def setModel(self):
        c1 = oim.oimUD()
        self.model = oim.oimModel(c1)
        c1.params['d'].set(min=0, max=20)
        c1.params['f'].set(free=False)


class pionierBin(baseTest):
    name = "PIONIER PT + SHIFT UD"
    # PIONIER data generated with ASPRO have somme bias
    truth = np.array([0.83, 5, 8.67, 3])
    chi2r0 = 3.973  # and Chi2r is not close to 1
    nwalker = 10
    nstep = 2000
    fdata = "ASPRO_PIONIER_1PT_0.2UD3_sep3_pa30.fits"
    filt = filt1

    def setModel(self):
        c1 = oim.oimPt()
        c2 = oim.oimUD()
        self.model = oim.oimModel(c1, c2)
        c1.params['f'].set(min=0, max=1)
        c2.params['x'].set(free=True, min=-10, max=10)
        c2.params['y'].set(free=True, min=-10, max=10)
        c2.params['d'].set(max=10)
        c2.params['f'] = oim.oimParamNorm(c1.params['f'])


class gravityUDAndShiftedGauss(baseTest):
    name = "GRAVITY UD + SHIFT EGAUSS"
    truth = np.array([0.2, 2, 0.867, 0.5, 2, 80, 10])
    nwalker = 20
    nstep = 2000
    fdata = "ASPRO_GRAVITY_0.2UD2_0.8GE10_apl2_pa80_sep1_pa60.fits"

    def setModel(self):
        c1 = oim.oimUD()
        c2 = oim.oimEGauss()
        self.model = oim.oimModel(c1, c2)
        c1.params['f'].set(min=0, max=1)
        c1.params['d'].set(min=0, max=10)
        c2.params['x'].set(free=True, min=-10, max=10)
        c2.params['y'].set(free=True, min=-10, max=10)
        c2.params['pa'].set(min=0, max=180)
        c2.params['fwhm'].set(max=50)
        c2.params['elong'].set(max=5)
        c2.params['f'] = oim.oimParamNorm(c1.params['f'])


class matissePtAndEIRing(baseTest):
    name = "MATISSE PT + EIRING"
    truth = np.array([0.3, 1.5, 50, 15])
    nwalker = 10
    nstep = 4000
    fdata = "ASPRO_MATISSE_0.3PT_0.7ER15_apl1.5_pa=50.fits"

    def setModel(self):
        c1 = oim.oimPt()
        c2 = oim.oimEIRing()
        self.model = oim.oimModel(c1, c2)
        c1.params['f'].set(min=0, max=1)
        c2.params['pa'].set(min=0, max=180)
        c2.params['d'].set(max=50)
        c2.params['elong'].set(max=5)
        c2.params['f'] = oim.oimParamNorm(c1.params['f'])


class matissePtAndERing(baseTest):
    name = "MATISSE PT + ERING"
    truth = np.array([0.3, 1.5, 50, 15, 0])
    nwalker = 20
    nstep = 4000
    chi2r0 = 1.17
    fdata = "ASPRO_MATISSE_0.3PT_0.7ER15_apl1.5_pa=50.fits"

    def setModel(self):
        c1 = oim.oimPt()
        c2 = oim.oimERing2()
        self.model = oim.oimModel(c1, c2)
        c1.params['f'].set(min=0, max=1)
        c2.params['pa'].set(min=0, max=180)
        c2.params['d'].set(max=50)
        c2.params['w'].set(max=10)
        c2.params['elong'].set(max=5)
        c2.params['f'] = oim.oimParamNorm(c1.params['f'])


# %%
tests = [pionierUD(), pionierBin(), gravityUDAndShiftedGauss(),
         matissePtAndEIRing(), matissePtAndERing()]
ntest = len(tests)
for itest, testi in enumerate(tests):
    if itest >= 0:
        res, chi2, dt = testi.compute(progress=True)
        pprint("{}/{}) {} => res={} chi2={} (dt={:.1f}s)".format(itest + 1, ntest, testi.name, res, chi2, dt))
        testi.makePlots(prefix=name+"_"+str(itest+1)+"_")
