# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""
import numpy as np
from astropy import units as units

from ..oimComponent import oimComponentRadialProfile
from ..oimParam import oimParam, _standardParameters

class oimExpRing(oimComponentRadialProfile):
    name = "Exponential Ring"
    shortname = "ExpR"

    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**(_standardParameters["d"]))
        self.params["fwhm"] = oimParam(**(_standardParameters["fwhm"]))
        self.params["dim"] = oimParam(**(_standardParameters["dim"]))

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None  # np.array([0.5,1])*1e-6
        # self._r = np.arange(0, self._dim)*pixSize

        # Finally call the _eval function that allow the parameters to be processed
        self._eval(**kwargs)

    def _radialProfileFunction(self, r, wl, t):

        r0 = self.params["d"](wl, t)/2
        fwhm = self.params["fwhm"](wl, t)
        I = np.nan_to_num((r > r0).astype(float) *
                          np.exp(-0.692*np.divide(r-r0, fwhm)), nan=0)
        return I

    @property
    def _r(self):
        if False:
            fwhm_max = np.max(self.params["fwhm"](self._wl, self._t))
            r0_max = np.max(self.params["d"](self._wl, self._t))/2
        else:
            fwhm_max = self.params["fwhm"](1e99)
            r0_max = self.params["d"](1e99)
        rmax = r0_max+8*fwhm_max
        return np.linspace(0, 1,  self.params["dim"].value)*rmax

    @_r.setter
    def _r(self, r):
        pass
