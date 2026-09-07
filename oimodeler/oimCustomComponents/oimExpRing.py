# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""

import numpy as np

from ..oimComponent import oimComponentRadialProfile
from ..oimParam import oimParam


class oimExpRing(oimComponentRadialProfile):
    name = "Exponential Ring"
    shortname = "ExpR"

    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(base="d")
        self.params["fwhm"] = oimParam(base="fwhm")
        self.params["dim"] = oimParam(base="dim")
        self._eval(**kwargs)

    @property
    def r(self):
        if False:
            fwhm_max = np.max(self.fwhm(self._wl, self._t))
            r0_max = np.max(self.d(self._wl, self._t)) / 2
        else:
            fwhm_max, r0_max = self.fwhm(1e99), self.d(1e99)

        rmax = r0_max + 8 * fwhm_max
        self._r = np.linspace(0, 1, self.dim.value) * rmax
        return self._r

    def _radialProfileFunction(self, r, wl, t):
        r0, fwhm = self.d(wl, t) / 2, self.fwhm(wl, t)
        return (r > r0) * np.exp(-0.692 * np.divide(r - r0, fwhm))
