# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:39:36 2021

@author: Ama
"""

import numpy as np
from astropy import units as units
from scipy.constants import c, h, k
from astorpy.constants import M_sun
import astropy.units as u

from ..oimComponent import oimComponentRadialProfile
from ..oimParam import oimParam


class oimTempGradient(oimComponentRadialProfile):
    """A ring defined by a radial temperature profile in r^q and a radial dust
    surface density profile in r^p

    Parameters
    ----------
    rin : float
        Inner radius of the disk [au]
    rout : float
        Outer radius of the disk [au]
    Tin : float
        Inner radius temperature [K]
    Mdust : float
        Mass of the dusty disk [M_sun]
    q : float
        Power-law exponent for the temperature profile
    p : float
        Power-law exponent for the dust surface density profile
    kappa_abs : float
        Dust mass absorption coefficient [cm2.g-1]
    dist : float
        Distance of the star [pc]

    Attributes
    ----------
    params : Dict[str, oimParam]
        Dictionary of parameters
    _t : np.array
        Array of time values
    _wl : np.array
        Array of wavelength values
    _r : np.array

    Methods
    -------
    _eval(**kwargs)
        Evaluate the model's parameters
    _radialProfileFunction(r, wl, t)
        Radial profile function
    """
    name = "Temerature Gradient"
    shortname = "TempGrad"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0,
                                      description="Inner radius of the disk [au]")
        self.params["rout"] = oimParam(name="rout", value=0,
                                       description="Outer radius of the disk [au]")
        self.params["Tin"] = oimParam(name="Tin", value=0,
                                      description="Inner radius temperature [K]")
        self.params["Mdust"] = oimParam(name="Mdust", value=0,
                                        description="Mass of the dusty disk [M_sun]")
        self.params["q"] = oimParam(name="q", value=0,
                                    description="Power-law exponent for the temperature profile")
        self.params["p"] = oimParam(name="p", value=0,
                                    description="Power-law exponent for the dust surface density profile")
        self.params["kappa_abs"] = oimParam(name="kappa_abs", value=0,
                                            description="Dust mass absorption coefficient [cm2.g-1]")
        self.params["dist"] = oimParam(name="dist", value=0,
                                       description="Distance of the star [pc]")

        self._t = np.array([0, 1e+99])  # constant value <=> static model
        self._wl = None  # np.array([0.5,1])*1e-6

        # NOTE: Finally call the _eval function that allow the parameters to be processed
        self._eval(**kwargs)

    def _radialProfileFunction(self, r, wl, t):
        dist = self.params["dist"](wl, t)
        rin = self.params["rin"](wl, t)
        rout = self.params["rout"](wl, t)
        q = self.params["q"](wl, t)
        p = self.params["p"](wl, t)
        Tin = self.params["Tin"](wl, t)
        Mdust = self.params["Mdust"](wl, t)*M_sun*1e+3
        kappa_abs = self.params["kappa_abs"](wl, t)

        rin_mas = 1e+3*rin/dist
        rout_mas = 1e+3*rout/dist

        # NOTE: Temperature radial profile
        T_profile = Tin*np.divide(r, rin_mas)**q

        # NOTE: Surface density radial profile
        rin_cm = rin*u.au.to(units.cm)
        rout_cm = rout*u.au.to(units.cm)
        if p == -2:
            sigma_in = Mdust/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f = ((rout_cm/rin_cm)**(p+2)-1)/(p+2)
            sigma_in = Mdust/(2.*np.pi*f*rin_cm**2)
        sigma_profile = sigma_in*np.divide(r, rin_mas)**p

        # NOTE: BB emission radialprofile
        BB_profile = (2*h*c/wl**5)*np.divide(1.,
                                             np.exp((h*c)/(wl*k*T_profile))-1)

        # NOTE: Emissivity factor
        emiss_factor = 1-np.exp(-sigma_profile*kappa_abs)

        # NOTE: Final intensity profile
        return np.nan_to_num(np.logical_and(r > rin_mas, r < rout_mas).astype(int)*BB_profile*emiss_factor, nan=0)

    @property
    def _r(self):
        if False:
            rout = self.params["rout"](self._wl, self._t)
        else:
            rout = self.params["rout"](1e99)
        rmax = 4*rout
        r = np.linspace(0, 1, self._dim)*1e+3*rmax / \
            self.params["dist"](self._wl, self._t)
        return r
