import astropy.units as u
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam, _standardParameters


class oimRadialPowerLaw(oimComponentImage):
    """A 2D-radial power law distribution

    It accounts for elongation and rotation
    """
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._t, self._wl = [None]*2

        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["pixSize"] = oimParam(**_standardParameters["pixSize"])
        self.params["p"] = oimParam(name="p", value=0, description="Power-law exponent")
        self._eval(**kwargs)

    def _getInternalGrid(self, simple, flatten, wl, t):
        """Creates the internal grid"""
        wl0 = np.sort(np.unique(wl)) if self._wl is None else self._wl
        t0 = np.sort(np.unique(t)) if self._t is None else self._t
        dim = self.params["dim"].value

        pix = self.params["pixSize"](wl, t)*self.params["pixSize"].unit.to(u.mas)
        xy = np.linspace(-0.5, 0.5, dim, endpoint=False)*pix*dim

        if simple:
            return wl, t, xy, xy
        else:
            nwl, nt = map(lambda x: np.array(x).flatten().size, (wl0, t0))

            xx, yy = np.meshgrid(xy,xy)
            x_arr = np.tile(xx[None,None,:,:], (nt,nwl, 1, 1))
            y_arr = np.tile(yy[None,None,:,:], (nt,nwl, 1, 1))
            wl_arr = np.tile(wl[None,:,None,None], (nt, 1, dim, dim))
            t_arr = np.tile(t[:,None,None,None],(1,nwl, dim, dim))

            if flatten:
                return t_arr.flatten(), wl_arr.flatten(),\
                        x_arr.flatten(), y_arr.flatten()
            else:
                return t_arr, wl_arr, x_arr, y_arr

    # def _internalImage(self):
    #     """The internal image that is calculated"""
    #     r=np.linspace(0,(self.params["dim"].value-1)*pix,self.params["dim"].value)
    #     radius = np.linspace(0, self.params["dim"].value, endpoint=False)

    #     computing the pixelSize based on the internal image size and the polar diameter
    #     self.params["pixSize"] = 1.5*dpole/dim*units.mas.to(units.rad)
        # return

    def _imageFunction(self, xx, yy, wl, t):
        """The function describing the component's image"""
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        r, p = np.sqrt(xx**2+yy**2), self.params["p"](wl, t)
        return np.nan_to_num(np.logical_and(r>rin, r<rout).astype(int)*np.divide(r, rin)**p, nan=0)

    def getInternalImage(self, wl, t):
        res = self._internalImage()

        if res is None:
            t_arr, wl_arr, x_arr, y_arr = self._getInternalGrid(simple=False, flatten=False,
                                                                wl=wl, t=t)
            res = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)


        for it in range(res.shape[0]):
            for iwl in range(res.shape[1]):
                res[it, iwl, :, :] = res[it, iwl, :, :]/np.sum(res[it, iwl, :, :])
        return res


    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        """Carries out the FFT"""
        wl = wl if wl is not None else ucoord*0
        t = t if t is not None else ucoord*0

        im0 = self.getInternalImage(wl, t)
        im = self._padImage(im0)
        pix = self.params["pixSize"](wl, t)*self.params["pixSize"].unit.to(u.mas)
        tr = self._ftTranslateFactor(ucoord, vcoord, wl, t)*self.params["f"](wl, t)

        if self._allowExternalRotation == True:
            pa_rad = (self.params["pa"](wl, t))*\
                        self.params["pa"].unit.to(u.rad)
            co = np.cos(pa_rad)
            si = np.sin(pa_rad)
            fxp = ucoord*co-vcoord*si
            fyp = ucoord*si+vcoord*co
            vcoord = fyp
            if self.elliptic == True:
                ucoord = fxp/self.params["elong"](wl, t)
            else:
                ucoord = fxp

        wl0 = np.sort(np.unique(wl)) if self._wl is None else self._wl
        t0 = np.sort(np.unique(t)) if self._t is None else self._t

        if self.FTBackend.check(self.FTBackendData, im, pix, wl0,
                                t0, ucoord, vcoord, wl, t) == False:
            self.FTBackendData = self.FTBackend.prepare(im, pix, wl0, t0,
                                                        ucoord, vcoord, wl, t)

        vc = self.FTBackend.compute(self.FTBackendData, im, pix, wl0,
                                    t0, ucoord, vcoord, wl, t)
        return vc*tr*self.params["f"](wl, t)


class oimAsymmetricRadialPowerLaw(oimRadialPowerLaw):
    """An asymmetrically azimuthally modulated 2D-radial power law distribution

    It accounts for elongation and rotation as well
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.params["a"] = oimParam(name="a", value=0, description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(name="phi", value=0, description="Azimuthal modulation angle")
        self._eval(**kwargs)

    def _azimuthal_modulation(self, x_arr, y_arr, wl, t):
        """Calculates the azimuthal modulation"""
        # TEST: Is it the other way around (y_arr, x_arr)?
        polar_angle = np.arctan2(x_arr, y_arr)
        return (1 + self.params["a"](wl, t)*np.cos(polar_angle-self.params["phi"](wl, t)))

    def _imageFunction(self, xx, yy, wl, t):
        """The function describing the component's image"""
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        r, p = np.sqrt(xx**2+yy**2), self.params["p"](wl, t)
        img = np.nan_to_num(np.logical_and(r>rin, r<rout).astype(int)*np.divide(r, rin)**p, nan=0)
        img *= self._azimuthal_modulation(xx, yy, wl, t)
        return img
