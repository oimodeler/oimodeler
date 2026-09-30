# -*- coding: utf-8 -*-
"""Components defined in Fourier or image planes"""

from __future__ import annotations

import copy
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy import units
from astropy.io import fits
from numpy.typing import ArrayLike, NDArray
from scipy import interpolate
from scipy.special import j0, jv

from . import __dict__ as oimDict
from .oimExtinction import ExtinctionMixIn
from .oimOptions import MAS2RAD, RAD2MAS, oimOptions
from .oimParam import (
    _standardParameters,
    oimInterp,
    oimParam,
    oimParamInterpolator,
    oimParamLinker,
    oimParamNorm,
)
from .oimUtils import (
    _deserialize_function,
    _pickle,
    _serialize_function,
    _unpickle,
    attach_methods,
    getWlFromFitsImageCube,
    pad_image,
)

EXEMPTED_KEYS: list[str] = ["asymmetric", "compute_sigma0", "modulation"]


# TODO: Move somewhere else
def getFourierComponents():
    """A function to get the list of all available components deriving from the
    oimComponentFourier class

    Returns
    -------
    res : list
        list of all available components deriving from the oimComponentFourier class.
    """
    fnames = dir()
    res = []
    for f in fnames:
        try:
            if issubclass(oimDict[f], oimComponentFourier):
                res.append(f)
        except:
            pass
    return res


class EllipticalMixIn:
    """Adds ellipticity to an `oimComponent <oimodeler.oimComponent.oimComponent>`."""

    elliptic = False
    flat = False

    def _setup_mixins(self, kwargs) -> None:
        """Adds ellipticity to component initialisation."""
        super()._setup_mixins(kwargs)

        elliptic = kwargs.get("elliptic", self.elliptic)
        if any(x in kwargs for x in ["cosi", "elong", "pa"]) or elliptic:
            self.elliptic = True
            self.params["pa"] = oimParam(base="pa")
            if "cosi" in kwargs or kwargs.get("flat", self.flat):
                self.flat = True
                self.params["cosi"] = oimParam(base="cosi")
            else:
                self.params["elong"] = oimParam(base="elong")

    def _apply_elliptical(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        plane: str,
        wl: NDArray[np.floating] | None = None,
        t: NDArray[np.floating] | None = None,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Applies an elliptical transform to the coordinates axes ``x`` and ``y``.

        The coordinates are first rotated by the position angle.
        The elongation is then applied by either by either multiplying
        (``plane="image"``) or dividing (``plane="fourier"``) the ``x`` axis.
        """
        xtr, ytr = x, y
        if not self.elliptic:
            return x, y

        elong = 1 / self.cosi(wl, t) if self.flat else self.elong(wl, t)
        pa_rad = self.pa.qty(wl, t).to(units.rad).value
        co, si = np.cos(pa_rad), np.sin(pa_rad)
        xtr, ytr = x * co - y * si, x * si + y * co

        if plane == "image":
            xtr *= elong
        elif plane == "fourier":
            xtr /= elong
        else:
            raise ValueError(
                'Only "image" and "fourier" are valid arguments for plane.'
            )

        return xtr, ytr


# TODO: Implement these attaches differently or not at all?
@attach_methods({"pickle": _pickle, "unpickle": classmethod(_unpickle)})
class oimComponent:
    """The oimComponent class is the abstract parent class for all types of
    components that can be added to an oimModel.

    It has a similar interface as the oimModel and allows to compute images
    (or image cubes for wavelength-dependent or time-dependent models)
    and complex-coherent fluxes for a vector of ``(t,wl,u,v)`` coordinates.

    Parameters
    ----------
    x : float or oimInterp
        x pos of the component (mas). Defaults to ``0``.
    y : float or oimInterp
        y pos of the component (mas). Defaults to ``0``.
    f : float or oimInterp
        Flux (ratio) of the component. Defaults to ``1``.

    Attributes
    ----------
    name : str
        Name of the component.
    shortname : str
        Short name for the component.
    description : str
        Description of the component.
    params : dict of str to oimParam
        Dictionary of the component parameters.
    x : oimParam
        x pos of the component (mas).
    y : oimParam
        y pos of the component (mas).
    f : oimParam
        Flux (ratio) of the component.
    """

    _firstInit = True
    name = "Generic component"
    shortname = "GenComp"
    description = "Class from which all components are derived"

    def __init__(self, **kwargs):
        """Create and initialize an instance of the oimComponent class."""
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static

        self.params = {}
        self.params["x"] = oimParam(**_standardParameters["x"])
        self.params["y"] = oimParam(**_standardParameters["y"])
        self.params["f"] = oimParam(**_standardParameters["f"])
        self._setup_mixins(kwargs)
        self._eval(**kwargs, checkParam=False)

    def _setup_mixins(self, kwargs) -> None:
        """Set up mixins. That is, additional capabilities for the component
        (e.g. ellipticity, extinction)."""

    def _paramstr(self):
        txt = []
        for paramname, param in self.params.items():
            if isinstance(param, oimParam):
                if isinstance(param, oimParamInterpolator):
                    # TODO: Have a string for each oimParamInterpolator
                    txt.append(f"{param.name}={param.__class__.__name__}")
                else:
                    if np.abs(param.value) < 1e-2 and param.value != 0:
                        value_string = f"{param.value:.2e}"
                    else:
                        value_string = f"{param.value:.2f}"

                    txt.append(f"{param.name}={value_string}")

            elif isinstance(param, oimParamNorm) or isinstance(
                param, oimParamLinker
            ):
                txt.append(f"{paramname}={param.__class__.__name__}")

        return " ".join(txt)

    def __str__(self):
        return self.name + ": " + self._paramstr()

    def __repr__(self):
        return (
            f"{self.__class__.__name__} at "
            f"{str(hex(id(self)))}: {self._paramstr()}"
        )

    @property
    def _wl(self) -> np.ndarray:
        """Gets the wavelengths."""
        return self.__wl

    # NOTE: .__wl or .__t is not the cleanest approach for serialisation (MBS).
    # Python adds, for instance, _oimComponent__wl automatically, making it harder
    # to directly serialise
    @_wl.setter
    def _wl(self, value: Any) -> np.ndarray | None:
        """Sets the wavelengths."""
        if value is None:
            self.__wl = None
            return None

        if isinstance(value, (np.ndarray, tuple, list)):
            value = value
        elif isinstance(value, u.Quantity):
            if not isinstance(value.value, (np.ndarray, tuple, list)):
                value = [value]
        else:
            value = [value]
        self.__wl = np.array(value)

    @property
    def _t(self) -> np.ndarray:
        """Gets the times."""
        return self.__t

    @_t.setter
    def _t(self, value: Any) -> np.ndarray | None:
        """Sets the times."""
        if value is None:
            self.__t = None
            return None

        if isinstance(value, (np.ndarray, tuple, list)):
            value = value
        elif isinstance(value, u.Quantity):
            if not isinstance(value.value, (np.ndarray, tuple, list)):
                value = [value]
        else:
            value = [value]
        self.__t = np.array(value)

    def _eval(self, checkParam=True, **kwargs) -> None:
        for key, value in kwargs.items():
            if key in self.params.keys():
                if isinstance(value, oimInterp):
                    if not isinstance(self.params[key], oimParamInterpolator):
                        self.params[key] = value.type(
                            self.params[key], **value.kwargs
                        )
                else:
                    self.params[key].value = value
                    if isinstance(value, u.Quantity):
                        self.params[key].value = value.value
                        self.params[key].unit = value.unit

            elif checkParam and key not in EXEMPTED_KEYS:
                warnings.warn(f"{key} not a parameter of {self.name}: ignored")

        for key in self.params.keys():
            prop = property(
                lambda self, k=key: self.params[k],
                lambda self, v, k=key: self.params.__setitem__(k, v),
            )
            setattr(type(self), key, prop)

    def getComplexCoherentFlux(
        self,
        u: ArrayLike,
        v: ArrayLike,
        wl: ArrayLike | None = None,
        t: ArrayLike | None = None,
    ) -> NDArray[np.floating]:
        """Compute and return the complex coherent flux for an array of u,v
        (and optionally wavelength and time ) coordinates

        Parameters
        ----------
        u : array_like
            Spatial coordinate u (cycles/rad).
        v : array_like
            Spatial coordinate v (cycles/rad).
        wl : array_like, optional
            Wavelength (m). Defaults to ``None``.
        t :  array_like, optional
            Time (mjd). Defaults to ``None``.

        Returns
        -------
        NDArray[np.floating]
            The complex coherent flux.
        """
        return np.array(u) * 0

    def getImage(
        self,
        dim: int,
        pixSize: float,
        wl: ArrayLike | None = None,
        t: ArrayLike | None = None,
    ) -> NDArray[np.floating]:
        """Compute an image or image cube (if wavelength and time are given).

        Parameters
        ----------
        dim : integer
            Image dimension (pixels).
        pixSize : float
            Pixel angular size (rad).
        wl : integer, list or numpy array, optional
             Wavelength (m). Defaults to ``None``.
        t :  integer, list or numpy array, optional
            Time (mjd). Defaults to ``None``.

        Returns
        -------
        image : NDArray[np.floating]
        """
        return np.zeros((dim, dim))

    def _ftTranslateFactor(self, ucoord, vcoord, wl, t) -> np.ndarray:
        x = self.params["x"](wl, t) * self.params["x"].unit.to(units.rad)
        y = self.params["y"](wl, t) * self.params["y"].unit.to(units.rad)
        return np.exp(-2 * 1j * np.pi * (ucoord * x + vcoord * y))

    def _directTranslate(self, x, y, wl, t):
        return x - self.params["x"](wl, t), y - self.params["y"](wl, t)

    def getNonRegularImage(self, xx, yy, wl=None, t=None):
        """Compute and return a non-regular image function at the xx, yy and
        optional wl and t coordinates)"""
        return 0 * xx

    def serialize(self, skip_copy: bool = False) -> dict[str, Any]:
        """Serializes the oimComponent.

        Parameters
        ----------
        skip_copy : bool, optional
            If ``True`` skips the top-level deepcopy of oimComponent.
            Sub-level deepcopies (e.g. oimParam) are skipped by default.
            Defaults to ``False``.
        """
        ser = {"params": {}, "other": {}}
        params = self.params
        if not skip_copy:
            params = copy.deepcopy(params)

        for name, param in params.items():
            ser["params"][name] = param.serialize(skip_copy=True)

        # TODO: Does this also need a deepcopy?
        ser["other"] = {
            k: v
            for k, v in self.__class__.__dict__.items()
            if not (k.startswith("_") or isinstance(v, (property, Callable)))
        }

        for key, value in vars(self).items():
            # TODO: This might not work for sub-sub components.
            # HACK: Solution for Python renaming private (e.g. __wl) components
            # with the prefix of the class they are private in.
            key = key.replace("_oimComponent_", "")
            if key == "params":
                continue

            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif isinstance(value, Callable):
                value = _serialize_function(value)

            ser["other"][key] = value

        return ser

    @classmethod
    def deserialize(cls, ser: dict[str, Any]) -> "oimComponent":
        """Deserializes into an oimComponent."""
        cls = copy.deepcopy(cls)
        # HACK: This makes sure that that things like "elliptic", "flat",
        # "extincted", or any future additions are read in first as they
        # set/enable/remove certain parameters
        for key, value in ser["other"].items():
            if isinstance(value, list):
                value = np.array(value)

            setattr(cls, key, value)

        comp = cls()
        # TODO: Merge this and the above loop if possible?
        for key, value in ser["other"].items():
            if isinstance(value, str) and value.startswith("fn::"):
                value = _deserialize_function(value)
            if isinstance(value, list):
                value = np.array(value)

            setattr(comp, key, value)

        for key, value in ser["params"].items():
            comp.params[key] = oimParam.deserialize(value)
            prop = property(
                lambda self, k=key: self.params[k],
                lambda self, v, k=key: self.params.__setitem__(k, v),
            )
            setattr(type(comp), key, prop)

        return comp

    def _fov(self, wl=None, t=None):
        return 0

    def getFOV(self, wl=None, t=None):

        fov0 = np.array(self._fov(wl, t)).max()

        fovX = np.array([-fov0 / 2, fov0 / 2])
        fovY = np.array([-fov0 / 2, fov0 / 2])

        fovX += self.params["x"].value
        fovY += self.params["y"].value

        return [*fovX, *fovY]

    def _solidAngle(self, wl=None, t=None):
        return None

    def getSolidAngle(self, wl=None, t=None):
        sa = self._solidAngle(wl, t)
        try:
            return sa
        except:
            raise NotImplementedError(
                f"getSolidAngle not implement for Class {self.__class__.__name__}"
            )


class oimComponentFourier(ExtinctionMixIn, EllipticalMixIn, oimComponent):
    """Class for all components analytically defined in the Fourier plane.

    Notes
    -----
    Inherits from the `oimComponent`. Has ellipticity and extinction support
    via the `EllipticalMixIn` and `ExtinctionMixIn` classes.

    Implements translation in direct and Fourier space, `getImage` from the
    Fourier definition of the object, ellipticity (i.e. flattening).
    Children classes should only implement the `_visFunction` and `_imageFunction`
    """

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        fxp, fyp = self._apply_elliptical(ucoord, vcoord, "fourier", wl, t)
        return (
            self._visFunction(fxp, fyp, np.hypot(fxp, fyp), wl, t)
            * self._ftTranslateFactor(ucoord, vcoord, wl, t)
            * self.params["f"](wl, t)
            * self._apply_extinction(wl)
        )

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        return ucoord * 0

    def getImage(self, dim, pixSize, wl=None, t=None):
        t = np.array(t).flatten()
        nt = t.size
        wl = np.array(wl).flatten()
        nwl = wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)
        x_arr, y_arr = self._apply_elliptical(
            x_arr, y_arr, "image", wl_arr, t_arr
        )
        image = (
            self._imageFunction(
                x_arr.reshape(dims),
                y_arr.reshape(dims),
                wl_arr.reshape(dims),
                t_arr.reshape(dims),
            )
            * self._apply_extinction(wl)[np.newaxis, :, np.newaxis, np.newaxis]
        )

        tot = np.sum(image, axis=(2, 3))

        # TODO: No loop for normalization
        for it, ti in enumerate(t):
            for iwl, wli in enumerate(wl):
                if tot[it, iwl] != 0:
                    image[it, iwl, :, :] = (
                        image[it, iwl, :, :]
                        / tot[it, iwl]
                        * self.params["f"](wli, ti)
                    )
        return image

    def getNonRegularImage(self, xx, yy, wl=None, t=None):
        xx, yy = self._directTranslate(xx, yy, wl, t)
        xx, yy = self._apply_elliptical(xx, yy, "image", wl, t)
        return (
            self._imageFunction(xx, yy, wl, t)
            * self._apply_extinction(wl)[np.newaxis, :, np.newaxis, np.newaxis]
        )

    def _imageFunction(self, xx, yy, wl, t):
        raise ValueError(
            f"image function not implemented for {self.shortname}.\n"
            "Use the fromFT=True option to get a model image"
            " from the inverse Fourier Transform"
        )


class oimComponentImage(ExtinctionMixIn, EllipticalMixIn, oimComponent):
    """Base class for components define in 2D : x,y (regular grid) in the image plan"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.interpFillValue = 0

        # NOTE: In rad
        self._pixSize = 0
        self._allowExternalRotation = True
        self.normalizeImage = True
        self.params["pa"] = oimParam(**_standardParameters["pa"])
        self.params["dim"] = oimParam(**_standardParameters["dim"])

        if "FTBackend" in kwargs:
            self.FTBackend = kwargs["FTBackend"]()
        else:
            self.FTBackend = oimOptions.ft.backend.active()

        self.FTBackendData = None
        self._eval(**kwargs, checkParam=False)

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if wl is None:
            wl = ucoord * 0
        if t is None:
            t = ucoord * 0

        im = self.getInternalImage(wl, t)
        im = pad_image(im)

        if self._pixSize != 0:
            pix = self._pixSize
        else:
            pix = self.getPixelSize()

        tr = self._ftTranslateFactor(ucoord, vcoord, wl, t)

        fxp, fyp = ucoord, vcoord
        if self._allowExternalRotation:
            pa_rad = (self.params["pa"](wl, t)) * self.params["pa"].unit.to(
                units.rad
            )
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = ucoord * co - vcoord * si
            fyp = ucoord * si + vcoord * co

            if self.elliptic:
                if self.flat:
                    fxp *= self.params["cosi"](wl, t)
                else:
                    fxp /= self.params["elong"](wl, t)

        if self._wl is None:
            wl0 = np.unique(wl)
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.unique(t)
        else:
            t0 = self._t

        if not (
            self.FTBackend.check(
                self.FTBackendData, im, pix, wl0, t0, fxp, fyp, wl, t
            )
        ):

            self.FTBackendData = self.FTBackend.prepare(
                im, pix, wl0, t0, fxp, fyp, wl, t
            )

        vc = self.FTBackend.compute(
            self.FTBackendData, im, pix, wl0, t0, fxp, fyp, wl, t
        )

        return vc * tr * self.params["f"](wl, t) * self._apply_extinction(wl)

    def getImage(self, dim, pixSize, wl=None, t=None):
        if wl is None:
            wl = 0
        if t is None:
            t = 0

        t = np.array(t).flatten()
        nt = t.size
        wl = np.array(wl).flatten()
        nwl = wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)
        # TODO: Figure out how to handle this with the overarching framework
        if self._allowExternalRotation:
            pa_rad = (self.params["pa"](wl_arr, t_arr)) * self.params[
                "pa"
            ].unit.to(units.rad)

            xp = x_arr * np.cos(pa_rad) - y_arr * np.sin(pa_rad)
            yp = x_arr * np.sin(pa_rad) + y_arr * np.cos(pa_rad)

            x_arr, y_arr = xp, yp
            if self.elliptic:
                if self.flat:
                    x_arr /= self.params["cosi"](wl_arr, t_arr)
                else:
                    x_arr *= self.params["elong"](wl_arr, t_arr)

        im0 = self._internalImage()
        if im0 is None:
            im = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)
        else:
            im0 = np.swapaxes(im0, -2, -1)
            grid = self._getInternalGrid()
            coord = np.transpose(np.array([t_arr, wl_arr, x_arr, y_arr]))

            im = interpolate.interpn(
                grid,
                im0,
                coord,
                bounds_error=False,
                fill_value=self.interpFillValue,
            )
            f0 = np.sum(im0)
            f = np.sum(im)
            im = im / f * f0

        im = (
            im.reshape(dims)
            * self._apply_extinction(wl)[np.newaxis, :, np.newaxis, np.newaxis]
        )

        # TODO: No loop for normalization
        if self.normalizeImage:
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = (
                            im[it, iwl, :, :]
                            / tot[it, iwl]
                            * self.params["f"](wli, ti)
                        )
        return im

    def getInternalImage(self, wl=None, t=None):
        res = self._internalImage()

        if res is None:
            t_arr, wl_arr, x_arr, y_arr = self._getInternalGrid(
                simple=False, wl=wl, t=t
            )
            res = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)

        # TODO: No loop for normalization
        if self.normalizeImage:
            for it in range(res.shape[0]):
                for iwl in range(res.shape[1]):
                    res[it, iwl, :, :] = res[it, iwl, :, :] / np.sum(
                        res[it, iwl, :, :]
                    )

        return res

    def _internalImage(self):
        return

    def _imageFunction(self, xx, yy, wl, t):
        image = xx * 0 + 1
        return image

    def _getInternalGrid(self, simple=True, flatten=False, wl=None, t=None):
        if self._wl is None:
            wl0 = np.unique(wl)
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.unique(t)
        else:
            t0 = self._t

        dim = self.params["dim"](wl, t)
        if self._pixSize != 0:
            pix = self._pixSize * RAD2MAS
        else:
            pix = self.getPixelSize() * RAD2MAS

        v = np.linspace(-0.5, 0.5, dim)
        xy = v * pix * dim

        if simple:
            return t0, wl0, xy, xy

        else:
            t = np.array(t0).flatten()
            nt = t.size
            wl = np.array(wl0).flatten()
            nwl = wl.size
            xx, yy = np.meshgrid(xy, xy)
            x_arr = np.tile(xx[None, None, :, :], (nt, nwl, 1, 1))
            y_arr = np.tile(yy[None, None, :, :], (nt, nwl, 1, 1))
            wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
            t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

            if flatten:
                return (
                    t_arr.flatten(),
                    wl_arr.flatten(),
                    x_arr.flatten(),
                    y_arr.flatten(),
                )
            else:
                return t_arr, wl_arr, x_arr, y_arr

    def getPixelSize(self, mas=False):
        raise ValueError(
            "getPixelSize Method not implemented"
            " while self._pixSize = "
            f"{self._pixSize}"
        )

    def _fov(self, wl=None, t=None):
        return self.getPixelSize() * RAD2MAS * self.params["dim"].value


class oimComponentRadialProfile(
    ExtinctionMixIn, EllipticalMixIn, oimComponent
):
    """Base class for components defined by a radial profile."""

    asymmetric = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static
        self.normalizeImage = True
        self.params["dim"] = oimParam(base="dim")
        self._r, self._dr = None, None
        self._eval(**kwargs, checkParam=False)

    @property
    def r(self) -> None | NDArray[np.float64]:
        """Gets the radial profile (mas)."""
        return self._r

    @property
    def dr(self) -> None | NDArray[np.float64]:
        """Gets the integration weights (mas)."""
        if self._dr is None:
            self._dr = np.gradient(self.r)

        return self._dr

    def _getInternalGrid(self, simple=True, flatten=False, wl=None, t=None):
        wl0 = np.unique(wl) if self._wl is None else self._wl
        t0 = np.unique(t) if self._t is None else self._t

        r = self.r
        if r is None:
            pix = self._pixSize * RAD2MAS
            r = np.linspace(0, self.dim.value - 1, self.dim.value) * pix

        if simple:
            return r, wl, t
        else:
            nt = np.array(t0).flatten().size
            nwl = np.array(wl0).flatten().size
            nr = r.flatten().size

            r_arr = np.tile(r[None, None, :], (nt, nwl, 1))
            wl_arr = np.tile(wl[None, :, None], (nt, 1, nr))
            t_arr = np.tile(t[:, None, None], (1, nwl, nr))

            if flatten:
                return t_arr.flatten(), wl_arr.flatten(), r_arr.flatten()
            else:
                return t_arr, wl_arr, r_arr

    def _internalRadialProfile(self):
        return None

    def _radialProfileFunction(self, r=None, wl=None, t=None):
        return 0

    def getInternalRadialProfile(self, wl, t):
        res = self._internalRadialProfile()
        if res is None:
            t_arr, wl_arr, r_arr = self._getInternalGrid(
                simple=False, wl=wl, t=t
            )
            res = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        return res

    def getImage(self, dim, pixSize, wl=None, t=None):
        wl, t = 0 if wl is None else wl, 0 if t is None else t
        t, wl = np.array(t).flatten(), np.array(wl).flatten()
        nt, nwl = t.size, wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)
        x_arr, y_arr = self._apply_elliptical(
            x_arr, y_arr, "image", wl_arr, t_arr
        )
        r_arr = np.hypot(x_arr, y_arr)
        im = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        im = np.nan_to_num(
            im.reshape(dims)
            * self._apply_extinction(wl)[np.newaxis, :, np.newaxis, np.newaxis]
        )

        if self.normalizeImage:
            # TODO: No loop for normalization
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = (
                            im[it, iwl, :, :] / tot[it, iwl] * self.f(wli, ti)
                        )

        return im

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        wl = ucoord * 0 if wl is None else wl
        t = ucoord * 0 if t is None else t

        # TODO: Performance: Move the `np.unique` lines into ``oimData``
        wl0, idx_wl = np.unique(wl, return_inverse=True)
        t0, idx_t = np.unique(t, return_inverse=True)
        (ucoord0, vcoord0), idx_uvcoord = np.unique(
            np.vstack((ucoord, vcoord)), return_inverse=True, axis=1
        )
        ucoord0, vcoord0 = self._apply_elliptical(
            ucoord0, vcoord0, "fourier", wl0, t0
        )
        Ir0 = self.getInternalRadialProfile(wl0, t0)
        r = self.r * MAS2RAD
        kr = (
            2.0
            * np.pi
            * r[:, np.newaxis]
            * np.hypot(ucoord0, vcoord0)[np.newaxis, :]
        )
        kernel = j0(kr)

        # FIXME: Not yet tested for correct output values.
        if self.asymmetric:
            psi = np.arctan2(ucoord0, vcoord0)
            for i in range(1, self.modulation + 1):
                skwi = getattr(self, f"skw{i}")(wl, t)
                skwPai = getattr(self, f"skwPa{i}").qty(wl, t).to(u.rad).value
                kernel += (
                    (-1j) ** i * skwi * np.cos(i * (psi - skwPai)) * jv(i, kr)
                )

        kernel *= (2 * np.pi * r * self.dr * MAS2RAD)[:, np.newaxis]

        # TODO: Grid is overcomputed: (nwl * nuv[m]) < (nwl * nuv[cycle/rad])
        vc0 = Ir0 @ kernel * 1e23 + 0j
        vc0 *= (
            self._ftTranslateFactor(
                ucoord0[np.newaxis, np.newaxis],
                vcoord0[np.newaxis, np.newaxis],
                wl0[np.newaxis, :, np.newaxis],
                t0[:, np.newaxis, np.newaxis],
            )
            * self.f(wl0, t0)
            * self._apply_extinction(wl0)[np.newaxis, :, np.newaxis]
        )

        # FIXME: Test if correct for ``(Ir0.shape[0] = nt0) != 1``
        return vc0[idx_t if Ir0.shape[0] != 1 else 0, idx_wl, idx_uvcoord]


class oimComponentFitsImage(oimComponentImage):
    """Component load load images or chromatic-cubes from fits files"""

    name = "Fits Image Component"
    shortname = "Fits_Comp"

    def __init__(self, fitsImage=None, useinternalPA=False, **kwargs):
        super().__init__(**kwargs)
        if fitsImage:
            self.loadImage(fitsImage, useinternalPA=useinternalPA)

        self.params["pa"] = oimParam(**_standardParameters["pa"])
        self.params["scale"] = oimParam(**_standardParameters["scale"])
        self._eval(**kwargs)

    def loadImage(self, fitsImage, useinternalPA=False):
        if isinstance(fitsImage, (str, Path)):
            try:
                im = fits.open(fitsImage)[0]
            except:
                raise TypeError("Not a valid fits file")

        elif isinstance(fitsImage, fits.hdu.hdulist.HDUList):
            im = fitsImage[0]
        elif isinstance(fitsImage, fits.hdu.image.PrimaryHDU):
            im = fitsImage

        self._header = im.header
        dims = self._header["NAXIS"]
        if dims < 2:
            raise TypeError(
                "oimComponentFitsImage require 2D images or "
                "3D chromatic-image-cubes"
            )

        dimx, dimy = self._header["NAXIS1"], self._header["NAXIS2"]
        if dimx != dimy:
            raise TypeError("Current version only works with square images")
        self._dim = dimx

        pixX, pixY = self._header["CDELT1"], self._header["CDELT2"]
        if pixX != pixY:
            raise TypeError(
                "Current version only works with the same pixels"
                " scale in x and y dimension"
            )
        if "CUNIT1" in self._header:
            try:
                unit0 = units.Unit(self._header["CUNIT1"])
            except:
                unit0 = units.rad

        else:
            unit0 = units.rad
        self._pixSize0 = pixX * unit0.to(units.rad)

        if "CROTA1" in self._header:
            pa0 = self._header["CROTA1"]
        elif "CROTA2" in self._header:
            pa0 = self._header["CROTA2"]

        # TODO: Check units
        if useinternalPA:
            self.params["pa"].value = pa0

        # NOTE: Adding the time dimension (nt,nwl,ny,nx)
        if dims == 3:
            self._wl = getWlFromFitsImageCube(self._header, units.m)
            self._image = im.data[None, :, :, :]
        else:
            self._image = im.data[None, None, :, :]
            self._wl = np.array([0])

    def _internalImage(self):
        self.params["dim"].value = self._dim
        self._pixSize = self._pixSize0 * self.params["scale"].value
        return self._image

    def getPixelSize(self, mas=False):
        self._pixSize = self._pixSize0 * self.params["scale"].value
        return self._pixSize * (RAD2MAS * mas + (not mas))
