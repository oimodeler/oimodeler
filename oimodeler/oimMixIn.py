from __future__ import annotations

import inspect

import numpy as np
from astropy import units
from numpy.typing import NDArray

from .oimExtinction import extlaw_FitzIndeb
from .oimParam import _standardParameters, oimParam


class EllipticalMixIn:
    """Adds ellipticity to an
    :class:`oimComponent <oimodeler.oimComponent.oimComponent>`."""

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

    def _apply_elliptical_xy(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        wl: NDArray[np.floating] | None = None,
        t: NDArray[np.floating] | None = None,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Applies an elliptical transform to the coordinates axes ``x`` and ``y``.

        The coordinates are first rotated by the position angle. The elongation is
        then applied by multiplying ``x`` axis.
        """
        if not self.elliptic:
            return x, y

        elong = 1 / self.cosi(wl, t) if self.flat else self.elong(wl, t)
        pa_rad = self.pa.qty(wl, t).to(units.rad).value
        co, si = np.cos(pa_rad), np.sin(pa_rad)
        return (x * co - y * si) * elong, x * si + y * co

    def _apply_elliptical_uv(
        self,
        u: NDArray[np.floating],
        v: NDArray[np.floating],
        wl: NDArray[np.floating] | None = None,
        t: NDArray[np.floating] | None = None,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Applies an elliptical transform to the coordinates axes ``u`` and ``v``.

        The coordinates are first rotated by the position angle. The elongation is
        then applied by dividing the ``u`` axis.
        """
        if not self.elliptic:
            return u, v

        elong = 1 / self.cosi(wl, t) if self.flat else self.elong(wl, t)
        pa_rad = self.pa.qty(wl, t).to(units.rad).value
        co, si = np.cos(pa_rad), np.sin(pa_rad)
        return (u * co - v * si) / elong, u * si + v * co


class ExtinctionMixIn:
    """Adds extinction to an
    :class:`oimComponent <oimodeler.oimComponent.oimComponent>`."""

    extincted = False

    def _setup_mixins(self, kwargs) -> None:
        """Adds extinction to component initialisation."""
        super()._setup_mixins(kwargs)

        if "extlaw" in kwargs or kwargs.get("extincted", False):
            self.extincted = True
            self.extargs = []
            self.extlaw = kwargs.get("extlaw", extlaw_FitzIndeb)

            for extarg in inspect.getfullargspec(self.extlaw).args[1:]:
                self.extargs.append(extarg)
                self.params[extarg] = oimParam(
                    **_standardParameters.get(extarg, {"name": extarg})
                )

        # TODO: Remove this after some time after it is standard behaviour
        elif "A_V" in kwargs:
            raise NotImplementedError(
                "Extinction must now be defined by specifying extlaw or extincted, "
                "instead only A_V"
            )

    def _apply_extinction(
        self,
        wl: NDArray[np.floating] | None = None,
    ) -> NDArray[np.floating]:
        """Apply extinction (wavelength dependent if ``wl is not None``)."""
        if not self.extincted:
            return np.array([1.0])

        extinction = self.extlaw(
            wl, *[self.params[extarg].value for extarg in self.extargs]
        )
        return 10 ** (-0.4 * extinction)
