from __future__ import annotations

import inspect

import numpy as np
from astropy import units
from numpy.typing import NDArray

from .oimExtinction import extlaw_FitzIndeb
from .oimParam import _standardParameters, oimParam


class EllipticalMixIn:
    """Mixin that adds elliptical geometry to an
    :class:`oimComponent <oimodeler.oimComponent.oimComponent>` or its subclasses.

    This mixin enables elliptical coordinate transformations by adding a position
    angle and an elongation parameter to the component. It supports two ways of
    defining the ellipticity:

    - **Elongation-based:** The ``elong`` parameter directly specifies the
      elongation factor.
    - **Inclination-based:** The ``cosi`` parameter specifies the cosine of the
      inclination angle, and the elongation factor is computed as ``1 / cosi``.

    Ellipticity is enabled automatically if ``elliptic=True`` or if any of
    ``pa``, ``elong``, or ``cosi`` are passed during initialization.

    Class Attributes
    ----------------
    elliptic : bool
        Whether elliptical geometry is enabled by default. Defaults to ``False``.
    flat : bool
        Whether to use the inclination-based parameterization. If ``True``,
        the component uses ``cosi`` instead of ``elong``. Defaults to ``False``.

    Notes
    -----
    During component initialization, the mixin adds the following parameters
    to ``self.params`` when ellipticity is enabled:

    - ``pa``: Position angle, represented by an :class:`oimParam <oimodeler.oimParam.oimParam>`.
    - ``elong``: Elongation factor, represented by an :class:`oimParam <oimodeler.oimParam.oimParam>`,
      when using the elongation-based parameterization.
    - ``cosi``: Cosine of the inclination angle, represented by an
      :class:`oimParam <oimodeler.oimParam.oimParam>`, when using the
      inclination-based parameterization.

    The inclination-based parameterization is selected if ``cosi`` is
    provided or if ``flat=True`` is specified.

    The coordinate transformation methods apply the same position-angle
    rotation to both spatial and Fourier coordinates, while accounting for
    the inverse scaling between these coordinate spaces.
    """

    elliptic = False
    flat = False

    def _setup_mixins(self, kwargs) -> None:
        """Configure elliptical geometry during component initialization.

        Parameters
        ----------
        kwargs : dict
            Component initialization arguments. The presence of ``pa``,
            ``elong``, or ``cosi``, or a truthy ``elliptic`` argument,
            enables elliptical geometry. The ``flat`` argument selects
            the inclination-based parameterization.
        """
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
        """Transform spatial coordinates to account for ellipticity.

        The coordinates are first rotated by the position angle ``pa``.
        The rotated x-coordinate is then scaled by the elongation factor,
        while the y-coordinate remains unchanged.

        The elongation factor is either ``elong`` or ``1 / cosi``,
        depending on the selected parameterization.

        Parameters
        ----------
        x, y : NDArray[np.floating]
            Spatial coordinates to transform.
        wl : NDArray[np.floating], optional
            Wavelength values passed to the component parameters. Defaults to ``None``.
        t : NDArray[np.floating], optional
            Time values passed to the component parameters. Defaults to ``None``.

        Returns
        -------
        xtr, ytr : tuple of NDArray[np.floating]
            Transformed spatial coordinates. If ellipticity is disabled,
            ``x`` and ``y`` are returned.
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
        """Transform Fourier coordinates to account for ellipticity.

        The Fourier coordinates are first rotated by the position angle
        ``pa``. The rotated u-coordinate is then divided by the elongation
        factor, while the v-coordinate remains unchanged.

        This inverse scaling is consistent with the corresponding
        transformation in spatial coordinates.

        The elongation factor is either ``elong`` or ``1 / cosi``,
        depending on the selected parameterization.

        Parameters
        ----------
        u, v : NDArray[np.floating]
            Spatial coordinates to transform.
        wl : NDArray[np.floating], optional
            Wavelength values passed to the component parameters. Defaults to ``None``.
        t : NDArray[np.floating], optional
            Time values passed to the component parameters. Defaults to ``None``.

        Returns
        -------
        utr, vtr : tuple of NDArray[np.floating]
            Transformed Fourier coordinates. If ellipticity is disabled,
            `u`` and ``v`` are returned.
        """
        if not self.elliptic:
            return u, v

        elong = 1 / self.cosi(wl, t) if self.flat else self.elong(wl, t)
        pa_rad = self.pa.qty(wl, t).to(units.rad).value
        co, si = np.cos(pa_rad), np.sin(pa_rad)
        return (u * co - v * si) / elong, u * si + v * co


class ExtinctionMixIn:
    r"""Mixin that adds wavelength-dependent extinction to an
    :class:`oimComponent <oimodeler.oimComponent.oimComponent>`
    or any of its subclasses.

    Extinction is applied as a multiplicative attenuation factor to the
    component's emission. The attenuation is calculated from an extinction
    law, which specifies the extinction as a function of wavelength and
    one or more model parameters.

    Extinction can be enabled by providing an extinction law through
    ``extlaw`` or by setting ``extincted=True`` during component
    initialization. If no extinction law is explicitly provided, the
    default law ``extlaw_FitzIndeb`` is used.

    Class Attributes
    ----------------
    extincted : bool
        Whether extinction is enabled by default. Defaults to ``False``.

    Notes
    -----
    When extinction is enabled, the mixin:

    - Stores the selected extinction law in ``self.extlaw``.
    - Identifies the extinction law's parameters and stores their names
      in ``self.extargs``.
    - Adds an :class:`oimParam <oimodeler.oimParam.oimParam>` for each
      extinction-law parameter to ``self.params``.

    The extinction law is expected to accept wavelength as its first
    argument, followed by its model parameters. The law returns extinction
    in magnitudes, which is converted into a multiplicative attenuation
    factor according to

    .. math::

        T(\lambda) = 10^{-0.4 A(\lambda)}

    where :math:`A(\lambda)` is the extinction in magnitudes and
    :math:`T(\lambda)` is the transmitted fraction of the emitted flux.

    Warnings
    --------
    Providing ``A_V`` without specifying ``extlaw`` or ``extincted`` is
    no longer supported and raises a :class:`NotImplementedError`.
    """

    extincted = False

    def _setup_mixins(self, kwargs) -> None:
        """Configure extinction during component initialization.

        Parameters
        ----------
        kwargs : dict
            Component initialization arguments. Extinction is enabled if
            ``extlaw`` is provided or ``extincted=True`` is specified.
            The ``extlaw`` argument selects the extinction law. If omitted,
            ``extlaw_FitzIndeb`` is used when extinction is enabled.

        Raises
        ------
        NotImplementedError
            If ``A_V`` is provided without ``extlaw`` or ``extincted``.
            Specifying ``A_V`` alone is deprecated and no longer supported.
        """
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
        """Compute the multiplicative attenuation due to extinction.

        The selected extinction law is evaluated at the supplied
        wavelengths using the current values of its model parameters.
        The resulting extinction, expressed in magnitudes, is converted
        into a multiplicative transmission factor.

        Parameters
        ----------
        wl : NDArray[np.floating], optional
            Wavelengths at which to evaluate the extinction law.
            Defaults to ``None``.

        Returns
        -------
        NDArray[np.floating]
            Multiplicative attenuation factor computed as
            ``10**(-0.4 * extinction)``. If extinction is disabled,
            returns an array containing a single value of ``1.0``.

        Notes
        -----
        An attenuation factor of ``1.0`` corresponds to no attenuation.
        """
        if not self.extincted:
            return np.array([1.0])

        extinction = self.extlaw(
            wl, *[self.params[extarg].value for extarg in self.extargs]
        )
        return 10 ** (-0.4 * extinction)
