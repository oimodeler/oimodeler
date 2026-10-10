from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy import interpolate

FITZINDEB = np.genfromtxt(
    Path(__file__).parent / "extlaws" / "FitzIndeb_3.1_VOSA.dat", unpack=True
)
FITZINDEBSPLINE = interpolate.splrep(FITZINDEB[0] / 1e10, FITZINDEB[1], s=1)


def extlaw_FitzIndeb(
    wavelength: float | NDArray[np.floating], A_V: float = 10.0
) -> float | NDArray[np.floating]:
    r"""Compute extinction using the :cite:alp:`Fitzpatrick1999` law with the infrared
    extension of :cite:alp:`Indebetouw2005`.

    The extinction curve is evaluated using a spline representation of the
    combined extinction law. The resulting curve is normalized to the
    visual extinction :math:`A_V`.

    This implementation is valid for wavelengths 0.02-1000 μm.

    Parameters
    ----------
    wavelength : float or NDArray[np.floating]
        Wavelength(s) (m) at which to evaluate the extinction law.
    A_V : float, optional
        Visual extinction in magnitudes. Sets the normalization of the
        extinction curve. Defaults to ``10.0``.

    Returns
    -------
    float or NDArray[np.floating]
        Extinction at the specified wavelength(s) (mag).

    Notes
    -----
    The extinction curve is obtained from `VOSA
    <https://svo2.cab.inta-csic.es/theory/vosa>`_.
    """
    kappa = interpolate.splev(wavelength, FITZINDEBSPLINE, der=0)
    return A_V * (kappa / 211.4)


def extlaw_Cardelli89(
    wavelength: float | NDArray[np.floating],
    A_V: float = 10.0,
    R_V: float = 3.1,
) -> float | NDArray[np.floating]:
    r"""Compute extinction using the :cite:alp:`Cardelli1989` law.

    The extinction is calculated using the piecewise polynomial
    parameterization of the extinction curve as a function of inverse
    wavelength. The curve is normalized by the visual extinction
    :math:`A_V` and parameterized by :math:`R_V`, the ratio of total to
    selective visual extinction.

    This implementation is valid for wavelengths 0.125-3.5 μm.

    Parameters
    ----------
    wavelength : float or NDArray[np.floating]
        Wavelength(s) (m) at which to evaluate the extinction law.
    A_V : float, optional
        Visual extinction in magnitudes. Defaults to ``10.0``.
    R_V : float, optional
        Ratio of total to selective visual extinction,
        :math:`R_V=A_V/E(B-V)`. Controls the shape of the extinction
        curve. Defaults to ``3.1``.

    Returns
    -------
    float or NDArray[np.floating]
        Extinction at the specified wavelength(s) (mag).

        .. math::

            A(\lambda)=A_V\left[a(x)+\frac{b(x)}{R_V}\right],

        where :math:`x = 1/\lambda` when the wavelength is expressed in
        micrometres, and :math:`a(x)` and :math:`b(x)` are the
        wavelength-dependent functions defined by the extinction law.
    """
    x = 1e6 / wavelength
    if np.isscalar(x):
        x = np.array([x])

    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # NOTE: Extend the extinction law smoothly to longer wavelengths; assume the
    # power-law behavior just continues on smoothly
    idx = x <= 1.1
    if np.any(idx):
        a[idx] = 0.574 * x[idx] ** 1.61
        b[idx] = -0.527 * x[idx] ** 1.61

    idx = (1.1 <= x) & (x <= 3.3)
    if np.any(idx):
        y = x[idx] - 1.82
        a[idx] = (
            1.0
            + 0.17699 * y
            - 0.50447 * y**2
            - 0.02427 * y**3
            + 0.72085 * y**4
            + 0.01979 * y**5
            - 0.77530 * y**6
            + 0.32999 * y**7
        )
        b[idx] = (
            1.41338 * y
            + 2.28305 * y**2
            + 1.07233 * y**3
            - 5.38434 * y**4
            - 0.62251 * y**5
            + 5.30260 * y**6
            - 2.09002 * y**7
        )

    idx = (3.3 <= x) & (x <= 8.0)
    if np.any(idx):
        F_a = np.zeros(idx.sum())
        F_b = np.zeros(idx.sum())
        idx2 = (8.0 >= x[idx]) & (x[idx] >= 5.9)
        if np.any(idx2):
            z = x[idx][idx2] - 5.9
            F_a[idx2] = -0.04473 * z**2 - 0.009779 * z**3
            F_b[idx2] = 0.2130 * z**2 + 0.1207 * z**3
        a[idx] = (
            1.752
            - 0.316 * x[idx]
            - 0.104 / ((x[idx] - 4.67) ** 2 + 0.341)
            + F_a
        )
        b[idx] = (
            -3.090
            + 1.825 * x[idx]
            + 1.206 / ((x[idx] - 4.62) ** 2 + 0.263)
            + F_b
        )

    idx = (8.0 <= x) & (x <= 10.0)
    if np.any(idx):
        z = x[idx] - 8.0
        a[idx] = -1.073 - 0.628 * z + 0.137 * z**2 - 0.070 * z**3
        b[idx] = 13.670 + 4.257 * z - 0.420 * z**2 + 0.374 * z**3

    return (a + b / R_V) * A_V
