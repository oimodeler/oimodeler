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
    """Optical Extinction law [1]_, with infrared extension [2]_.

    Appropriate for 0.02-1000 µm.

    Parameters
    ----------
    wavelength : float or NDArray[np.floating]
        Wavelength(s) to compute the extinction for.
    A_V : float, optional
        Defaults to ``10.0``.

    Returns
    -------
    float or NDArray[np.floating]

    Notes
    -----
    Obtained from `VOSA <https://svo2.cab.inta-csic.es/theory/vosa>`_.

    References
    ----------
    .. [1] Fitzpatrick, "Correcting for the Effects of Interstellar Extinction",
       PASP, Volume 111, id. 755, 63-75 pp., (1999).

    .. [2] Indebetouw et al., "The Wavelength Dependence of Interstellar Extinction
       from 1.25 to 8.0 Microns", ApJ, Volume 619, id. 2, 931-938 pp. (2005).
    """
    kappa = interpolate.splev(wavelength, FITZINDEBSPLINE, der=0)
    return A_V * (kappa / 211.4)


def extlaw_Cardelli89(
    wavelength: float | NDArray[np.floating],
    A_V: float = 10.0,
    R_V: float = 3.1,
) -> float | NDArray[np.floating]:
    """Extinction law [1]_.

    Appropriate for 0.125-3.5 µm.

    Parameters
    ----------
    wavelength : float or NDArray[np.floating]
        Wavelength(s) to compute the extinction for.
    A_V : float, optional
        Defaults to ``10.0``.
    R_V : float, optional
        Defaults to ``3.1``.

    Returns
    -------
    float or NDArray[np.floating]

    References
    ----------
    .. [1] Cardelli et al, "The Relationship between Infrared,
    Optical, and Ultraviolet Extinction", ApJ, Volume 345, 245-256 pp., (1989).
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
