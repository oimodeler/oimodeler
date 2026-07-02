from pathlib import Path
from typing import Union

import numpy as np
from scipy import interpolate

FITZINDEB = np.genfromtxt(
    Path(__file__).parent / "extlaws" / "FitzIndeb_3.1_VOSA.dat", unpack=True
)
FITZINDEBSPLINE = interpolate.splrep(FITZINDEB[0] / 1e10, FITZINDEB[1], s=1)


def extlaw_FitzIndeb(
    wavelength: Union[float, np.ndarray], A_V: float = 10.0
) -> Union[float, np.ndarray]:
    """Extinction law."""
    kappa = interpolate.splev(wavelength, FITZINDEBSPLINE, der=0)
    return A_V * (kappa / 211.4)

def extlaw_Cardelli(wavelength, A_V=10.0, R_V=3.1):

    x = 1e6 / wavelength

    if np.isscalar(x):
        x = np.array([x])

    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # Extend the extinction law smoothly to longer wavelengths; assume the
    # power-law behavior just continues on smoothly
    idx = (x <= 1.1)
    if np.any(idx):
        a[idx] = 0.574 * x[idx]**1.61
        b[idx] = -0.527 * x[idx]**1.61

    idx = (1.1 <= x) & (x <= 3.3)
    if np.any(idx):
        y = x[idx] - 1.82
        a[idx] = 1.0 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
        b[idx] = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7

    idx = (3.3 <= x) & (x <= 8.0)
    if np.any(idx):
        F_a = np.zeros(idx.sum())
        F_b = np.zeros(idx.sum())
        idx2 = (8.0 >= x[idx]) & (x[idx] >= 5.9)
        if np.any(idx2):
            z = x[idx][idx2] - 5.9
            F_a[idx2] = -0.04473 * z**2 - 0.009779 * z**3
            F_b[idx2] = 0.2130 * z**2 + 0.1207 * z**3
        a[idx] = 1.752 - 0.316 * x[idx] - 0.104/((x[idx] - 4.67)**2 + 0.341) + F_a
        b[idx] = -3.090 + 1.825 * x[idx] + 1.206/((x[idx] - 4.62)**2 + 0.263) + F_b

    idx = (8.0 <= x) & (x <= 10.0)
    if np.any(idx):
        z = x[idx] - 8.0
        a[idx] = -1.073 - 0.628 * z + 0.137 * z**2 - 0.070 * z**3
        b[idx] = 13.670 + 4.257 * z - 0.420 * z**2 + 0.374 * z**3

    return (a + b/R_V) * A_V
