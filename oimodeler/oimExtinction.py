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

