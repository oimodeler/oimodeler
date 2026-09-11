"""Set global options of the oimodeler software."""

from types import SimpleNamespace

import astropy.constants as const
import astropy.units as u

# NOTE: Physical constants
CGS = SimpleNamespace(
    H=const.h.cgs.value, C=const.c.cgs.value, K_B=const.k_B.cgs.value
)
SI = SimpleNamespace(
    H=const.h.value,
    C=const.c.value,
    K_B=const.k_B.value,
    SIGMA_SB=const.sigma_sb.value,
)

# NOTE: Conversion constants
ARCSEC2RAD: float = u.arcsec.to(u.rad)
MAS2RAD: float = u.mas.to(u.rad)
M2AU: float = u.m.to(u.au)
RAD2MAS: float = 1 / MAS2RAD

# NOTE: Fourier transform settings
backend = SimpleNamespace(active=None, available=[])
fftw = SimpleNamespace(initialized=False)
ft = SimpleNamespace(backend=backend, binning=None, padding=4, fftw=fftw)

grid = SimpleNamespace(type="logarithmic")
model = SimpleNamespace(grid=grid)
general = SimpleNamespace(warning=True)

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.

oimOptions = SimpleNamespace(ft=ft, model=model, general=general)
