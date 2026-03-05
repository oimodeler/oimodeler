"""Set global options of the oimodeler software."""

from types import SimpleNamespace

import astropy.constants as const

# NOTE: Fourier transform settings
# TODO: Remove this definition (overlapping namespace with astropy; bad form)
constants = SimpleNamespace(
    h=const.h.value,
    c=const.c.value,
    kB=const.k_B.value,
    sigma_sb=const.sigma_sb.value,
    cgs=SimpleNamespace(
        h=const.h.cgs.value, c=const.c.cgs.value, kB=const.k_B.cgs.value
    ),
)
backend = SimpleNamespace(active=None, available=[])
fftw = SimpleNamespace(initialized=False)
ft = SimpleNamespace(backend=backend, binning=None, padding=4, fftw=fftw)

grid = SimpleNamespace(type="linear")
model = SimpleNamespace(grid=grid)

general =  SimpleNamespace(warning=True)

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.

oimOptions = SimpleNamespace(ft=ft, model=model,general=general)
