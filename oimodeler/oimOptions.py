"""Set global options of the oimodeler software."""
from types import SimpleNamespace

# NOTE: Fourier transform settings
backend = SimpleNamespace(active=None, available=[])
fftw = SimpleNamespace(initialized=False)
ft = SimpleNamespace(
        backend=backend, binning=None,
        padding=1, fftw=fftw)

grid = SimpleNamespace(type="linear")
model = SimpleNamespace(grid=grid)

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.
oimOptions = SimpleNamespace(ft=ft, model=model)
