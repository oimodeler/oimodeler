"""Set global options of the oimodeler software."""
from types import SimpleNamespace

import astropy.units as u


# NOTE: Here is a list of standard parameters to be used when defining new components
standard_parameters = SimpleNamespace(
    x={"name": "x", "value": 0, "description": "x position", "unit": u.mas, "free": False},
    y={"name": "y", "value": 0, "description": "y position", "unit": u.mas, "free": False},
    f={"name": "f", "value": 1, "description": "flux", "unit": u.one, "mini": 0, "maxi": 1},
    fwhm={"name": "fwhm", "value": 0, "description": "FWHM", "unit": u.mas, "mini": 0},
    d={"name": "d", "value": 0, "description": "Diameter", "unit": u.mas, "mini": 0},
    din={"name": "din", "value": 0, "description": "Inner Diameter", "unit": u.mas, "mini": 0},
    dout={"name": "dout", "value": 0, "description": "Outer Diameter", "unit": u.mas, "mini": 0},
    elong={"name": "elong", "value": 1, "description": "Elongation Ratio", "unit": u.one, "mini": 1},
    pa={"name": "pa", "value": 0, "description": "Major-axis Position angle", "unit": u.deg, "mini": -180, "maxi": 180},
    skw={"name": "skw", "value": 0, "description": "Skewedness", "unit": u.one, "mini": 0, "maxi": 1},
    skwPa={"name": "skwPa", "value": 0, "description": "Skewedness Position angle", "unit": u.deg, "mini": -180, "maxi": 180},
    pixSize={"name": "pixSize", "value": 0, "description": "Pixel Size", "unit": u.mas, "free": False, "mini": 0},
    dim={"name": "dim", "value": 128, "description": "Dimension in pixels", "unit": u.one, "free": False, "mini": 1},
    wl={"name": "wl", "value": 0, "description": "Wavelength", "unit": u.m, "mini": 0},
    mjd={"name": "mjd", "value": 0, "description": "MJD", "unit": u.day},
    scale={"name": "scale", "value": 1, "description": "Scaling Factor", "unit": u.one},
    index={"name": "index", "value": 1, "description": "Index", "unit": u.one},
    fov={"name": "fov", "value": 0, "description": "The interferometric field of view", "unit": u.mas, "free": False, "mini": 0},
    amp={"name": "amplitude", "value": 1, "description": "Amplitude", "unit": u.one},
    p={"name": "p", "value": 0, "description": "Power-law Exponent", "unit": u.one},
)

# NOTE: Sets the available interpolators for oimodeler. If strings are provided,
# the `oimInterp` looks through `oimParam` in order to find the class.
# To overwrite provide class variables
interpolators = SimpleNamespace(
        wl="oimParamInterpolatorWl",
        time="oimParamInterpolatorTime",
        GaussWl="oimParamGaussianWl",
        GaussTime="oimParamGaussianTime",
        mGaussWl="oimParamMultipleGaussianWl",
        mGaussTime="oimParamMultipleGaussianTime",
        multiParam="oimParamMultiWl",
        cosTime="oimParamCosineTime",
        polyWl="oimParamPolynomialWl",
        polyTime="oimParamPolynomialTime",
        powerlawWl="oimParamPowerLawWl",
        powerlawTime="oimParamPowerLawTime",
        rangeWl="oimParamLinearRangeWl",
        templateWl="oimParamLinearTemplateWl",
        tempWl="oimParamLinearTemperatureWl",
        starWl="oimParamLinearStarWl"
)


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
