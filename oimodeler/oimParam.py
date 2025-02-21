# -*- coding: utf-8 -*-
"""The oimParam.py module contains the definition of main model parameter class
:func:`oimParam <oimodeler.oimParam.oimParam>`, as well as parameter linkers,
normalizers and interpolators.
"""
import operator
import sys
from functools import reduce
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import astropy.constants as const
import astropy.units as u
import numpy as np
import toml
from astropy.modeling import models
from numpy.typing import ArrayLike
from scipy.interpolate import interp1d


def load_toml(toml_file: Path) -> Dict[str, Any]:
    """Loads a toml file into a dictionary."""
    with open(toml_file, "r") as file:
        dictionary = toml.load(file)

    for value in dictionary.values():
        if "unit" in value:
            if value["unit"] == "one":
                value["unit"] = u.one
            else:
                value["unit"] = u.Unit(value["unit"])

    return dictionary


_standardParameters: Dict[str, Any] = load_toml(
    Path(__file__).parent / "config" / "standard_parameters.toml"
)
_interpolators: Dict[str, Any] = load_toml(
    Path(__file__).parent / "config" / "interpolators.toml"
)["interpolators"]

_operators = {
    "add": operator.add,
    "+": operator.add,
    "sub": operator.sub,
    "-": operator.sub,
    "mul": operator.mul,
    "*": operator.mul,
    "div": operator.truediv,
    "/": operator.truediv,
    "floordiv": operator.floordiv,
    "//": operator.floordiv,
    "mod": operator.mod,
    "%": operator.mod,
    "pow": operator.pow,
    "**": operator.pow,
}


class oimParam:
    """Class of model parameters.

    Parameters
    ----------
    name : string, optional
        Name of the Parameter. The default is None.
    value : float, optional
        Value of the parameter. The default is None.
    mini : float, optional
        Mininum value allowed for the parameter. The default is -1*np.inf.
    maxi : float, optional
        maximum value allowed for the parameter. The default is np.inf.
    description : string, optional
        Description of the parameter. The default is "".
    unit : 1 or astropy.unit, optional
        Unit of the parameter. The default is astropy.units.one
    free : bool, optional
        Determines if the parameter is to be fitted. The default is None
    """

    def __init__(
        self,
        name=None,
        value=None,
        mini=-np.inf,
        maxi=np.inf,
        description="",
        unit=u.one,
        free=True,
        error=0,
    ):
        """Initialize a new instance of the oimParam class."""
        self.name = name
        self.value = value
        self.error = error
        self.min = mini
        self.max = maxi
        self.free = free
        self.description = description
        self.unit = unit

    def set(self, **kwargs):
        for key, value in kwargs.items():
            try:
                self.__dict__[key] = value
            except NameError:
                print("Note valid parameter : {}".format(value))

    def __call__(self, wl=None, t=None):
        """The call function will be useful for wavelength or time dependent
        parameters. In a simple oimParam it only return the parameter value
        """
        return self.value

    def __str__(self):
        """String (print) representation of the oimParam class"""
        try:
            return "oimParam {} = {} \xB1 {} {} range=[{},{}] {} ".format(
                self.name,
                self.value,
                self.error,
                self.unit.to_string(),
                self.min,
                self.max,
                "free" if self.free else "fixed",
            )
        except:
            return "oimParam is {}".format(type(self))

    def __repr__(self):
        """String (console) representation of the oimParam class"""
        try:
            return "oimParam at {} : {}={} \xB1 {} {} range=[{},{}] free={} ".format(
                hex(id(self)),
                self.name,
                self.value,
                self.error,
                self.unit.to_string(),
                self.min,
                self.max,
                self.free,
            )
        except:
            return "oimParam at {} is  {}".format(hex(id(self)), type(self))


class oimParamLinker:
    """Class to directly link two oimParam."""

    def __init__(
        self,
        param: oimParam,
        operator: str = "add",
        fact: Union[float, oimParam, List[Union[float, oimParam]]] = 0,
    ) -> None:
        """

        Parameters
        ----------
        param : oimParam
            the oimParam to link with.
        operator : str, optional
            the operator to use. All python operators are available (case-insensitive) either spelt out like
            "add" or with the symbol like "+". The default is "add".
        fact : list of float or list of oimParam or oimParam or float, optional
            The value used for the operation. Can be a list or a single value of float or an oimParam.
            The default is 0.
        """

        self.param = param
        self.fact = fact if isinstance(fact, (tuple, list, np.ndarray)) else [fact]
        self.op = _operators[operator.lower()]
        self.free = False

    @property
    def unit(self):
        return self.param.unit

    def __call__(self, wl=None, t=None):
        values = [
            self.param(wl, t),
            *[
                val(wl, t) if isinstance(val, oimParam) else val
                for val in self.fact
            ],
        ]
        return reduce(self.op, values)


class oimParamNorm:
    """
    Class to normalize a list of oimParam

    Example :
    p2 = oimParamNorm([p0,p1)]

    The value of p2 will always be 1-p2-p1
    """

    def __init__(self, params, norm=1):
        if not isinstance(params, list):
            self.params = [params]

        self.norm = norm
        self.free = False

    @property
    def unit(self):
        return self.params.unit

    def __call__(self, wl=None, t=None):

        res = self.norm
        for p in self.params:
            res -= p(wl, t)
        return res


class oimInterp:
    """Macro to directly create an oimParamInterpolator

    Example :

    .. code-block:: python

        g = oim.oimGauss(fwhm=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 8]))

    This will create a gaussiand component and remplace the parameter fwhm by an
    instance of the oimParamInterpolatorWl class (i.e., the wavelength linear
    interpolator. The custom interpolator the reference to the interpolator
    need to be added to the  :func:`_interpolators <oimodeler.oimParam._interpolators>`
    dictionnary.

    Parameters
    ----------
    name : str
        Keyname for the interpolators registered in the _interpolators
        dictionary.
    **kwargs : dict
        Parameters from the create oimParamInterpolator-derived class.

    Attributes
    ----------
    kwargs : dict
        Parameters from the create oimParamInterpolator-derived class.
    type : oimParamInterpolator
        A param interpolator contained in the _interpolators dictionary.
        For the local definition in the `oimParam` module, strings can be used.
        To redefine the dictionary elements from outside, use the class
        variables, otherwise the local definition will be used.
    """

    def __init__(self, name, **kwargs):
        self.kwargs = kwargs
        self.type = _interpolators[name]
        if isinstance(self.type, str):
            self.type = getattr(sys.modules[__name__], self.type)


class oimParamInterpolator(oimParam):
    def __init__(self, param, **kwargs):
        self.name = param.name
        self.description = param.description
        self.unit = param.unit
        self.param0 = param
        self._init(param, **kwargs)

    def _init(self, param, **kwargs):
        pass

    def _interpFunction(self, wl, t):
        return 0

    def __call__(self, wl=None, t=None):
        return self._interpFunction(wl, t)

    def _getParams(self):
        pass

    @property
    def params(self):
        params0 = self._getParams()

        params = []
        for pi in params0:
            if pi not in params and not isinstance(pi, oimParamLinker):
                params.append(pi)
        return params

    def set(self, **kwargs):
        print(
            "Warning set method is not defined for all interpolators. use getFreeParametezrs instead"
        )
        for _ in self.params:
            for key, value in kwargs.items():
                try:
                    self.__dict__[key] = value
                except NameError:
                    print("Not valid parameter : {}".format(value))


class oimParamInterpolatorKeyframes(oimParamInterpolator):
    def _init(
        self,
        param,
        dependence="wl",
        keyframes=[],
        keyvalues=[],
        kind="linear",
        fixedRef=True,
        extrapolate=False,
        **kwargs,
    ):
        self.dependence = dependence
        self.fixedRef = fixedRef
        self.kind = kind
        self.extrapolate = extrapolate

        self.keyframes = []
        self.keyvalues = []

        for kf in keyframes:
            self.keyframes.append(oimParam(**_standardParameters[dependence]))
            self.keyframes[-1].value = kf

        for kv in keyvalues:
            pi = oimParam(
                name=param.name,
                value=kv,
                mini=param.min,
                maxi=param.max,
                description=param.description,
                unit=param.unit,
                free=param.free,
                error=param.error,
            )
            self.keyvalues.append(pi)

    def _interpFunction(self, wl, t):
        if self.dependence == "wl":
            var = wl
        else:
            var = t
        values = np.array([pi() for pi in self.keyvalues])
        keyframes = np.array([pi() for pi in self.keyframes])

        if self.extrapolate:
            fill_value = "extrapolate"
            bounds_error = None
        else:
            fill_value = (values[0], values[-1])
            bounds_error = False

        return interp1d(
            keyframes,
            values,
            fill_value=fill_value,
            kind=self.kind,
            bounds_error=bounds_error,
        )(var)

    def _getParams(self):
        params = []
        if not self.fixedRef:
            params.extend(self.keyframes)

        params.extend(self.keyvalues)
        return params


class oimParamMultiWl(oimParamInterpolatorKeyframes):
    """Class of model parameters for multiple wavelengths.

    Notes
    -----
    The difference of this class to the interpolators is that
    there will be no interpolation done, but the parameters will be
    returned for the different wavelengths, at which they are specified.
    Will raise an error if there are missing wavelengths.
    """

    def _init(self, param, wl=[], values=[], **kwargs):
        super()._init(
            param, dependence="wl", keyframes=wl, keyvalues=values, **kwargs
        )

    def _interpFunction(self, wl, t):
        """Returns the parameter's value for the given wavelength."""
        keyframes = np.array([param.value for param in self.keyframes])
        values = np.array([param.value for param in self.keyvalues])
        new_values = np.zeros(wl.shape)

        for index, wavelength in enumerate(wl):
            if wavelength in keyframes:
                value_index = np.where(wavelength == keyframes)[0][0]
                if len(wl.shape) == 1:
                    new_values[index] = values[value_index]
                else:
                    new_values[:, index, :, :] = values[value_index]

        return new_values


class oimParamInterpolatorWl(oimParamInterpolatorKeyframes):
    def _init(self, param, wl=[], values=[], **kwargs):
        super()._init(
            param, dependence="wl", keyframes=wl, keyvalues=values, **kwargs
        )


class oimParamInterpolatorTime(oimParamInterpolatorKeyframes):
    def _init(self, param, mjd=[], values=[], **kwargs):
        super()._init(
            param, dependence="mjd", keyframes=mjd, keyvalues=values, **kwargs
        )


class oimParamCosineTime(oimParamInterpolator):
    def _init(self, param, T0=0, P=1, values=[0, 1], x0=None, **kwargs):
        self.asymmetric = False
        self.T0 = oimParam(
            name="T0", value=T0, description="Start", unit=u.day
        )
        self.P = oimParam(name="P", value=P, description="Period", unit=u.day)

        self.values = []

        for kv in values:
            pi = oimParam(
                name=param.name,
                value=kv,
                mini=param.min,
                maxi=param.max,
                description=param.description,
                unit=param.unit,
                free=param.free,
                error=param.error,
            )
            self.values.append(pi)

        if x0 is not None:
            self.x0 = oimParam(
                name="x0", value=x0, description="Inflection point", unit=u.one
            )
            self.asymmetric = True

    def _interpFunction(self, wl, t):
        normt = np.divmod((t - self.T0.value) / self.P.value, 1)[1]
        if self.asymmetric:
            normt = interp1d([0, self.x0(), 1], [0, 0.5, 1], kind="slinear")(
                normt
            )

        return (np.cos(normt * 2 * np.pi) + 1) / 2 * (
            self.values[0]() - self.values[1]()
        ) + (self.values[1]())

    def _getParams(self):
        params = []
        params.extend([self.T0, self.P])
        try:
            params.append(self.x0)
        except:
            pass
        params.extend(self.values)
        return params


class oimParamGaussian(oimParamInterpolator):
    def _init(
        self, param, dependence="wl", val0=0, value=0, x0=0, fwhm=0, **kwargs
    ):
        self.dependence = dependence
        self.x0 = oimParam(**_standardParameters[dependence])
        self.x0.name = "x0"
        self.x0.description = "x0"
        self.x0.value = x0
        self.fwhm = oimParam(**_standardParameters[dependence])
        self.fwhm.name = "fwhm"
        self.fwhm.description = "fwhm"
        self.fwhm.value = fwhm

        self.val0 = oimParam(
            name=param.name,
            value=val0,
            mini=param.min,
            maxi=param.max,
            description=param.description,
            unit=param.unit,
            free=param.free,
            error=param.error,
        )

        self.value = oimParam(
            name=param.name,
            value=value,
            mini=param.min,
            maxi=param.max,
            description=param.description,
            unit=param.unit,
            free=param.free,
            error=param.error,
        )
        self.value.value

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t

        return self.val0() + (self.value() - self.val0()) * np.exp(
            -2.77 * (var - self.x0()) ** 2 / self.fwhm() ** 2
        )

    def _getParams(self):
        return [self.x0, self.fwhm, self.val0, self.value]


class oimParamGaussianWl(oimParamGaussian):
    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(
            param, dependence="wl", val0=val0, value=value, x0=x0, fwhm=fwhm
        )


class oimParamGaussianTime(oimParamGaussian):
    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(
            param, dependence="mjd", val0=val0, value=value, x0=x0, fwhm=fwhm
        )


class oimParamMultipleGaussian(oimParamInterpolator):

    def _init(
        self,
        param,
        dependence="wl",
        val0=0,
        values=[],
        x0=[],
        fwhm=[],
        **kwargs,
    ):
        try:
            if len(values) != len(x0) and len(values) != len(fwhm):
                raise TypeError(
                    "values, x0 and fwhm should have the same length"
                )
        except:
            print("values, x0 and fwhm should have the same length")

        n = len(values)
        self.x0 = []
        self.fwhm = []
        self.values = []
        for i in range(n):
            self.x0.append(oimParam(**_standardParameters[dependence]))
            self.x0[-1].name = "x0"
            self.x0[-1].description = "x0"
            self.x0[-1].value = x0[i]
            self.fwhm.append(oimParam(**_standardParameters[dependence]))
            self.fwhm[-1].name = "fwhm"
            self.fwhm[-1].description = "fwhm"
            self.fwhm[-1].value = fwhm[i]
            self.values.append(
                oimParam(
                    name=param.name,
                    value=values[i],
                    mini=param.min,
                    maxi=param.max,
                    description=param.description,
                    unit=param.unit,
                    free=param.free,
                    error=param.error,
                )
            )

        self.dependence = dependence

        self.val0 = oimParam(
            name=param.name,
            value=val0,
            mini=param.min,
            maxi=param.max,
            description=param.description,
            unit=param.unit,
            free=param.free,
            error=param.error,
        )

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        val = self.val0()
        for i in range(len(self.x0)):
            val += (self.values[i]() - self.val0()) * np.exp(
                -2.77 * (var - self.x0[i]()) ** 2 / self.fwhm[i]() ** 2
            )
        return val

    def _getParams(self):
        params = []
        params.append(self.val0)
        params.extend(self.x0)
        params.extend(self.fwhm)
        params.extend(self.values)
        return params


class oimParamMultipleGaussianWl(oimParamMultipleGaussian):
    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(
            param, dependence="wl", val0=val0, values=values, x0=x0, fwhm=fwhm
        )


class oimParamMultipleGaussianTime(oimParamMultipleGaussian):
    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(
            param, dependence="mjd", val0=val0, values=values, x0=x0, fwhm=fwhm
        )


class oimParamPolynomial(oimParamInterpolator):
    def _init(
        self, param, dependence="wl", order=2, coeffs=None, x0=None, **kwargs
    ):
        self.dependence = dependence

        if x0 is None:
            self.x0 = 0
        else:
            self.x0 = x0

        if coeffs is None:
            coeffs = [0] * (order + 1)

        self.coeffs = []
        for ci in coeffs:
            pi = oimParam(
                name=param.name,
                value=ci,
                mini=param.min,
                maxi=param.max,
                description=param.description,
                unit=param.unit,
                free=param.free,
                error=param.error,
            )
            self.coeffs.append(pi)

        self.params.extend(self.coeffs)
        if not (x0 is None):
            self.params.append(x0)

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        var = var - self.x0
        c = np.flip([ci() for ci in self.coeffs])
        return np.poly1d(c)(var)

    def _getParams(self):
        return self.coeffs


class oimParamPolynomialWl(oimParamPolynomial):
    def _init(self, param, order=2, coeffs=None, x0=None, **kwargs):
        super()._init(
            param, dependence="wl", order=order, coeffs=coeffs, x0=x0
        )


class oimParamPolynomialTime(oimParamPolynomial):
    def _init(self, param, order=2, coeffs=None, x0=None):
        super()._init(
            param, dependence="mjd", order=order, coeffs=coeffs, x0=x0
        )


class oimParamPowerLaw(oimParamInterpolator):
    """Power-law interpolation, i.e. A*(x/x0)**p."""

    def _init(self, param, dependence, x0, A, p, **kwargs):
        self.dependence = dependence

        self.x0 = oimParam(**_standardParameters[dependence])
        self.x0.name = "x0"
        self.x0.description = "Power-law reference value (wl0 or t0)"
        self.x0.value = x0
        self.x0.free = False

        self.A = oimParam(**_standardParameters["scale"])
        self.A.name = "A"
        self.A.description = "Power-law scale factor"
        self.A.value = A
        self.A.free = True

        self.p = oimParam(**_standardParameters["index"])
        self.p.name = "p"
        self.p.description = "Power-law index"
        self.p.value = p
        self.p.free = True

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            x = wl
        elif self.dependence == "mjd":
            x = t
        else:
            raise NotImplementedError(
                'No support for interpolation along "%s"' % self.dependence
            )

        return self.A() * (x / self.x0()) ** self.p()

    def _getParams(self):
        return [self.x0, self.A, self.p]


class oimParamPowerLawTime(oimParamPowerLaw):
    def _init(self, param, x0, A, p):
        super()._init(param, dependence="t", x0=x0, A=A, p=p)


class oimParamPowerLawWl(oimParamPowerLaw):
    def _init(self, param, x0, A, p):
        super()._init(param, dependence="wl", x0=x0, A=A, p=p)


class oimParamLinearRangeWl(oimParamInterpolator):
    def _init(
        self, param, wlmin=2e-6, wlmax=3e-6, values=[], kind="linear", **kwargs
    ):
        self.kind = kind

        n = len(values)
        self.wlmin = oimParam(**_standardParameters["wl"])
        self.wlmin.name = "wlmin"
        self.wlmin.description = "Min of wl range"
        self.wlmin.value = wlmin
        self.wlmin.free = False

        self.wlmax = oimParam(**_standardParameters["wl"])
        self.wlmax.name = "wlmax"
        self.wlmax.description = "Max of the wl range"
        self.wlmax.value = wlmax
        self.wlmax.free = False

        self.values = []

        for i in range(n):
            self.values.append(
                oimParam(
                    name=param.name,
                    value=values[i],
                    mini=param.min,
                    maxi=param.max,
                    description=param.description,
                    unit=param.unit,
                    free=param.free,
                    error=param.error,
                )
            )

    def _interpFunction(self, wl, t):

        vals = np.array([vi.value for vi in self.values])
        nwl = vals.size
        wl0 = np.linspace(self.wlmin.value, self.wlmax.value, num=nwl)
        return interp1d(wl0, vals, kind=self.kind, fill_value="extrapolate")(
            wl
        )

    def _getParams(self):
        params = []
        params.extend(self.values)
        params.append(self.wlmin)
        params.append(self.wlmax)
        return params


class oimParamLinearTemplateWl(oimParamInterpolator):
    def _init(
        self,
        param,
        wl0=2e-6,
        dwl=1e-9,
        f_contrib=1.0,
        values=[],
        kind="linear",
        **kwargs,
    ):
        self.kind = kind

        n = len(values)
        self.wl0 = oimParam(**_standardParameters["wl"])
        self.wl0.name = "wl0"
        self.wl0.description = "Initial wl of the range"
        self.wl0.value = wl0
        self.wl0.free = False

        self.dwl = oimParam(**_standardParameters["wl"])
        self.dwl.name = "dwl"
        self.dwl.description = "wl step in range"
        self.dwl.value = dwl
        self.dwl.free = False

        self.f_contrib = oimParam(**_standardParameters["f"])
        self.f_contrib.name = "f_contrib"
        self.f_contrib.description = (
            "Flux contribution at the wavelength"
            " corresponding to the maximum of the"
            " associated emission template"
        )
        self.f_contrib.value = f_contrib
        self.f_contrib.free = False

        self.values = []

        for i in range(n):
            self.values.append(
                oimParam(
                    name=param.name,
                    value=values[i],
                    mini=param.min,
                    maxi=param.max,
                    description=param.description,
                    unit=param.unit,
                    free=param.free,
                    error=param.error,
                )
            )

    def _interpFunction(self, wl, t):
        vals = np.array([vi.value for vi in self.values])
        nwl = vals.size
        wl0 = np.linspace(
            self.wl0.value, self.wl0.value + self.dwl.value * nwl, num=nwl
        )
        interp_template = interp1d(
            wl0, vals, kind=self.kind, fill_value="extrapolate"
        )(wl)
        interp_template_norm = interp_template / interp_template.max()
        return interp_template_norm * self.f_contrib.value

    def _getParams(self):
        params = []
        params.append(self.f_contrib)
        return params


class oimParamLinearTemperatureWl(oimParamInterpolatorKeyframes):
    """Interpolates/Calculates a blackbody distribution
    for different wavelengths for the input temperatures
    (via Planck's law).

    Parameters
    ----------
    param : oimParam
        The parameter that is to be calculated/interpolated.
    temperature : int or float or numpy.ndarray
        The temperature(s) to be calculated.

    Attributes
    ----------
    temp : oimParam
    """

    def _init(
        self,
        param: oimParam,
        temp: Union[int, float, ArrayLike],
        solid_angle: Union[int, float, ArrayLike, oimParam],
        **kwargs,
    ) -> None:
        """The subclass's constructor."""
        self.temp = oimParam(
            name="T",
            value=temp,
            unit=u.K,
            free=True,
            mini=0,
            maxi=3000,
            description="The temperature",
        )
        self.solid_angle = solid_angle

    def _getParams(self):
        """Gets the parameters of the interpolator."""
        return [self.temp]

    def _interpFunction(self, wl: np.ndarray, t: np.ndarray):
        """Calculates a temperature and wavelength dependent blackbody
        distribution via Planck's law.

        Parameters
        ----------
        wl : numpy.ndarray
            Wavelengths [m].
        t : numpy.ndarray
            Times.

        Returns
        -------
        blackbody_distribution : astropy.units.Jy
            The star's flux.
        """
        bb = models.BlackBody(temperature=self.temp(wl, t))
        return bb(wl * u.m).to(u.erg / u.cm**2 / u.Hz / u.s / u.mas**2).value


class oimParamLinearStarWl(oimParamInterpolator):
    """Interpolates/Calculates the stellar flux for distinct wavelengths.

    Parameters
    ----------
    param : oimParam
        The parameter that is to be calculated/interpolated.
    temperature : array_like
        The temperature distribution to be calculated at different wavelengths.
    distance : int or float
        The distance to the star [pc].
    luminostiy : int or float
        The star's luminosity [Lsun].
    **kwargs : dict

    Attributes
    ----------
    stellar_radius : astropy.units.m
    stellar_angular_radius : astropy.units.mas
    dist : oimParam
        An oimParam containing the distance to the star.
    lum : oimParam
        An oimParam containing the star's luminosity.
    """

    def _init(
        self,
        param: oimParam,
        temp: Union[int, float, ArrayLike],
        dist: Union[int, float],
        lum: Union[int, float],
        radius: Optional[Union[int, float]] = None,
        **kwargs: Dict,
    ) -> None:
        """The subclass's constructor."""
        self._stellar_radius = None
        self._stellar_angular_radius = None
        self.temp = oimParam(
            name="T",
            value=temp,
            unit=u.K,
            free=False,
            description="The Temperature",
        )
        self.dist = oimParam(
            name="dist",
            value=dist,
            unit=u.pc,
            free=False,
            description="Distance to the star",
        )
        self.lum = oimParam(
            name="lum",
            value=lum,
            unit=u.Lsun,
            free=False,
            description="The star's luminosity",
        )
        self.radius = oimParam(
            name="stellar radius",
            value=radius,
            unit=u.R_sun,
            free=False,
            description="The star's radius",
        )

    @property
    def stellar_radius(self) -> u.m:
        """Calculates the stellar radius.

        Returns
        -------
        stellar_radius : astropy.units.m
            The star's radius.
        """
        if self.radius.value is not None:
            return self.radius.value
        if self._stellar_radius is None:
            luminosity = (self.lum.value * self.lum.unit).to(u.W)
            self._stellar_radius = np.sqrt(
                luminosity
                / (
                    4
                    * np.pi
                    * const.sigma_sb
                    * (self.temp.value * self.temp.unit) ** 4
                )
            )
        return self._stellar_radius

    @property
    def stellar_radius_angular(self) -> u.mas:
        r"""Calculates the parallax from the stellar radius and the distance to
        the object.

        Returns
        -------
        stellar_radius_angular : astropy.units.mas
            The parallax of the stellar radius.

        Notes
        -----
        The formula for the angular diameter $ \delta = \frac{d}{D} $ is used.
        This produces an output in radians.
        """
        if self._stellar_angular_radius is None:
            distance = self.dist.value * self.dist.unit
            if self.radius.value is not None:
                self._stellar_angular_radius = (
                    (self.radius.value * self.radius.unit).to(u.m)
                    / distance.to(u.m)
                    * u.rad
                )
            else:
                self._stellar_angular_radius = (
                    self.stellar_radius.to(u.m) / distance.to(u.m) * u.rad
                )
        return self._stellar_angular_radius

    def _getParams(self):
        """Gets the parameters of the interpolator."""
        return [self.temp]

    def _calc_spectral_radiance(self, wl: u.m) -> np.ndarray:
        """Calculates a temperature and wavelength dependent blackbody
        distribution via Planck's law.


        Parameters
        ----------
        wl : astropy.units.m
            Wavelengths [m].

        Returns
        -------
        blackbody_distribution : astropy.units.Jy
            The star's flux.
        """
        plancks_law = models.BlackBody(
            temperature=self.temp.value * self.temp.unit
        )
        return plancks_law(wl).to(u.erg / (u.cm**2 * u.Hz * u.s * u.rad**2))

    def _interpFunction(self, wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates the stellar flux from its distance and radius and
        interpolates it to the input wavelengths,

        This function is not time dependent.

        Parameters
        ----------
        wl : numpy.ndarray
            Wavelengths [m].
        t : numpy.ndarray
            Times.

        Returns
        -------
        stellar_flux : np.ndarray
            The star's flux [Jy].
        """
        spectral_radiance = self._calc_spectral_radiance(wl * u.m)
        spectral_radiance = (
            (np.pi * spectral_radiance * self.stellar_radius_angular**2)
            .to(u.Jy)
            .value
        )
        return spectral_radiance
