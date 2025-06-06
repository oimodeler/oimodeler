# -*- coding: utf-8 -*-
"""The oimParam.py module contains the definition of main model parameter class
:func:`oimParam <oimodeler.oimParam.oimParam>`, as well as parameter linkers,
normalizers and interpolators.
"""
import operator
import sys
from functools import reduce
from pathlib import Path
from typing import Any, Dict, List, Union

import astropy.units as u
import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import interp1d

from .oimOptions import constants as const
from .oimUtils import blackbody, linear_to_angular, load_toml

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
    value : int or float, optional
        Value of the parameter. The default is None.
    mini : int or float, optional
        Mininum value allowed for the parameter. The default is -1*np.inf.
    maxi : int or float, optional
        maximum value allowed for the parameter. The default is np.inf.
    description : string, optional
        Description of the parameter. The default is "".
    unit : astropy.unit.Quantity, optional
        Unit of the parameter. The default is astropy.units.one
    free : bool, optional
        Determines if the parameter is to be fitted. The default is None
    error : int or float, optional
        The error of the parameter. The default is 0.
    """

    def __init__(
        self,
        name: Union[str, None] = None,
        value: Union[int, float, None] = None,
        mini: Union[int, float] = -np.inf,
        maxi: Union[int, float] = np.inf,
        description: str = "",
        unit: u.Quantity = u.one,
        free: bool = True,
        error: Union[int, float] = 0,
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
            return "oimParam {} = {} \xb1 {} {} range=[{},{}] {} ".format(
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
            return "oimParam at {} : {}={} \xb1 {} {} range=[{},{}] free={} ".format(
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
        fact: Union[
            int, float, oimParam, List[Union[int, float, oimParam]]
        ] = 0,
    ) -> None:
        """

        Parameters
        ----------
        param : .oimParam
            the oimParam to link with.
        operator : str, optional
            the operator to use. All python operators are available (case-insensitive) either spelt out like
            "add" or with the symbol like "+". The default is "add".
        fact : list of int, float, or .oimParam or int, float, or .oimParam, optional
            The value used for the operation. Can be a list or a single value of float or an oimParam.
            The default is 0.
        """

        self.param = param
        self.fact = (
            fact if isinstance(fact, (tuple, list, np.ndarray)) else [fact]
        )
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
        else:
            self.params = params

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
    interpdescription = "Generic interpolator (to be subclassed)"
    interpparams = []

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
            "Warning set method is not defined for all interpolators.\n"
            "Consider using getFreeParameters method instead"
        )
        params = self.params
        for pi in params:
            for key, value in kwargs.items():
                try:
                    pi.__dict__[key] = value
                except NameError:
                    print("Not valid parameter : {}".format(value))


class oimParamInterpolatorKeyframes(oimParamInterpolator):
    interpdescription = "Interpolation between keyframes in wl or mjd"
    interparams = [
        "keyframes",
        "keyvalues",
        "kind",
        "dependence",
        "fixedRef",
        "extrapolate",
    ]

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


class oimParamInterpolatorWl(oimParamInterpolatorKeyframes):
    interpdescription = "Interpolation between keyframes in wl"
    interparams = ["wl", "values", "kind", "fixedRef", "extrapolate"]

    def _init(self, param, wl=[], values=[], **kwargs):
        super()._init(
            param, dependence="wl", keyframes=wl, keyvalues=values, **kwargs
        )


class oimParamInterpolatorTime(oimParamInterpolatorKeyframes):
    interpdescription = "Interpolation between keyframes in mjd"
    interparams = ["mjd", "values", "kind", "fixedRef", "extrapolate"]

    def _init(self, param, mjd=[], values=[], **kwargs):
        super()._init(
            param, dependence="mjd", keyframes=mjd, keyvalues=values, **kwargs
        )


class oimParamCosineTime(oimParamInterpolator):
    interpdescription = "Assymetrical cosine time interpolator"
    interparams = ["T0", "P", "values", "x0"]

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
    interpdescription = "Generic gaussian interpolator in wl or mjd"
    interparams = ["dependence", "val0", "value", "x0", "fwhm"]

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
    interpdescription = " Gaussian interpolator in wl"
    interparams = ["val0", "value", "x0", "fwhm"]

    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(
            param, dependence="wl", val0=val0, value=value, x0=x0, fwhm=fwhm
        )


class oimParamGaussianTime(oimParamGaussian):
    interpdescription = "Gaussian interpolator in mjd"
    interparams = ["val0", "value", "x0", "fwhm"]

    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(
            param, dependence="mjd", val0=val0, value=value, x0=x0, fwhm=fwhm
        )


class oimParamMultipleGaussian(oimParamInterpolator):
    interpdescription = "Generic multiple gaussian interpolator in wl or mjd"
    interparams = ["dependence", "val0", "values", "x0", "fwhm"]

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
    interpdescription = "Multiple Gaussian interpolator in wl"
    interparams = ["val0", "values", "x0", "fwhm"]

    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(
            param, dependence="wl", val0=val0, values=values, x0=x0, fwhm=fwhm
        )


class oimParamMultipleGaussianTime(oimParamMultipleGaussian):
    interpdescription = "Multiple Gaussian interpolator in mjd"
    interparams = ["val0", "values", "x0", "fwhm"]

    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(
            param, dependence="mjd", val0=val0, values=values, x0=x0, fwhm=fwhm
        )


class oimParamPolynomial(oimParamInterpolator):
    interpdescription = "Generic polynomial interpolation in wl or mjd"
    interparams = ["dependence", "order", "coeffs", "x0"]

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
    interpdescription = "Polynomial interpolation in wl"
    interparams = ["order", "coeffs", "x0"]

    def _init(self, param, order=2, coeffs=None, x0=None, **kwargs):
        super()._init(
            param, dependence="wl", order=order, coeffs=coeffs, x0=x0
        )


class oimParamPolynomialTime(oimParamPolynomial):
    interpdescription = "Polynomial interpolation in mjd"
    interparams = ["order", "coeffs", "x0"]

    def _init(self, param, order=2, coeffs=None, x0=None):
        super()._init(
            param, dependence="mjd", order=order, coeffs=coeffs, x0=x0
        )


class oimParamPowerLaw(oimParamInterpolator):
    interpdescription = "Powerlaw interpolation in wl or mjd"
    """Power-law interpolation, i.e. A*(x/x0)**p."""
    interparams = ["dependence", "x0", "A", "p"]

    def _init(self, param, dependence="wl", x0=0, A=0, p=1, **kwargs):
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
    interpdescription = "Powerlaw interpolation in mjd"
    interparams = ["x0", "A", "p"]

    def _init(self, param, x0=0, A=0, p=1):
        super()._init(param, dependence="mjd", x0=x0, A=A, p=p)


class oimParamPowerLawWl(oimParamPowerLaw):
    interpdescription = "Powerlaw interpolation in wl"
    interparams = ["x0", "A", "p"]

    def _init(self, param, x0=0, A=0, p=1):
        super()._init(param, dependence="wl", x0=x0, A=A, p=p)


class oimParamLinearRangeWl(oimParamInterpolator):
    interpdescription = "Linear range interpolation in wl"
    interparams = ["wlmin", "wlmax", "values"]

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
    interpdescription = (
        "Interpolation in wl using external regular grid template"
    )
    interparams = ["wl0", "dwl", "f_contrib", "values", "kind"]

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
    """Calculates a temperature and wavelength dependent blackbody
    distribution via Planck's law.

    Parameters
    ----------
    param : .oimParam
        The parameter that is to be calculated (interpolated).
    temperature : int or float
        The temperature (K).
    solid_angle : int, float, or .oimParam
        The solid angle of the object or an oimParam containing the solid angle (mas).

    Attributes
    ----------
    temp : .oimParam
    solid_angle : int, float or .oimParam
    """

    interpdescription = "Blackbody in wl for given temperature"
    interparams = ["temp", "solid_angle"]

    def _init(
        self,
        param: oimParam,
        temp: Union[int, float] = 0,
        solid_angle: Union[int, float, oimParam] = 0,
        **kwargs,
    ) -> None:
        """The subclass's constructor."""
        self.temp = oimParam(
            name="temp",
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
            Wavelengths (m).
        t : numpy.ndarray
            Times (mjd).

        Returns
        -------
        blackbody_distribution : numpy.ndarray
            The star's flux (Jy).
        """
        if isinstance(self.solid_angle, oimParam):
            solid_angle = self.solid_angle(wl, t)
        else:
            solid_angle = self.solid_angle

        return (
            blackbody(self.temp(wl, t), const.c / wl)
            / u.rad.to(u.mas) ** 2
            * solid_angle**2
            * 1e23
        )


# TODO: Add extinction to this model at some point
class oimParamLinearStarWl(oimParamInterpolator):
    """Calculates the stellar flux for different wavelengths.

    Calculates the stellar radius none is provided from the luminosity.

    Parameters
    ----------
    param : oimParam
        The parameter that is to be calculated (interpolated).
    temp : array_like
        The temperature (K).
    dist : int or float
        The distance to the star (pc).
    lum : int or float, optional
        The star's luminosity (Lsun).
    radius : int or float, optional
        The star's radius (Rsun).

    Attributes
    ----------
    stellar_radius : astropy.units.m
        The stellar radius (m).
    stellar_angular_radius : astropy.units.mas
        The angular stellar radius (mas).
    dist : .oimParam
        An oimParam containing the distance to the star (pc).
    lum : .oimParam
        An oimParam containing the star's luminosity (Lsun).
    """

    interpdescription = (
        "Blackbody in wl for given Teff, dist, and radius or lum"
    )
    interparams = ["temp", "dist", "lum", "radius"]

    def _init(
        self,
        param: oimParam,
        temp: Union[int, float, ArrayLike] = 0,
        dist: Union[int, float] = 0,
        lum: Union[int, float, None] = 0,
        radius: Union[int, float, None] = None,
        **kwargs: Dict,
    ) -> None:
        self._angular_stellar_radius = None
        self.temp = oimParam(
            name="temp",
            value=temp,
            unit=u.K,
            free=False,
            description="The star's effective temperature",
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
            name="radius",
            value=radius,
            unit=u.R_sun,
            free=False,
            description="The star's radius",
        )
        if self.lum.value is None and self.radius.value is None:
            raise ValueError(
                "Either luminosity or radius must be provided to calculate stellar flux."
            )

    @property
    def angular_stellar_radius(self) -> float:
        """Calculates the angular stellar radius.

        Returns
        -------
        angular_stellar_radius : float
            The star's radius (mas).
        """
        if self._angular_stellar_radius is not None:
            return self._angular_stellar_radius

        if self.radius.value is not None:
            stellar_radius = self.radius.value * self.radius.unit.to(u.au)
        else:
            luminosity = self.lum.value * self.lum.unit.to(u.W)
            stellar_radius = np.sqrt(
                luminosity / (4 * np.pi * const.sigma_sb * self.temp.value**4)
            ) * u.m.to(u.au)

        self._angular_stellar_radius = (
            linear_to_angular(stellar_radius, self.dist.value) * 1e3
        )
        return self._angular_stellar_radius

    def _getParams(self):
        """Gets the parameters of the interpolator."""
        return [self.temp]

    def _interpFunction(self, wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates the stellar flux from its distance and radius at
        the specified wavelengths.

        Parameters
        ----------
        wl : numpy.ndarray
            Wavelengths (m).
        t : numpy.ndarray
            Times (mjd).

        Returns
        -------
        stellar_flux : np.ndarray
            The star's flux (Jy).
        """
        return (
            blackbody(self.temp(wl, t), const.c / wl)
            / u.rad.to(u.mas) ** 2
            * np.pi
            * self.angular_stellar_radius**2
            * 1e23
        )
