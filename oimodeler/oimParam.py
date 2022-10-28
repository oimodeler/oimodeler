# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 11:21:58 2022

@author: Ame
"""

import numpy as np
from astropy import units as units
import numbers
from scipy.interpolate import interp1d

#import oimodeler as oim

###############################################################################


class oimParam(object):
    """
    Class of model parameters
    """

    def __init__(self, name=None, value=None, mini=-1*np.inf, maxi=np.inf,
                 description="", unit=1, free=True, error=0):
        """
        Create and initiliaze a new instance of the oimParam class
        Parameters
        ----------
        name : string, optional
            name of the Parameter. The default is None.
        value : float, optional
            value of the parameter. The default is None.
        mini : float, optional
            mininum value allowed for the parameter. The default is -1*np.inf.
        maxi : float, optional
            maximum value allowed for the parameter. The default is np.inf.
        description : string, optional
            A description of the parameter. The default is "".
        unit : 1 or astropy.unit, optional
            unit of the parameter. The default is 1.

        Returns
        -------
        None.

        """
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
        """ The call function will be useful for wavelength or time dependent
        parameters. In a simple oimParam it only return the parameter value
        """
        return self.value

    def __str__(self):
        try:
            return "oimParam {} = {} \xB1 {} {} range=[{},{}] {} ".format(self.name,
                                                                          self.value, self.error, self.unit.to_string(), self.min, self.max, 'free' if self.free else 'fixed')
        except:
            return "oimParam is {}".format(type(self))

    def __repr__(self):
        try:
            return "oimParam at {} : {}={} \xB1 {} {} range=[{},{}] free={} ".format(hex(id(self)), self.name,
                                                                                     self.value, self.error, self.unit.to_string(), self.min, self.max, self.free)
        except:
            return "oimParam at {} is  {}".format(hex(id(self)), type(self))


###############################################################################

class oimInterpWl(object):
    """
    Structure for creating oimParamInterpWl directly in oimParam defintion
    """

    def __init__(self, wl=[], value=None):
        self.wl = wl
        self.value = value

###############################################################################


class oimParamInterpWl(oimParam):
    def __init__(self, param, interpWl):

        self.name = param.name
        self.description = param.description
        self.unit = param.unit

        value = interpWl.value
        wl = interpWl.wl
        nwl = len(wl)

        if value == None:
            value = [self.value]*nwl
        elif isinstance(value, numbers.Number):
            value = [value]*nwl
        else:
            if len(value) != nwl:
                raise TypeError("wl and val should have the same length :"
                                "len(x)={}, len(y)={}".format(len(wl), len(value)))

        self.params = []
        self._wl = wl
        self._nwl = nwl

        for i in range(nwl):

            pi = oimParam(name=param.name, value=value[i], mini=param.min,
                          maxi=param.max, description=param.description,
                          unit=param.unit, free=param.free, error=param.error)
            self.params.append(pi)

        self.value = self.params

    def __call__(self, wl=None, t=None):
        values = np.array([pi.value for pi in self.params])
        return np.interp(wl, self.wl, values, left=values[0], right=values[-1])

    @property
    def wl(self):
        return self._wl

    @wl.setter
    def wl(self, _wl):
        nwl = len(_wl)
        if nwl == self._nwl:
            self._wl = np.array(_wl)
        else:
            raise TypeError("Can't modify number of key wls in oimParamInterpWl "
                            "after creation. Consider creating a new parameter")
###############################################################################


class oimInterpTime(object):
    """
    Structure for creating oimParamInterpTime directly in oimParam defintion
    """

    def __init__(self, t=[], value=None):
        self.t = t
        self.value = value
###############################################################################


class oimParamInterpTime(oimParam):
    def __init__(self, param, interpTime):

        self.name = param.name
        self.description = param.description
        self.unit = param.unit

        value = interpTime.value
        t = interpTime.t
        nt = len(t)

        if value == None:
            value = [self.value]*nt
        elif isinstance(value, numbers.Number):
            value = [value]*nt
        else:
            if len(value) != nt:
                raise TypeError("nt and val should have the same length :"
                                "len(x)={}, len(y)={}".format(len(nt), len(value)))

        self.params = []
        self._t = t
        self._nt = nt

        for i in range(nt):

            pi = oimParam(name=param.name, value=value[i], mini=param.min,
                          maxi=param.max, description=param.description,
                          unit=param.unit, free=param.free, error=param.error)
            self.params.append(pi)

        self.value = self.params

    def __call__(self, wl=None, t=None):
        values = np.array([pi.value for pi in self.params])
        return np.interp(t, self.t, values, left=values[0], right=values[-1])

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, _t):
        nt = len(_t)
        if nt == self._nt:
            self._t = np.array(_t)
        else:
            raise TypeError("Can't modify number of key wls in oimParamInterpTime "
                            "after creation. Consider creating a new parameter")


###############################################################################
class oimParamLinker(object):
    def __init__(self, param, operator="add", fact=0):
        self.param = param
        self.fact = fact

        self.op = None
        self._setOperator(operator)
        self.free = False

    @property
    def unit(self):
        return self.param.unit

    def _setOperator(self, operator):
        if operator == "add":
            self.op = self._add
        elif operator == "mult":
            self.op = self._mult

    def _add(self, val):
        return val + self.fact

    def _mult(self, val):
        return val * self.fact

    def __call__(self, wl=None, t=None):
        return self.op(self.param.__call__(wl, t))

###############################################################################


class oimParamNorm(object):
    def __init__(self, params, norm=1):

        if type(params) == list:
            self.params = params
        else:
            self.params = [params]

        self.norm = norm

        self.free = False

    @property
    def unit(self):
        return self.param.unit

    def __call__(self, wl=None, t=None):

        return self.norm - np.sum([p.__call__(wl, t) for p in self.params])


###############################################################################

class oimParamInterpolator(oimParam):
    def __init__(self, param, **kwargs):

        self.name = param.name
        self.description = param.description
        self.unit = param.unit

        self._init(param, **kwargs)
        self.value = self.params

    def _init(self, param, **kwargs):
        pass

    def _interpFunction(self, wl, t):
        return 0

    def __call__(self, wl=None, t=None):
        return self._interpFunction(wl, t)

    @property
    def params(self):
        params0=self._getParams()
        
        params=[]
        for pi in params0:
            if not(pi in params) and not(isinstance(pi, oimParamLinker)):
                params.append(pi)
        return params
        
       

###############################################################################

class oimParamInterpolatorKeyframes(oimParamInterpolator):

    def _init(self, param, dependence="wl", keyframes=[], keyvalues=[], **kwargs):

        self.dependence = dependence

        self.keyframes = []
        self.keyvalues = []

        for kf in keyframes:
            self.keyframes.append(oimParam(**_standardParameters[dependence]))
            self.keyframes[-1].value = kf

        for kv in keyvalues:
            pi = oimParam(name=param.name, value=kv, mini=param.min,
                          maxi=param.max, description=param.description,
                          unit=param.unit, free=param.free, error=param.error)
            self.keyvalues.append(pi)

        self.params.extend(self.keyframes)
        self.params.extend(self.keyvalues)

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        values = np.array([pi() for pi in self.keyvalues])
        keyframes = np.array([pi() for pi in self.keyframes])
        return np.interp(var, keyframes, values, left=values[0], right=values[-1])
    
    

    def _getParams(self):
        params=[]
        params.extend(self.keyframes)
        params.extend(self.keyvalues)
        return params

###############################################################################


class oimParamInterpolatorWl(oimParamInterpolatorKeyframes):
    def _init(self, param, wl=[], values=[], **kwargs):
        super()._init(param, dependence="wl", keyframes=wl, keyvalues=values)

###############################################################################


class oimParamInterpolatorTime(oimParamInterpolatorKeyframes):
    def _init(self, param, mjd=[], values=[], **kwargs):
        super()._init(param, dependence="mjd", keyframes=mjd, keyvalues=values)


###############################################################################

class oimParamCosineTime(oimParamInterpolator):

    def _init(self, param, T0=0, P=1, values=[0, 1], x0=None, **kwargs):

        self.assymetric = False
        self.T0 = oimParam(name="T0", value=T0,
                           description="Start", unit=units.day)
        self.P = oimParam(name="P", value=P,
                          description="Period", unit=units.day)

        self.values = []

        for kv in values:
            pi = oimParam(name=param.name, value=kv, mini=param.min,
                          maxi=param.max, description=param.description,
                          unit=param.unit, free=param.free, error=param.error)
            self.values.append(pi)

        #self.params.append(self.T0)
        #self.params.append(self.P)
        if x0 != None:
            self.x0 = oimParam(name="x0", value=x0,
                               description="Inflection point", unit=units.one)
            #self.params.append(self.x0)
            self.assymetric = True

        #self.params.extend(self.values)

    def _interpFunction(self, wl, t):
        normt = np.divmod((t-self.T0.value)/self.P.value, 1)[1]
        if self.assymetric == True:
            normt = interp1d([0, self.x0(), 1], [
                             0, 0.5, 1], kind='slinear')(normt)

        return (np.cos(normt*2*np.pi)+1)/2*(self.values[0]() -
                                            self.values[1]())+(self.values[1]())

    def _getParams(self):
        params=[]
        params.extend([self.T0,self.P])
        try:
            params.append(self.x0)
        except:
            pass
        params.extend(self.values)
        return params

###############################################################################

class oimParamGaussian(oimParamInterpolator):

    def _init(self, param, dependence="wl", offset=0, height=0, x0=0, fwhm=0, **kwargs):

        self.dependence = dependence

        self.x0 = oimParam(**_standardParameters[dependence])
        self.x0.name = "x0"
        self.x0.description = "x0"
        self.x0.value = x0
        self.fwhm = oimParam(**_standardParameters[dependence])
        self.fwhm.name = "fwhm"
        self.fwhm.description = "fwhm"
        self.fwhm.value = fwhm

        self.offset = oimParam(name=param.name, value=offset, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)

        self.height = oimParam(name=param.name, value=height, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)

        self.params.extend([self.x0, self.fwhm, self.offset, self.height])

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t

        return self.offset.value+self.height.value *\
            np.exp(-2.77*(var-self.x0())**2/self.fwhm()**2)
            
    def _getParams(self):
        return [self.x0, self.fwhm, self.offset, self.height]

###############################################################################


class oimParamMultipleGaussian(oimParamInterpolator):

    def _init(self, param, dependence="wl", offset=0, height=[], x0=[], fwhm=[], **kwargs):

        try:
            if len(height) != len(x0) and len(height) != len(fwhm):
                raise TypeError(
                    "height, x0 and fwhm should have the same length")
        except:
            print("height, x0 and fwhm should have the same length")

        n = len(height)
        self.x0 = []
        self.fwhm = []
        self.height = []
        for i in range(n):
            self.x0.append(oimParam(**_standardParameters[dependence]))
            self.x0[-1].name = "x0"
            self.x0[-1].description = "x0"
            self.x0[-1].value = x0[i]
            self.fwhm.append(oimParam(**_standardParameters[dependence]))
            self.fwhm[-1].name = "fwhm"
            self.fwhm[-1].description = "fwhm"
            self.fwhm[-1].value = fwhm[i]
            self.height.append(oimParam(name=param.name, value=height[i], mini=param.min,
                                        maxi=param.max, description=param.description,
                                        unit=param.unit, free=param.free, error=param.error))

        self.dependence = dependence

        self.offset = oimParam(name=param.name, value=offset, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)

        self.params.append(self.offset)
        self.params.extend(self.x0)
        self.params.extend(self.fwhm)
        self.params.extend(self.height)

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        val=self.offset.value
        for i in range(len(self.x0)):
            val+=self.height[i]() *np.exp(-2.77*(var-self.x0[i]())**2/self.fwhm[i]()**2)
        return val
    
    def _getParams(self):
        params=[]
        params.append(self.offset)
        params.extend(self.x0)
        params.extend(self.fwhm)
        params.extend(self.height)
        return params
    
###############################################################################


class oimParamGaussianWl(oimParamGaussian):
    def _init(self, param, offset=0, height=0, x0=0, fwhm=0, **kwargs):
        super()._init(param, dependence="wl", offset=offset, height=height, x0=x0, fwhm=fwhm)

###############################################################################


class oimParamGaussianTime(oimParamGaussian):
    def _init(self, param, mjd=[], values=[], **kwargs):
        super()._init(param, dependence="mjd", keyframes=mjd, keyvalues=values)


###############################################################################


# Here is a list of standard parameters to be used when defining new components
_standardParameters = {
    "x": {"name": "x", "value": 0, "description": "x position", "unit": units.mas, "free": False},
    "y": {"name": "y", "value": 0, "description": "y position", "unit": units.mas, "free": False},
    "f": {"name": "f", "value": 1, "description": "flux", "unit": units.one},
    "fwhm": {"name": "fwhm", "value": 0, "description": "FWHM", "unit": units.mas},
    "d": {"name": "d", "value": 0, "description": "Diameter", "unit": units.mas},
    "din": {"name": "din", "value": 0, "description": "Inner Diameter", "unit": units.mas},
    "dout": {"name": "dout", "value": 0, "description": "Outer Diameter", "unit": units.mas},
    "elong": {"name": "elong", "value": 1, "description": "Elongation Ratio", "unit": units.one},
    "pa": {"name": "pa", "value": 0, "description": "Major-axis Position angle", "unit": units.deg},
    "skw": {"name": "skw", "value": 0, "description": "Skewedness", "unit": units.one},
    "skwPa": {"name": "skwPa", "value": 0, "description": "Skewedness Position angle", "unit": units.deg},
    "pixSize": {"name": "pixSize", "value": 0.1, "description": "Pixel Size", "unit": units.mas},
    "dim": {"name": "dim", "value": 128, "description": "Dimension in pixel", "unit": units.one},
    "wl": {"name": "wl", "value": 0, "description": "Wavelength", "unit": units.m, "mini": 0},
    "mjd": {"name": "mjd", "value": 0, "description": "MJD", "unit": units.day}
}
