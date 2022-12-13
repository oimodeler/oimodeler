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
        
        res = self.norm
        for p in self.params:
            res -= p(wl, t)
        return res


###############################################################################

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

    def _init(self, param, dependence="wl", keyframes=[], keyvalues=[],
              kind="linear",fixedRef=True, extrapolate=False,**kwargs):

        self.dependence = dependence
        self.fixedRef = fixedRef
        self.kind=kind
        self.extrapolate = extrapolate

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
        

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        values = np.array([pi() for pi in self.keyvalues])
        keyframes = np.array([pi() for pi in self.keyframes])
        
        if self.extrapolate==True:
            fill_value="extrapolate"
            bounds_error=None
        else:
            fill_value = (values[0],values[-1])
            bounds_error=False
        #return np.interp(var, keyframes, values, left=values[0], right=values[-1])
        return interp1d( keyframes, values, fill_value=fill_value,
                        kind=self.kind,bounds_error=bounds_error)(var)
    
    

    def _getParams(self):
        params=[]
        if self.fixedRef==False:
            params.extend(self.keyframes)
        params.extend(self.keyvalues)
        return params

###############################################################################


class oimParamInterpolatorWl(oimParamInterpolatorKeyframes):
    def _init(self, param, wl=[], values=[], **kwargs):
        super()._init(param, dependence="wl", keyframes=wl, keyvalues=values,**kwargs)

###############################################################################


class oimParamInterpolatorTime(oimParamInterpolatorKeyframes):
    def _init(self, param, mjd=[], values=[], **kwargs):
        super()._init(param, dependence="mjd", keyframes=mjd, keyvalues=values,**kwargs)


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

    def _init(self, param, dependence="wl", val0=0, value=0, x0=0, fwhm=0, **kwargs):
        self.dependence = dependence
        self.x0 = oimParam(**_standardParameters[dependence])
        self.x0.name = "x0"
        self.x0.description = "x0"
        self.x0.value = x0
        self.fwhm = oimParam(**_standardParameters[dependence])
        self.fwhm.name = "fwhm"
        self.fwhm.description = "fwhm"
        self.fwhm.value = fwhm

        self.val0 = oimParam(name=param.name, value=val0, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)

        self.value = oimParam(name=param.name, value=value, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)
        self.value.value
    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t

        return self.val0()+(self.value()-self.val0()) \
            *np.exp(-2.77*(var-self.x0())**2/self.fwhm()**2)
            
    def _getParams(self):
        return [self.x0, self.fwhm, self.val0, self.value]

    
###############################################################################

class oimParamGaussianWl(oimParamGaussian):
    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(param, dependence="wl", val0=val0, value=value, x0=x0, fwhm=fwhm)

###############################################################################

class oimParamGaussianTime(oimParamGaussian):
    def _init(self, param, val0=0, value=0, x0=0, fwhm=0, **kwargs):
        super()._init(param, dependence="mjd", val0=val0, value=value, x0=x0, fwhm=fwhm)

###############################################################################

class oimParamMultipleGaussian(oimParamInterpolator):

    def _init(self, param, dependence="wl", val0=0, values=[], x0=[], fwhm=[], **kwargs):
        try:
            if len(values) != len(x0) and len(values) != len(fwhm):
                raise TypeError(
                    "values, x0 and fwhm should have the same length")
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
            self.values.append(oimParam(name=param.name, value=values[i], mini=param.min,
                                        maxi=param.max, description=param.description,
                                        unit=param.unit, free=param.free, error=param.error))

        self.dependence = dependence

        self.val0 = oimParam(name=param.name, value=val0, mini=param.min,
                               maxi=param.max, description=param.description,
                               unit=param.unit, free=param.free, error=param.error)


    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        val=self.val0()
        for i in range(len(self.x0)):
            val+=(self.values[i]()-self.val0()) \
                *np.exp(-2.77*(var-self.x0[i]())**2/self.fwhm[i]()**2)
        return val
    
    def _getParams(self):
        params=[]
        params.append(self.val0)
        params.extend(self.x0)
        params.extend(self.fwhm)
        params.extend(self.values)
        return params
###############################################################################

class oimParamMultipleGaussianWl(oimParamMultipleGaussian):
    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(param, dependence="wl", val0=val0, values=values, x0=x0, fwhm=fwhm)

###############################################################################

class oimParamMultipleGaussianTime(oimParamMultipleGaussian):
    def _init(self, param, val0=0, values=[], x0=[], fwhm=[], **kwargs):
        super()._init(param, dependence="mjd", val0=val0, values=values, x0=x0, fwhm=fwhm)
        
###############################################################################

class oimParamPolynomial(oimParamInterpolator):

    def _init(self, param, dependence="wl", order=2,coeffs=None, x0=None,**kwargs):
        self.dependence = dependence
        
        if x0 is None:
            self.x0=0
        else:
            self.x0=x0
        
        if coeffs is None:
            coeffs = [0]*(order+1)
        
        self.coeffs=[]
        for ci in coeffs:
            pi = oimParam(name=param.name, value=ci, mini=param.min,
                          maxi=param.max, description=param.description,
                          unit=param.unit, free=param.free, error=param.error)
            self.coeffs.append(pi)

        self.params.extend(self.coeffs)
        if not(x0 is None):
            self.params.append(x0)

    def _interpFunction(self, wl, t):

        if self.dependence == "wl":
            var = wl
        else:
            var = t
        var=var-self.x0
        c=np.flip([ci() for ci in self.coeffs])
        return np.poly1d(c)(var)
    
    

    def _getParams(self):
        return self.coeffs

###############################################################################

class oimParamPolynomialWl(oimParamPolynomial):
    def _init(self, param, order=2,coeffs=None,x0=None, **kwargs):
        super()._init(param, dependence="wl", order=order,coeffs=coeffs,x0=x0)

###############################################################################

class oimParamPolynomialTime(oimParamPolynomial):
    def _init(self, param, order=2,coeffs=None,x0=None):
        super()._init(param, dependence="mjd", order=order,coeffs=coeffs,x0=x0)
        
###############################################################################        
class oimParamLinearRangeWl(oimParamInterpolator):

    def _init(self, param, wlmin=2e-6, wlmax=3e-6,values=[], kind="linear",**kwargs):

        self.kind=kind

        n = len(values)
        self.wlmin = (oimParam(**_standardParameters["wl"]))
        self.wlmin.name="wlmin"
        self.wlmin.description="Min of wl range"
        self.wlmin.value=wlmin
        self.wlmin.free=False

        self.wlmax = (oimParam(**_standardParameters["wl"]))
        self.wlmax.name="wlmax"
        self.wlmax.description="Max of the wl range"
        self.wlmax.value=wlmax
        self.wlmax.free=False


        self.values = []

        for i in range(n):
            self.values.append(oimParam(name=param.name, value=values[i],
                                        mini=param.min,maxi=param.max,
                                        description=param.description,
                                        unit=param.unit, free=param.free,
                                        error=param.error))

    def _interpFunction(self, wl, t):
    
          vals=np.array([vi.value for vi in self.values])
          nwl=vals.size
          wl0=np.linspace(self.wlmin.value,self.wlmax.value,num=nwl)
          print(wl0)
          return interp1d (wl0,vals,kind=self.kind,fill_value="extrapolate")(wl)
      
    def _getParams(self):
        params=[]
        params.extend(self.values)
        params.append(self.wlmin)
        params.append(self.wlmax)
        return params    
            

###############################################################################
# List of interpolators defined in oimodels
_interpolator={"wl":oimParamInterpolatorWl,
                "time":oimParamInterpolatorTime,
                "GaussWl":oimParamGaussianWl,
                "GaussTime":oimParamGaussianTime,
                "mGaussWl":oimParamMultipleGaussianWl,
                "mGaussTime":oimParamMultipleGaussianTime,
                "cosTime":oimParamCosineTime,
                "polyWl":oimParamPolynomialWl,
                "polyTime":oimParamPolynomialTime,
                "rangeWl":oimParamLinearRangeWl}


###############################################################################

class oimInterp(object):
    def __init__(self,name, **kwargs):
        self.kwargs=kwargs
        self.type=_interpolator[name]

###############################################################################
# Here is a list of standard parameters to be used when defining new components
_standardParameters = {
    "x": {"name": "x", "value": 0, "description": "x position", "unit": units.mas, "free": False},
    "y": {"name": "y", "value": 0, "description": "y position", "unit": units.mas, "free": False},
    "f": {"name": "f", "value": 1, "description": "flux", "unit": units.one,"mini":0,"maxi":1},
    "fwhm": {"name": "fwhm", "value": 0, "description": "FWHM", "unit": units.mas,"mini":0},
    "d": {"name": "d", "value": 0, "description": "Diameter", "unit": units.mas,"mini":0},
    "din": {"name": "din", "value": 0, "description": "Inner Diameter", "unit": units.mas,"mini":0},
    "dout": {"name": "dout", "value": 0, "description": "Outer Diameter", "unit": units.mas,"mini":0},
    "elong": {"name": "elong", "value": 1, "description": "Elongation Ratio", "unit": units.one,"mini":1},
    "pa": {"name": "pa", "value": 0, "description": "Major-axis Position angle", "unit": units.deg,"mini":-180,"maxi":180},
    "skw": {"name": "skw", "value": 0, "description": "Skewedness", "unit": units.one,"mini":0,"maxi":1},
    "skwPa": {"name": "skwPa", "value": 0, "description": "Skewedness Position angle", "unit": units.deg,"mini":-180,"maxi":180},
    "pixSize": {"name": "pixSize", "value": 0.1, "description": "Pixel Size", "unit": units.mas,"mini":0},
    "dim": {"name": "dim", "value": 128, "description": "Dimension in pixel", "unit": units.one, "free": False,"mini":1},
    "wl": {"name": "wl", "value": 0, "description": "Wavelength", "unit": units.m, "mini": 0},
    "mjd": {"name": "mjd", "value": 0, "description": "MJD", "unit": units.day}
}
