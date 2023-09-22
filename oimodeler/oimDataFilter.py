# -*- coding: utf-8 -*-
"""Data filtering/modifying"""
import numpy as np
from .oimUtils import cutWavelengthRange,shiftWavelength, spectralSmoothing,\
    binWavelength, oifitsFlagWithExpression, \
    computeDifferentialError, setMinimumError, \
    getDataArrname,getDataType, _oimDataType, _oimDataTypeArr

###############################################################################
class oimDataFilterComponent(object):
    """Base class for data filter"""
    name = "Generic Filter"
    shortname = "Genfilt"
    description = "This is the class from which all filters derived"

    def __init__(self, **kwargs):
        self.params = {}

        self.params["targets"] = "all"
        self.params["arr"] = "all"

        self._eval(**kwargs)

    def _eval(self, **kwargs):
        for key, value in kwargs.items():
            if key in self.params.keys():
                self.params[key] = value

    def _filteringFunction(self, data):
        pass

    def applyFilter(self, data):
        if type(self.params["targets"]) != type([]):
            self.params["targets"] = [self.params["targets"]]

        if type(self.params["arr"]) != type([]):
            self.params["arr"] = [self.params["arr"]]

        if self.params["targets"] == ["all"]:
            idx = list(range(len(data)))
        else:
            idx = self.params["targets"]

        for datai in [data[i] for i in idx]:
            self._filteringFunction(datai)

###############################################################################
class oimRemoveArrayFilter(oimDataFilterComponent):
    """Simple filter removing arrays by type"""
    name = "Remove array by type Filter"
    shortname = "RemArrFilt"
    description = "Remove array by type Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data):

        for arri in self.params["arr"]:
            while len(np.where(np.array([t.name for t in data]) == arri)[0]) != 0:
                data.pop(arri)

###############################################################################
class oimWavelengthRangeFilter(oimDataFilterComponent):
    """Filter for cutting wavelength range"""
    name = "Wavelength range Filter"
    shortname = "WlRgFilt"
    description = "Wavelength range Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["wlRange"] = []
        self.params["addCut"] = []
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        cutWavelengthRange(data, wlRange=self.params["wlRange"],
                           addCut=self.params["addCut"])

###############################################################################
class oimDataTypeFilter(oimDataFilterComponent):
    """ """
    name = "Filtering by datatype"
    shortname = "DTFilt"
    description = "Filtering by datatype : VIS2DATA, VISAMP..."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["dataType"] = []
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        if type(self.params["dataType"]) != type([]):
            self.params["dataType"] = [self.params["dataType"]]

        for dtype in self.params["dataType"]:
            idx = np.where(np.array(_oimDataType) == dtype)[0]
            if idx.size == 1:
                dtypearr = _oimDataTypeArr[idx[0]]

                for datai in data:
                    if datai.name == dtypearr:
                        datai.data[dtype] *= 0
                        
###############################################################################
class oimKeepDataType(oimDataFilterComponent):
    """ """
    name = "Keep datatype filter"
    shortname = "KeepDTFilt"
    description = "Keep atatype that are listed: VIS2DATA, VISAMP..."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["dataType"] = []
        self.params["removeArrIfPossible"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        if type(self.params["dataType"]) != type([]):
            self.params["dataType"] = [self.params["dataType"]]

        dataType=self.params["dataType"]
        #dataType=["VISAMP","VISPHI","T3PHI"]
        arr0=np.array(["PRIMARY","OI_ARRAY","OI_WAVELENGTH","OI_TARGET"])
        arr2Keep=np.unique(np.array([getDataArrname(dti) for dti in dataType]))

        hduname=[hdu.name for hdu in data]

        arr2remove=[]
        for ihdu,hdunamei in enumerate(hduname):
            if not(hdunamei in arr0 or hdunamei in arr2Keep):
                   extver=data[ihdu].header["EXTVER"]
                   arr2remove.append((hdunamei,extver))
            elif hdunamei in arr2Keep:
                dataTypesi=getDataType(hdunamei)
                for dataTypeij in dataTypesi:
                    if not(dataTypeij in dataType):
                        data[ihdu].data[dataTypeij][:]= 0
        for arr2removei in arr2remove:
            data.pop(arr2removei)    


###############################################################################
class oimWavelengthShiftFilter(oimDataFilterComponent):
    """Filter for shifting wavelength"""
    name = "Shift Wavelength Filter"
    shortname = "WlShFilt"
    description = "Wavelength Shift Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["wlShift"] = 0
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        shiftWavelength(data, self.params["wlShift"])


###############################################################################
class oimWavelengthSmoothingFilter(oimDataFilterComponent):
    """Filter for Smoothing wavelength"""
    name = "Wavelength Smoothing Filter"
    shortname = "WlSmFilt"
    description = "Wavelength Smoothing Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["smoothPix"] = 2
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        spectralSmoothing(data, self.params["smoothPix"], 
                        cols2Smooth=self.params["arr"],
                        normalizeError=self.params["normalizeError"])
        
###############################################################################
class oimWavelengthBinningFilter(oimDataFilterComponent):
    """Filter for binning wavelength"""
    name = "Wavelength binning Filter"
    shortname = "WlBinFilt"
    description = "Wavelength binning Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["bin"] = 3
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        binWavelength(data, self.params["bin"],
                        normalizeError=self.params["normalizeError"])


###############################################################################
class oimFlagWithExpressionFilter(oimDataFilterComponent):
    """Flaging based on expression """
    name = "Flag With Expression filter"
    shortname = "FlagExprFilt"
    description = "Flaging based on expression"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["expr"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        oifitsFlagWithExpression(data, self.params["arr"],1,self.params["expr"],
                                 keepOldFlag=self.params["keepOldFlag"])

###############################################################################
class oimDiffErrFilter(oimDataFilterComponent):
    """Compute differential error from std of signal inside or outside a range """
    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = "Compute differential error from std of signal inside or outside a range"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["ranges"] = [[0,5]]
        self.params["rangeType"] = "index"
        self.params["excludeRange"] = False
        self.params["dataType"] ="VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data):

        computeDifferentialError(data,ranges=self.params["ranges"],
                         excludeRange=self.params["excludeRange"],
                         rangeType=self.params["rangeType"],
                         dataType=self.params["dataType"],
                         extver=[None])

###############################################################################
class oimSetMinErrFilter(oimDataFilterComponent):
    """Set minimum error on data in % for vis ans deg for phases"""
    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = "Compute differential error from std of signal inside or outside a range"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["values"] = 5
        self.params["dataType"] ="VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data):

        setMinimumError(data,self.params["dataType"],
                         self.params["values"], extver=None)


###############################################################################
class oimDataFilter(object):
    """Class for data filter stack"""

    def __init__(self, filters=[]):
        if isinstance(filters,oimDataFilterComponent):
            filters=[filters]
        self.filters = filters

    def applyFilter(self, data):
        for filt in self.filters:
            filt.applyFilter(data)
