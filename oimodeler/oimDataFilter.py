# -*- coding: utf-8 -*-
"""Data filtering/modifying"""
import numpy as np

from .oimData import _oimDataType, _oimDataTypeArr
from .oimUtils import cutWavelengthRange


class oimDataFilterComponent(object):
    """Base class for data filter"""
    name = "Generic Filter"
    shortname = "Gen filter"
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


class oimRemoveArrayFilter(oimDataFilterComponent):
    """Simple filter removing arrays by type"""
    name = "Remove array by type Filter"
    shortname = "Remove Arr"
    description = "Remove array by type Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data):

        for arri in self.params["arr"]:
            while len(np.where(np.array([t.name for t in data]) == arri)[0]) != 0:
                data.pop(arri)


class oimWavelengthRangeFilter(oimDataFilterComponent):
    """Filter for cutting wavelength range"""
    name = "Wavelength range Filter"
    shortname = "WlRange Filter"
    description = "Wavelength range Filter"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["wlRange"] = []
        self.params["addCut"] = []
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        cutWavelengthRange(data, wlRange=self.params["wlRange"],
                           addCut=self.params["addCut"])


class oimDataTypeFilter(oimDataFilterComponent):
    """ """
    name = "Filtering by datatype"
    shortname = "DataType Filter"
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


class oimDataFilter(object):
    """Class for data filter stack"""

    def __init__(self, filters=[]):
        self.filters = filters

    def applyFilter(self, data):
        for filt in self.filters:
            filt.applyFilter(data)
