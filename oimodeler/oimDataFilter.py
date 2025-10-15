# -*- coding: utf-8 -*-
"""Data filtering/modifying"""
import numpy as np
from fnmatch import fnmatch

from .oimUtils import (
    _oimDataType,
    _oimDataTypeArr,
    binWavelength,
    computeDifferentialError,
    cutWavelengthRange,
    getDataArrname,
    getDataType,
    intpBinWavelength,
    oifitsFlagWithExpression,
    oifitsKeepBaselines,
    oifitsKeepTelescopes,
    oifitsRemoveBaselines,
    oifitsRemoveTelescopes,
    setMinimumError,
    shiftWavelength,
    spectralSmoothing,
)


class oimDataFilterComponent:
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

    def __str__(self):
        txt=self.name 
        for key,value in self.params.items():
            txt+="\n"
            txt+= f"{key}:".ljust(10)
            txt+= f"{value}"
        return txt
        

class oimDataFilter:
    """Class for data filter stack"""

    def __init__(self, filters=[]):
        if isinstance(filters, oimDataFilterComponent):
            filters = [filters]
        self.filters = filters

    def applyFilter(self, data):
        for filt in self.filters:
            filt.applyFilter(data)


class oimRemoveArrayFilter(oimDataFilterComponent):
    """Simple filter removing arrays by type"""

    name = "Remove array by type Filter"
    shortname = "RemArrFilt"
    description = "Removing arrays by type: OI_VIS2, OI_T3..."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        for arri in self.params["arr"]:
            while (
                len(np.where(np.array([t.name for t in data]) == arri)[0]) != 0
            ):
                data.pop(arri)


class oimRemoveInsnameFilter(oimDataFilterComponent):
    """Simple filter removing arrays by type"""

    name = "Remove arrays by insname Filter"
    shortname = "RemInsFilt"
    description = "Remove array by insname"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["insname"] = None
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        to_remove = []
        insnameToRemove = self.params["insname"]
        if type(insnameToRemove)!=type([]):
            insnameToRemove=[insnameToRemove]
        #print(insnameToRemove)
        for di in data:
            if "INSNAME" in di.header:
                for insnamei in insnameToRemove:   
                    if fnmatch(di.header["INSNAME"],insnamei):
                        to_remove.append(di)

        
        for di in to_remove:
            data.pop(di)


class oimWavelengthRangeFilter(oimDataFilterComponent):
    """Filter for cutting wavelength range"""

    name = "Wavelength range Filter"
    shortname = "WlRgFilt"
    description = "Cutting wavelength range(s)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["wlRange"] = []
        self.params["addCut"] = []
        self.params["method"] = "cut"
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        if self.params["method"] == "cut":
            cutWavelengthRange(
                data,
                wlRange=self.params["wlRange"],
                addCut=self.params["addCut"],
            )
        else:
            expr = ""
            wlRange = np.array(self.params["wlRange"])
            if wlRange.ndim == 1:
                wlRange = wlRange.reshape((1, len(wlRange)))

            for wlRangei in wlRange:
                expr += (
                    f"((EFF_WAVE<{wlRangei[0]}) | (EFF_WAVE>{wlRangei[1]})) &"
                )
            expr = expr[:-1]
            oifitsFlagWithExpression(
                data, self.params["arr"], None, expr, keepOldFlag=True
            )


class oimDataTypeFilter(oimDataFilterComponent):
    """ """

    name = "Filtering by datatype"
    shortname = "DTFilt"
    description = "Filtering by datatype(s) : VIS2DATA, VISAMP..."

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

        dataType = self.params["dataType"]
        # dataType=["VISAMP","VISPHI","T3PHI"]
        arr0 = np.array(["PRIMARY", "OI_ARRAY", "OI_WAVELENGTH", "OI_TARGET"])
        arr2Keep = np.unique(
            np.array([getDataArrname(dti) for dti in dataType])
        )

        hduname = [hdu.name for hdu in data]

        arr2remove = []
        for ihdu, hdunamei in enumerate(hduname):
            if not (hdunamei in arr0 or hdunamei in arr2Keep):
                extver = data[ihdu].header.get("EXTVER", 1)
                arr2remove.append((hdunamei, extver))
            elif hdunamei in arr2Keep:
                dataTypesi = getDataType(hdunamei)
                for dataTypeij in dataTypesi:
                    if dataTypeij not in dataType:
                        data[ihdu].data[dataTypeij][:] = 0
        for arr2removei in arr2remove:
            data.pop(arr2removei)


class oimWavelengthShiftFilter(oimDataFilterComponent):
    """Filter for shifting wavelength"""

    name = "Shift Wavelength Filter"
    shortname = "WlShFilt"
    description = "Shifting wavelength table"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["wlShift"] = 0
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        shiftWavelength(data, self.params["wlShift"])


class oimWavelengthSmoothingFilter(oimDataFilterComponent):
    """Filter for Smoothing wavelength"""

    name = "Wavelength Smoothing Filter"
    shortname = "WlSmFilt"
    description = "Spectral smoothing "

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["smoothPix"] = 2
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        spectralSmoothing(
            data,
            self.params["smoothPix"],
            cols2Smooth=self.params["arr"],
            normalizeError=self.params["normalizeError"],
        )


class oimWavelengthBinningFilter(oimDataFilterComponent):
    """Filter for binning wavelength"""

    name = "Wavelength binning Filter"
    shortname = "WlBinFilt"
    description = "Spectral Binning"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["bin"] = None
        self.params["binGrid"] = None
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        binWavelength(
            data,
            binsize=self.params["bin"],
            normalizeError=self.params["normalizeError"],
        )


class oimWavelengthIntpBinFilter(oimDataFilterComponent):
    """Filter that bins the wavelength to a specified grid.
    It also interpolates at the edges of the bins, to ensure a minimum
    number of elements.

    Parameters
    ----------
    binGrid : array_like
        The grid to bin to.
    binWindow : array_like, optional
        The bin windows that correspond to the binGrid elements.
        If None, computes the bin windows from the distance between two
        elements in the binGrid. Default is None.
    resetFlags : bool, optional
        If True, resets the flags after binning. Default is True.
    averageError : bool, optional
        If True, forgoes error propagation and simply averages the errors
        for each bin. Default is False.
    spectralChannels : int, optional
        The number of channels of the set bin resolution. Will be used to
        calculate the divisor within the error propagation. Default is 1.0.

        .. math:: divisor = bin_elements / spectralChannels


    """

    name = "Wavelength Interpolation Binning Filter"
    shortname = "WlIntpBinFilt"
    description = (
        "Binning to wavelength grid with interpolation at window edges."
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["binGrid"] = None
        self.params["binWindow"] = None
        self.params["resetFlags"] = True
        self.params["averageError"] = False
        self.params["spectralChannels"] = 1
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        intpBinWavelength(
            data,
            self.params["binGrid"],
            binWindow=self.params["binWindow"],
            resetFlags=self.params["resetFlags"],
            averageError=self.params["averageError"],
            spectralChannels=self.params["spectralChannels"],
        )


class oimFlagWithExpressionFilter(oimDataFilterComponent):
    """Flaging based on expression"""

    name = "Flag With Expression filter"
    shortname = "FlagExprFilt"
    description = "Flagging based on boolean expressions"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["expr"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        oifitsFlagWithExpression(
            data,
            self.params["arr"],
            None,
            self.params["expr"],
            keepOldFlag=self.params["keepOldFlag"],
        )


class oimKeepBaselinesFilter(oimDataFilterComponent):
    """Select baselines to keep"""

    name = "Baseline selection filter"
    shortname = "KeepBaselinesFilt"
    description = "Selection based on baseline name(s)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["baselines"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        for arri in self.params["arr"]:
            oifitsKeepBaselines(
                data,
                arri,
                self.params["baselines"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimRemoveBaselinesFilter(oimDataFilterComponent):
    """Select baselines to remove"""

    name = "Baseline selection filter"
    shortname = "RemoveBAselineFilt"
    description = "Selection based on baseline name(s)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["baselines"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        for arri in self.params["arr"]:
            oifitsRemoveBaselines(
                data,
                arri,
                self.params["baselines"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimKeepTelescopesFilter(oimDataFilterComponent):
    """Select telescopes to keep"""

    name = "Telescopes selection filter"
    shortname = "KeepTelescopesFilt"
    description = "Selection based on telescope name(s)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["telescopes"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        for arri in self.params["arr"]:
            oifitsKeepTelescopes(
                data,
                arri,
                self.params["telescopes"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimRemoveTelescopesFilter(oimDataFilterComponent):
    """Select telescopes to remove"""

    name = "Telescope selection filter"
    shortname = "RemoveTelescopeFilt"
    description = "Selection based on telescope name(s)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["telescopes"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        for arri in self.params["arr"]:
            oifitsRemoveTelescopes(
                data,
                arri,
                self.params["telescopes"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimResetFlags(oimDataFilterComponent):
    """Unflag all data"""

    name = "Reset flags filter"
    shortname = "ResFlagsFilt"
    description = "Set all the flags to False"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        oifitsFlagWithExpression(
            data, self.params["arr"], None, "False", keepOldFlag=False
        )


class oimDiffErrFilter(oimDataFilterComponent):
    """Compute differential error from std of signal inside or outside a range"""

    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = (
        "Compute differential error from std inside or outside a range"
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["ranges"] = [[0, 5]]
        self.params["rangeType"] = "index"
        self.params["excludeRange"] = False
        self.params["dataType"] = "VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        computeDifferentialError(
            data,
            ranges=self.params["ranges"],
            excludeRange=self.params["excludeRange"],
            rangeType=self.params["rangeType"],
            dataType=self.params["dataType"],
            extver=[None],
        )


class oimSetMinErrFilter(oimDataFilterComponent):
    """Set minimum error on data in % for vis ans deg for phases"""

    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = "Set minimum error on data in % for vis ans deg for phases"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["values"] = 5
        self.params["dataType"] = "VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data):
        setMinimumError(
            data, self.params["dataType"], self.params["values"], extver=None
        )
