# -*- coding: utf-8 -*-
"""Data filtering/modifying"""

from __future__ import annotations

from fnmatch import fnmatch

import numpy as np

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
    """Base class for data filter.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Can be `"all"` or a
        list of indices corresponding to the list of input data. Defaults to `"all"`.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Generic Filter"
    shortname = "Genfilt"
    description = "This is the class from which all filters derived"

    def __init__(self, **kwargs) -> None:
        self.params = {}

        self.params["targets"] = "all"
        self.params["arr"] = "all"

        self._eval(**kwargs)

    def _eval(self, **kwargs) -> None:
        """Evaluates the ``kwargs`` and passes them to the ``self.params``
        dictionary if they exist in it."""
        for key, value in kwargs.items():
            if key in self.params:
                self.params[key] = value

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""

    def applyFilter(self, data) -> None:
        """Applies the filter to the data."""
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

    def __str__(self) -> str:
        txt = self.name
        for key, value in self.params.items():
            txt += "\n"
            txt += f"{key}:".ljust(10)
            txt += f"{value}"

        return txt


class oimDataFilter:
    """Class for data filter stack."""

    def __init__(
        self,
        filters: oimDataFilterComponent | list[oimDataFilterComponent] = [],
    ) -> None:
        if isinstance(filters, oimDataFilterComponent):
            filters = [filters]

        self.filters = filters

    def applyFilter(self, data) -> None:
        """Applies all filters in ``self.filters`` to the data."""
        for filt in self.filters:
            filt.applyFilter(data)


class oimRemoveArrayFilter(oimDataFilterComponent):
    """Filter that removes array(s)/table(s) from OIFITS file by name.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Remove Array Filter"
    shortname = "RemArrFilt"
    description = "Removing array(s)/table(s) by name: OI_VIS2, OI_T3..."

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        for arri in self.params["arr"]:
            while (
                len(np.where(np.array([t.name for t in data]) == arri)[0]) != 0
            ):
                data.pop(arri)


class oimRemoveInsnameFilter(oimDataFilterComponent):
    """Filter that removes arrays/tables from OIFITS file by ``"INSNAME"``.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    insname : str or list of str
        One or more `"INSNAME"` that is/are to be removed. Defaults to `None`.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Remove Insname Filter"
    shortname = "RemInsFilt"
    description = "Remove arrays/tables by insname(s)"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["insname"] = None
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        to_remove = []
        insnameToRemove = self.params["insname"]
        if type(insnameToRemove) != type([]):
            insnameToRemove = [insnameToRemove]

        for di in data:
            if "INSNAME" in di.header:
                for insnamei in insnameToRemove:
                    if fnmatch(di.header["INSNAME"], insnamei):
                        to_remove.append(di)

        for di in to_remove:
            data.pop(di)


class oimDataTypeFilter(oimDataFilterComponent):
    """Filter that sets column values to ``0`` by column name(s).

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    dataType : str or list of str
        The OIFITS datatypes(s)/column(s) to be removed.  Defaults to `[]`.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Data Type Filter"
    shortname = "DTFilt"
    description = (
        "Sets column names to 0 by colum name(s) : VIS2DATA, VISAMP..."
    )

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["dataType"] = []
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        if type(self.params["dataType"]) != type([]):
            self.params["dataType"] = [self.params["dataType"]]

        for dtype in self.params["dataType"]:
            idx = np.where(np.array(_oimDataType) == dtype)[0]
            if idx.size == 1:
                dtypearr = _oimDataTypeArr[idx[0]]

                for datai in data:
                    if datai.name == dtypearr:
                        datai.data[dtype] *= 0


class oimKeepDataTypeFilter(oimDataFilterComponent):
    """Filter that removes all columns except those specified by name.

    If no data columns remain in array(s)/table(s) after filtering and
    ``removeArrIfPossible=True` then they are removed as well.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    dataType : str or list of str
        The OIFITS datatype(s)/column(s) to be kept. Defaults to `[]`.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Keep Datatype filter"
    shortname = "KeepDTFilt"
    description = "Keep columns that are specified: VIS2DATA, VISAMP..."

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["dataType"] = []
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        if type(self.params["dataType"]) != type([]):
            self.params["dataType"] = [self.params["dataType"]]

        dataType = self.params["dataType"]
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


class oimFlagWithExpressionFilter(oimDataFilterComponent):
    """Filter that flags based on an expression.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    expr : str
        The expression to be applied. Defaults to `""`.
    keepOldFlag : bool, optional
        If `True`, the old flag(s) are kept. Defaults to `True`.

    See Also
    --------
    oimUtils.oifitsFlagWithExpression : Flag data with an expression.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Flag With Expression filter"
    shortname = "FlagWExprFilt"
    description = "Flags based on boolean expressions"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["expr"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        oifitsFlagWithExpression(
            data,
            self.params["arr"],
            None,
            self.params["expr"],
            keepOldFlag=self.params["keepOldFlag"],
        )


class oimWavelengthRangeFilter(oimDataFilterComponent):
    """Filter that cuts the wavelength range(s) of all data.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    wlRange : list of float
        The wavelength range after filtering. Defaults to `[]`.
    addCut : list of float, optional
    method : str, optional
        The method for . If `method="cut"`, the `.oimUtils.cutWavelengthRange` function
        is used, otherwise the `.oimUtils.oimFlagWithExpression`. Defaults to `"cut"`.

    See Also
    --------
    oimUtils.cutWavelengthRange : Cut the wavelength range of an OIFITS file.
    oimUtils.oifitsFlagWithExpression : Flag data with an expression.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Wavelength Range Filter"
    shortname = "WlRgFilt"
    description = "Cuts the wavelength range(s)"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["wlRange"] = []
        self.params["addCut"] = []
        self.params["method"] = "cut"
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
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


class oimWavelengthShiftFilter(oimDataFilterComponent):
    """Filter that shifts the wavelength table.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    wlShift : float
        The amount the wavelenght is shifted by. Defaults to `0`.

    See Also
    --------
    oimUtils.shiftWavelength : Shift the wavelength of an OIFITS file.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Wavelength Shift Filter"
    shortname = "WlShFilt"
    description = "Shifts the wavelength table"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["wlShift"] = 0
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        shiftWavelength(data, self.params["wlShift"])


class oimWavelengthSmoothingFilter(oimDataFilterComponent):
    """Filter for Smoothing wavelength.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    smoothPix : int, optional
        The kernel size of the smoothing. Defaults to `2`.
    normalizeError : bool, optional
        If `True`, the errors are normalised by the kernel size. Defaults to `True`.

    See Also
    --------
    oimUtils.spectralSmoothing : Smooth the spectral data of an OIFITS file.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Wavelength Smoothing Filter"
    shortname = "WlSmFilt"
    description = "Spectral smoothing "

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["smoothPix"] = 2
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        spectralSmoothing(
            data,
            self.params["smoothPix"],
            cols2Smooth=self.params["arr"],
            normalizeError=self.params["normalizeError"],
        )


class oimWavelengthBinningFilter(oimDataFilterComponent):
    """Filter for binning wavelength.

    Other Parameters
    ----------------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    bin : int
        The bin size. Defaults to `None`.
    normalizeError : bool, optional
        If `True`, the errors are normalised by the bin size. Defaults to `True`.

    See Also
    --------
    oimUtils.binWavelength : Bin the wavelength of an OIFITS file.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Wavelength Binning Filter"
    shortname = "WlBinFilt"
    description = "Spectral Binning"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["bin"] = None
        self.params["normalizeError"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        binWavelength(
            data,
            binSize=self.params["bin"],
            normalizeError=self.params["normalizeError"],
        )


class oimWavelengthIntpBinFilter(oimDataFilterComponent):
    """Filter that bins the wavelength to a specified grid. It also interpolates
    at the edges of the bins, ensuring a minimum number of elements.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    binGrid : array_like
        The grid that is to be achieved/binned to.
    binWindow : array_like, optional
        The bin windows that correspond to the `binGrid`. If `None`, the
        bin windows are computed from the distance between two elements in the
        `binGrid`. Defaults to `None`.
    resetFlags : bool, optional
        If `True`, resets all flags to `False` after binning. Defaults to `True`.
    averageError : bool, optional
        If `True`, forgoes the error propagation and simply averages the errors
        for each bin. Defaults to `False`.
    nSpecChannels : float, optional
        The number of spectral channels determined by the spectral resolution.
        Will be used to calculate the divisor within the error propagation.
        Defaults to `1.0`.

        .. math:: divisor = bin_elements / spectralChannels

    See Also
    --------
    oimUtils.intpBinWavelength : Bin the wavelength of an OIFITS file to a specified binGrid.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Wavelength Interpolation Binning Filter"
    shortname = "WlIntpBinFilt"
    description = (
        "Binning to wavelength grid with interpolation at window edges."
    )

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["binGrid"] = None
        self.params["binWindow"] = None
        self.params["resetFlags"] = True
        self.params["averageError"] = False
        self.params["nSpecChannels"] = 1.0

        # TODO: Remove this eventually. Just in place due to breaking change
        # after v0.9.8 and before next version release.
        if "spectralChannels" in kwargs:
            raise ValueError(
                "The kwarg 'spectralChannels' was renamed to 'nSpecChannels'."
            )

        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        intpBinWavelength(
            data,
            self.params["binGrid"],
            binWindow=self.params["binWindow"],
            resetFlags=self.params["resetFlags"],
            averageError=self.params["averageError"],
            nSpecChannels=self.params["nSpecChannels"],
        )


class oimKeepBaselinesFilter(oimDataFilterComponent):
    """Filter that removes all baselines except those specified by name.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    baselines : str or list of str
        The baseline(s) to be kept. Defaults to `""`.
    keepOldFlag : bool, optional
        If `True`, the old flag(s) are kept. Defaults to `True`.

    See Also
    --------
    oimUtils.oifitsKeepBaselines : Remove all baselines except those specified by name.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Keep Baseline Filter"
    shortname = "KeepBaseFilt"
    description = "Keeps baseline(s) specified by name"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["baselines"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        for arri in self.params["arr"]:
            oifitsKeepBaselines(
                data,
                arri,
                self.params["baselines"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimRemoveBaselinesFilter(oimDataFilterComponent):
    """Filter that removes all baselines specified by name.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    baselines : str or list of str
        The baseline(s) to be removed. Defaults to `""`.
    keepOldFlag : bool, optional
        If `True`, the old flag(s) are kept. Defaults to `True`.

    See Also
    --------
    oimUtils.oifitsRemoveBaselines : Remove all baselines specified by name.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Remove Baselines filter"
    shortname = "RemBaseFilt"
    description = "Removes baseline(s) specified by name"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["baselines"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        for arri in self.params["arr"]:
            oifitsRemoveBaselines(
                data,
                arri,
                self.params["baselines"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimKeepTelescopesFilter(oimDataFilterComponent):
    """Filter that removes all telescopes except those specified by name.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    telescopes : str or list of str
        The telescopes(s) to be kept. Defaults to `""`.
    keepOldFlag : bool, optional
        If `True`, the old flag(s) are kept. Defaults to `True`.

    See Also
    --------
    oimUtils.oifitsKeepTelescopes : Remove all telescopes except those specified by name.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Keep Telescopes Filter"
    shortname = "KeepTelFilt"
    description = "Keeps telescope(s) specified by name"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["telescopes"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        for arri in self.params["arr"]:
            oifitsKeepTelescopes(
                data,
                arri,
                self.params["telescopes"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimRemoveTelescopesFilter(oimDataFilterComponent):
    """Filter that removes all telescopes specified by name.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    telescopes : str or list of str
        The telescopes(s) to be removed. Defaults to `""`.
    keepOldFlag : bool, optional
        If `True`, the old flag(s) are kept. Defaults to `True`.

    See Also
    --------
    oimUtils.oifitsRemoveTelescopes : Remove all telescopes specified by name.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Remove Telescopes Filter"
    shortname = "RemTelFilt"
    description = "Removes telescope(s) specified by name"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["telescopes"] = ""
        self.params["keepOldFlag"] = True
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        for arri in self.params["arr"]:
            oifitsRemoveTelescopes(
                data,
                arri,
                self.params["telescopes"],
                keepOldFlag=self.params["keepOldFlag"],
            )


class oimResetFlagsFilter(oimDataFilterComponent):
    """Filter that unflags all data (i.e. sets flags to `False`).

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.

    See Also
    --------
    oimUtils.oifitsFlagWithExpression : Flag data with an expression.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Reset Flags Filter"
    shortname = "ResFlagsFilt"
    description = "Sets all flags to `False`"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        oifitsFlagWithExpression(
            data, self.params["arr"], None, "False", keepOldFlag=False
        )


# TODO: Finish this function's documentation.
class oimDiffErrFilter(oimDataFilterComponent):
    """Compute differential error from std of signal inside or outside a range.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    ranges : list of list of float, optional
        Defaults to `[[0, 5]]`.
    rangeType : str, optional
        Defaults to `"index"`.
    excludeRange : bool, optional
        Defaults to `False`.
    dataType : str or list of str, optional
        The OIFITS datatype(s)/column(s) to be kept. Defaults to `"VISPHI"`.

    See Also
    --------
    oimUtils.computeDifferentialError : Compute the differential error.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = (
        "Compute differential error from std inside or outside a range"
    )

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["ranges"] = [[0, 5]]
        self.params["rangeType"] = "index"
        self.params["excludeRange"] = False
        self.params["dataType"] = "VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        computeDifferentialError(
            data,
            ranges=self.params["ranges"],
            excludeRange=self.params["excludeRange"],
            rangeType=self.params["rangeType"],
            dataType=self.params["dataType"],
            extver=[None],
        )


class oimSetMinErrFilter(oimDataFilterComponent):
    """Set minimum error on data in % for vis and deg for phases.

    Other Parameters
    ----------
    targets : str or list of int, optional
        The targets that this filter is applied to. Either ``"all"`` or a
        list of indices corresponding to the list of input data. Defaults to ``"all"``.
    arr : str or list of str, optional
        The OIFITS array(s)/table(s) this filter is applied to. Can be `"all"` or a
        string or list of strings with (a) table name(s). Defaults to `"all"`.
    values : float or list of float
        The minimum error values corresponding to the column(s)/datatype(s). Defaults to `5`.
    dataType : str or list of str, optional
        The OIFITS datatype(s)/column(s) to be kept. Defaults to `"VISPHI"`.

    See Also
    --------
    oimUtils.setMinimumError : Set the minimum error of a given data type to a given value.

    Notes
    -----
    All keyword arguments are passed to the `self.params` dictionary, which is then used
    to pass it to the underlying filter class/function.
    """

    name = "Differential Error Filter"
    shortname = "DiffErrFilt"
    description = "Set minimum error on data in % for vis and deg for phases"

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.params["values"] = 5
        self.params["dataType"] = "VISPHI"
        self._eval(**kwargs)

    def _filteringFunction(self, data) -> None:
        """The filter applied to the data."""
        setMinimumError(
            data, self.params["dataType"], self.params["values"], extver=None
        )
