# -*- coding: utf-8 -*-
"""Data for optical interferometry"""
import os
from enum import IntFlag
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np
from astropy.io import fits

from .oimUtils import hdulistDeepCopy, _oimDataTypeArr
from .oimDataFilter import oimDataFilter, oimDataFilterComponent


def oimDataGetWl(data: fits.HDUList, array: fits.BinTableHDU,
                 dwl: Optional[bool] = True) -> np.ndarray:
    """Gets the wavelengths from the data.

    Parameters
    ----------
    data : astropy.io.fits.HDUList
        The data from which to get the wavelengths.
    array : astropy.io.fits.BinTableHDU
    dwl : bool, optional

    Returns
    -------
    wavelengths : numpy.ndarray
    """
    instrument = array.header["INSNAME"]
    oi_wavelengths = [arri for arri in data
                      if (arri.name == "OI_WAVELENGTH"
                          and arri.header["INSNAME"] == instrument)][0]
    if not dwl:
        return oi_wavelengths.data["EFF_WAVE"]
    return oi_wavelengths.data["EFF_WAVE"], oi_wavelengths.data["EFF_BAND"]


class oimDataType(IntFlag):
    """Data types for the oifits data."""
    NONE = 0
    VIS2DATA = 1
    VISAMP_ABS = 2
    VISAMP_DIF = 4
    VISAMP_COR = 8
    VISPHI_ABS = 16
    VISPHI_DIF = 32
    T3AMP = 64
    T3PHI = 128
    FLUXDATA = 256


def oimGetDataValErrAndTypeFlag(table: fits.BinTableHDU) -> Tuple[np.ndarray]:
    """Get the data, error and flag arrays from the oifits data.

    Parameters
    ----------
    table : astropy.io.fits.BinTableHDU
        The table from which to get the data, error and flag arrays.

    Returns
    -------
    dtype : oimDataType
        The type of the data.
    val : numpy.ndarray
        The data array.
    err : numpy.ndarray
        The error array.
    flag : numpy.ndarray
        The flag array.
    """
    dtype, val, err, flag = oimDataType(0), [], [], []
    if table.name == "OI_VIS2":
        nv2 = np.size(np.where(table.data["VIS2DATA"] != 0))
        if nv2 != 0:
            val.append(table.data["VIS2DATA"])
            err.append(table.data["VIS2ERR"])
            flag.append(table.data["FLAG"])
            dtype |= oimDataType.VIS2DATA

    if table.name == "OI_VIS":
        nvamp = np.size(np.where(table.data["VISAMP"] != 0))
        if nvamp != 0:
            val.append(table.data["VISAMP"])
            err.append(table.data["VISAMPERR"])
            flag.append(table.data["FLAG"])
            try:
                if table.header["AMPTYP"].lower() == "absolute":
                    dtype |= oimDataType.VISAMP_ABS
                elif table.header["AMPTYP"].lower() == "differential":
                    dtype |= oimDataType.VISAMP_DIF
                else:
                    dtype |= oimDataType.VISAMP_COR
            except:
                dtype |= oimDataType.VISAMP_ABS

        nvphi = np.size(np.where(table.data["VISPHI"] != 0))
        if nvphi != 0:
            val.append(table.data["VISPHI"])
            err.append(table.data["VISPHIERR"])
            flag.append(table.data["FLAG"])
            try:
                if table.header["PHITYP"].lower() == "absolute":
                    dtype |= oimDataType.VISPHI_ABS
                else:
                    dtype |= oimDataType.VISPHI_DIF
            except:
                dtype |= oimDataType.VISPHI_DIF

    if table.name == "OI_T3":
        t3amp = np.size(np.where(table.data["T3AMP"] != 0))
        if t3amp != 0:
            val.append(table.data["T3AMP"])
            err.append(table.data["T3AMPERR"])
            flag.append(table.data["FLAG"])
            dtype |= oimDataType.T3AMP
        t3phi = np.size(np.where(table.data["T3PHI"] != 0))
        if t3phi != 0:
            val.append(table.data["T3PHI"])
            err.append(table.data["T3PHIERR"])
            flag.append(table.data["FLAG"])
            dtype |= oimDataType.T3PHI

    if table.name == "OI_FLUX":
        try:
            nflx = np.size(np.where(table.data["FLUXDATA"] != 0))
            if nflx != 0:
                val.append(table.data["FLUXDATA"])
                err.append(table.data["FLUXERR"])
                flag.append(table.data["FLAG"])
                dtype |= oimDataType.FLUXDATA
        except:
            nflx = np.size(np.where(table.data["FLUX"] != 0))
            if nflx != 0:
                val.append(table.data["FLUX"])
                err.append(table.data["FLUXERR"])
                flag.append(table.data["FLAG"])
                dtype |= oimDataType.FLUXDATA
    return dtype, val, err, flag


def oimDataCheckData(table: fits.BinTableHDU) -> List[str]:
    cdata = []
    if np.size(np.where(table.data["FLAG"] == False)) != 0:
        if table.name == "OI_VIS2":
            nv2 = np.size(np.where(table.data["VIS2DATA"] != 0))
            if nv2 != 0:
                cdata.append("VIS2DATA")

        if table.name == "OI_VIS":
            nvamp = np.size(np.where(table.data["VISAMP"] != 0))
            if nvamp != 0:
                cdata.append("VISAMP")
            nvphi = np.size(np.where(table.data["VISPHI"] != 0))
            if nvphi != 0:
                cdata.append("VISPHI")

        if table.name == "OI_T3":
            t3amp = np.size(np.where(table.data["T3AMP"] != 0))
            if t3amp != 0:
                cdata.append("T3AMP")
            t3phi = np.size(np.where(table.data["T3PHI"] != 0))
            if t3phi != 0:
                cdata.append("T3PHI")

        if table.name == "OI_FLUX":
            try:
                nflx = np.size(np.where(table.data["FLUXDATA"] != 0))
                if nflx != 0:
                    cdata.append("FLUXDATA")
            except:
                nflx = np.size(np.where(table.data["FLUX"] != 0))
                if nflx != 0:
                    cdata.append("FLUX")
    return cdata


def oimDataGetVectCoord(data: fits.HDUList, arr: fits.BinTableHDU) -> Tuple[np.ndarray]:
    """Get the (u, v)-coordinates from the data.

    Parameters
    ----------
    data : astropy.io.fits.HDUList
        The data from which to get the (u, v)-coordinates.

    Returns
    -------
    u : numpy.ndarray
        The u-coordinates.
    v : numpy.ndarray
        The v-coordinates.
    wl : numpy.ndarray
        The wavelengths.
    dwl : numpy.ndarray
        The wavelength bandwidths.
    mjd : numpy.ndarray
        The modified Julian dates.
    nB : int
        The number of baselines.
    nwl : int
        The number of wavelengths.
    """
    wl, dwl = oimDataGetWl(data, arr)
    nwl = np.size(wl)

    # TODO: Here I assume that the mjd is the same for all baselines
    # This is the case for all know instruments, but it is not a requirement of
    # the oifits format. If the mjd are different, then we should have more data
    # passed to the simulator, i.e we need to compute the zero frequency
    # not only for all wl but also all mjd.

    mjd0 = arr.data["MJD"][0]
    mjd = np.outer(mjd0, np.ones(nwl)).flatten()

    # NOTE: Zero freq vector for vis normalization
    uv0 = wl*0

    if arr.name == "OI_T3":
        nB = np.shape(arr.data["T3PHI"])[0]*3+1
        um1, um2 = map(lambda x: arr.data[f"U{x}COORD"], (1, 2))
        vm1, vm2 = map(lambda x: arr.data[f"V{x}COORD"], (1, 2))
        um3, vm3 = um1+um2, vm1+vm2

        u1, u2, u3 = map(lambda x: np.outer(x, 1./wl).flatten(), (um1, um2, um3))
        v1, v2, v3 = map(lambda x: np.outer(x, 1./wl).flatten(), (vm1, vm2, vm3))
        u = np.concatenate((u1, u2, u3))
        v = np.concatenate((v1, v2, v3))

    elif arr.name == "OI_FLUX":
        try:
            nB = np.shape(arr.data["FLUXDATA"])[0]
        except:
            nB = np.shape(arr.data["FLUX"])[0]
        u = np.zeros(nB*nwl)
        v = np.zeros(nB*nwl)

    else:
        um, vm = arr.data["UCOORD"], arr.data["VCOORD"]
        nB = np.size(um)+1
        u, v = map(lambda x: np.outer(x, 1./wl).flatten(), (um, vm))

    if arr.name != "OI_FLUX":
        u = np.concatenate((uv0, u))
        v = np.concatenate((uv0, v))

    mjd = np.outer(np.ones(nB), mjd).flatten()
    wl = np.outer(np.ones(nB), wl).flatten()
    dwl = np.outer(np.ones(nB), dwl).flatten()
    return u, v, wl, dwl, mjd, nB, nwl


class oimData:
    def __init__(self, data: Any) -> None:
        self.data = data


def loadOifitsData(input: Union[str, Path, List[str],
                                List[Path], fits.HDUList, oimData],
                   mode: Optional[str] = "listOfHdlulist") -> oimData:
    """Return the oifits data from either filenames, already opened oifts or a
    oimData boject as either a list of hdlulist (default) or as a oimData
    object using the option mode="oimData".

    Parameters
    ----------
    input : string or pathlib.Path or list of str or list of pathlib.Path
            or astropy.io.fits.hdu.hdulist.HDUList or oimodeler.oimData
        The data to deal with. Can be a oimData object, a hdulist, a string
        representing a filename or a list of these kind of object
    mode : str, optional
        The type of the return data, either "listOfHdlulist" or "oimData"
        The default is "listOfHdlulist"

    Returns
    -------
    data : .oimData
    """
    if isinstance(input, oimData):
        if mode == "oimData":
            data = input
        else:
            data = input.data
    else:
        if isinstance(input, (fits.hdu.hdulist.HDUList, str, Path)):
            input = [input]

        if isinstance(input, list):
            data = []

            for elem in input:
                if isinstance(elem, fits.hdu.hdulist.HDUList):
                    data.append(elem)
                else:
                    try:
                        data.append(fits.open(elem))
                    except:
                        raise ValueError("The path does not exist or is not a"\
                                         " valid fits files")
        else:
            raise TypeError("Only oimData, hdulist, Path or string, or list of"\
                            " these kind of objects allowed ")

        if mode == "oimData":
             data = oimData(data)
    return data


class oimData(object):
    """A class to hold and manipulate data

    Parameters
    ----------
    dataOrFilename : any
        The data to add. Can be a filename, a list of filenames, a hdulist
        or a list of hdulist.
    filt : oimDataFilter, optional
        The filter to use. The default is None.
    """

    def __init__(self, dataOrFilename: Optional[Any] = None,
                 filt: Optional[oimDataFilter] = None) -> None:
        """Initialize the class with the data and the filter to use."""
        self._data = []
        self.dataInfo = []
        self.vect_u = None
        self.vect_v = None
        self.vect_wl = None
        self.vect_dwl = None
        self.vect_mjd = None

        self._prepared = False

        self._filter = filt
        self._useFilter = False
        self._filteredData = None
        self._filteredDataReady = False

        if dataOrFilename:
            self.addData(dataOrFilename)
        
        if filt is not None:
            self.setFilter(filt)

        self.prepareData()

    def __str__(self):
        """Return a string representation of the class."""
        txt = f"oimData containing {len(self.data)} file(s)\n"
        for ifile, fi in enumerate(self.dataInfo):
            try:
                fname = Path(self.data[ifile].filename()).name
            except:
                fname = "None:"

            txt += f"{fname}\n"
            for di in fi:
                txt += f"\t{di['arr']}{di['nB']}\n"

                for ddi in di["data"]:
                    txt += f"\t\t{ddi}\n"
        return txt

    @property
    def data(self) -> None:
        """Return the data."""
        if not self._useFilter or self._filter is None:
            return self._data
        else:
            if not self._filteredDataReady:
                self.applyFilter()
            return self._filteredData

    def addData(self, dataOrFilename, prepare: Optional[bool] = True) -> None:
        """Add data to the class.

        Parameters
        ----------
        dataOrFilename : any
            The data to add. Can be a filename, a list of filenames, a hdulist
            or a list of hdulist.
        prepare : bool, optional
            Whether to prepare the data or not. The default is True.
        """
        self._data.extend(loadOifitsData(dataOrFilename))
        
        self.prepared = False
        self._filteredDataReady = False
        
        if prepare == True:
            self.prepareData()

    # TODO: Make it possible to selectively remove data
    def removeData(self, dataOrIndex: Any) -> None:
        """Remove either all data or specified data from the class.

        Parameters
        ----------
        dataOrIndex : any
            The data to remove. Can be a hdulist, a filename, a list of these
            or an index to remove.
        """
        self._prepared = False
        self._filteredDataReady = False
        self.prepareData()

    def setFilter(self, filt: Optional[oimDataFilter] = None,
                  useFilter: Optional[bool] = True) -> None:
        """Set the filter to use.

        Parameters
        ----------
        filt : oimDataFilter, optional
            The filter to use. The default is None.
        useFilter : bool, optional
            Whether to use the filter or not. The default is True.
        """
        # NOTE: Check type of filter
        if isinstance(filt, (list, oimDataFilterComponent)):
            filt = oimDataFilter(filt)
        
        self._filter = filt
        self._filteredDataReady = False
        self.useFilter = useFilter

    def applyFilter(self) -> None:
        """Apply the used filter(s) to the data."""
        self._filteredData = []

        for data in self._data:
            self._filteredData.append(hdulistDeepCopy(data))

        if self._filter is not None:
            self._filter.applyFilter(self._filteredData)

        self._filteredDataReady = True
        self.prepareData()

    @property
    def useFilter(self) -> None:
        """Return whether the filter is used or not."""
        return self._useFilter

    @useFilter.setter
    def useFilter(self, val) -> None:
        """Set whether to use the filter or not."""
        self._useFilter = val
        if val:
            if not self._filteredDataReady:
                self.applyFilter()

    def _analyzeOIFitFile(self, data: List[fits.HDUList]) -> None:
        """Analyze the oifits file and get the data info.

        Parameters
        ----------
        data : list of astropy.io.fits.HDUList
            The data to analyze.
        """
        dataInfo = []
        for datai in data :
            dataInfoi=[]
            for iarr, arri in enumerate(datai):
                info = None
                if arri.name in _oimDataTypeArr:
                    info = {"arr": arri.name, "idx": iarr}
                    if arri.name == "OI_VIS2":
                        nB=np.shape(arri.data["VIS2DATA"])
                    if arri.name == "OI_VIS":
                        nB = np.shape(arri.data["VISAMP"])
                    if arri.name == "OI_T3":
                        nB = np.shape(arri.data["T3AMP"])
                    if arri.name == "OI_FLUX":
                        try:
                            nB = np.shape(arri.data["FLUXDATA"])
                        except Exception:
                            nB = np.shape(arri.data["FLUX"])
                    info["nB"] = nB
                    info["data"] = oimDataCheckData(arri)
                if info:
                    dataInfoi.append(info)
            dataInfo.append(dataInfoi)
        self.dataInfo = dataInfo

    def prepareData(self) -> None:
        """Prepare the data for further analysis."""
        self._analyzeOIFitFile(self.data)
        self.vect_u = np.array([])
        self.vect_v = np.array([])
        self.vect_wl = np.array([])
        self.vect_dwl = np.array([])
        self.vect_mjd = np.array([])

        self.struct_u = []
        self.struct_v = []
        self.struct_wl = []
        self.struct_dwl = []
        self.struct_mjd = []
        self.struct_nB = []
        self.struct_nwl = []
        self.struct_val = []
        self.struct_err = []
        self.struct_flag = []
        self.struct_arrNum = []
        self.struct_arrType = []
        self.struct_dataType = []

        for _, datai in enumerate(self.data):
            self.struct_u.append([])
            self.struct_v.append([])
            self.struct_wl.append([])
            self.struct_dwl.append([])
            self.struct_mjd.append([])
            self.struct_nB.append([])
            self.struct_nwl.append([])
            self.struct_arrNum.append([])
            self.struct_arrType.append([])
            self.struct_dataType.append([])
            self.struct_val.append([])
            self.struct_err.append([])
            self.struct_flag.append([])

            for iarr, arri in enumerate(datai):
                if arri.name in _oimDataTypeArr:
                    dataTypeFlag, val, err, flag = oimGetDataValErrAndTypeFlag(arri)

                    if dataTypeFlag != oimDataType.NONE:
                        u, v, wl, dwl, mjd, nB, nwl = oimDataGetVectCoord(
                            datai, arri)

                        self.vect_u = np.concatenate((self.vect_u, u))
                        self.vect_v = np.concatenate((self.vect_v, v))
                        self.vect_wl = np.concatenate((self.vect_wl, wl))
                        self.vect_dwl = np.concatenate((self.vect_dwl, dwl))
                        self.vect_mjd = np.concatenate((self.vect_mjd, mjd))

                        self.struct_u[-1].append(u)
                        self.struct_v[-1].append(v)
                        self.struct_wl[-1].append(wl)
                        self.struct_dwl[-1].append(dwl)
                        self.struct_mjd[-1].append(mjd)
                        self.struct_nB[-1].append(nB)
                        self.struct_nwl[-1].append(nwl)

                    else:
                        self.struct_u[-1].append(np.array([]))
                        self.struct_v[-1].append(np.array([]))
                        self.struct_wl[-1].append(np.array([]))
                        self.struct_dwl[-1].append(np.array([]))
                        self.struct_mjd[-1].append(np.array([]))
                        self.struct_nB[-1].append(0)
                        self.struct_nwl[-1].append(0)

                    self.struct_arrNum[-1].append(iarr)
                    self.struct_arrType[-1].append(arri.name)
                    self.struct_dataType[-1].append(dataTypeFlag)
                    self.struct_val[-1].append(val)
                    self.struct_err[-1].append(err)
                    self.struct_flag[-1].append(flag)
        self._prepared = True

    def writeto(self, filename: Optional[Union[str, Path]] = None,
                overwrite: Optional[bool] = False,
                directory: Optional[Union[str, Path]] = None) -> None:
        """Write the data to a file.

        Parameters
        ----------
        filename : str or pathlib.Path, optional
            The filename to write the data to. The default is None.
        overwrite : bool, optional
            Whether to overwrite the file or not. The default is False.
        directory : str or pathlib.Path, optional
            The directory to write the data to. The default is None.
        """
        ndata = len(self.data)
        for idata, datai in enumerate(self.data):
            if filename is not None:
                if ndata == 1:
                    filenamei = Path(filename)
                else:
                    filenamei = Path(f"{Path(filename).stem}_{idata}.fits")
            elif directory is not None and datai.filename() is not None:
                filenamei = Path(directory) / Path(datai.filename()).name
            elif datai.filename() is not None:
                filenamei = Path(datai.filename())
            else:
                raise TypeError("Can't save the data!")
            try:
                datai.writeto(filenamei, overwrite=overwrite)
            except:
                raise TypeError("Can't save the data!")
