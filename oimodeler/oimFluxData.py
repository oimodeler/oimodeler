# -*- coding: utf-8 -*-
"""Photometric and spectrascopic data wrappers"""

from astropy.io import fits
import numpy as np
from .oimUtils import (
    createOiFlux,
    createOiWavelength,
    createOiArray,
    createOiTargetFromSimbad,
)


class oimFluxData(fits.HDUList):
    """Base class for photometric and spectroscopic data wrappers into hdulist"""

    def __init__(self, target, wl, dwl, flx, flxerr,
                 insname="DUMMY", array="DUMMY", dateobs="1980-04-23T20:15:00",
                 mjd=None, int_time=None, sta_index=[0], unit="Jy", flag=None):
        super().__init__()

        oi_tar = self._createOiTarget(target)
        oi_arr = self._createOiArray(array)
        arrname = oi_arr.header["ARRNAME"]
        oi_wl = self._createOiWavelength(wl, dwl, insname)
        oi_flx = self._createOiFlux(
                flx, flxerr, arrname=arrname, insname=insname,
                mjd=mjd, int_time=int_time, sta_index=sta_index,
                unit=unit, flag=flag, dateobs=dateobs)

        self._setHdulist(None, oi_tar, oi_arr, oi_wl, oi_flx)

    def _setHdulist(self, primary, oi_tar, oi_arr, oi_wl, oi_flx):
        if not (primary):
            primary = fits.PrimaryHDU()
        self.append(primary)
        self.append(oi_tar)
        self.append(oi_arr)
        self.append(oi_wl)
        self.append(oi_flx)

    def _createOiArray(self, array):
        if array == "DUMMY":
            oi_arr = createOiArray(
                arrname="DUMMY", FRAME="GEOCENTRIC",
                arrayx=4580293.73, arrayy=572102.89,
                arrayz=4386963.85, tel_name=["DUMMY"],
                sta_name=["DUMMY"], sta_index=[0],
                DIAMETER=[0], STAXYZ=[0, 0, 0],
                FOV=[0], FOVTYPE=["FWHM"])
        else:
            oi_arr = array
        return oi_arr

    def _createOiTarget(self, target):
        if isinstance(target, str):
            oi_tar = createOiTargetFromSimbad(target)
        else:
            oi_tar = target
        return oi_tar

    def _createOiWavelength(self, wl, dwl, insname):
        wl = np.array(wl)
        nwl = np.size(wl)
        if np.size(dwl) == 1:
            dwl = np.full(nwl, dwl)
        return createOiWavelength(insname=insname, eff_wave=wl, eff_band=dwl)

    def _createOiFlux(self, flx, flxerr, arrname="DUMMY",
                      insname="DUMMY", mjd=None, int_time=None, sta_index=[1],
                      unit="Jy", flag=None, dateobs=None):
        nflx = flx.size
        if not (mjd):
            mjd = np.zeros(nflx)
        if not (int_time):
            int_time = np.zeros(nflx)
        if not (flag):
            flag = np.full((nflx), True, dtype=bool)

        return createOiFlux(
            dateobs=dateobs,
            arrname=arrname,
            insname=insname,
            target_id=[0],
            mjd=mjd,
            int_time=int_time,
            fluxdata=flx,
            fluxerr=flxerr,
            sta_index=sta_index,
            flag=flag,
            unit=unit,
        )
