# -*- coding: utf-8 -*-
"""Various utilities for optical interferometry"""
import astropy.units as units
import numpy as np
from astropy.coordinates import Angle
from astropy.io import fits
from astroquery.simbad import Simbad


def getBaselineName(oifits, hduname="OI_VIS2", length=False, angle=False,
                    extver=None, squeeze=True):
    """Gets the baseline names, i.e., telescopes names
    separated by minus sign, in an extension of a oifits file.

    By default it is reading the 'OI_VIS' extension

    Parameters
    ----------
    oifits: astropy.io.fits.hdu.hdulist.HDUList
        An oifits file structure already opened with astropy.io.fits
    hduname: str, optional
        The fits extension name. The default is "OI_VIS2".
    length: bool, optional
        Add baseline length to the returned result. The default is False.
    angle: bool, optional
        Add baseline position angle ((in deg)à=) to the returned result.
        The default is False

    Returns
    -------
    name:  python list of str
        The array containing the baseline names (or triplet) and optionally
        the baseline length and orientation.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header['ARRNAME'] for i in idx_arr])

    if extver != None:
        data_arrnames = [data[hduname, extver].header['ARRNAME']]
    else:
        idx = np.where(extnames == hduname)[0]
        data_arrnames = [data[i].header['ARRNAME'] for i in idx]

    name = []

    for idata, data_arrname in enumerate(data_arrnames):

        iarr = idx_arr[np.where(arrnames == data_arrname)[0][0]]

        stanames = data[iarr].data['STA_NAME']
        staindexes = data[iarr].data['STA_INDEX']

        staidx = data[idx[idata]].data['STA_INDEX']
        shape = np.shape(staidx)
        namei = []
        if length or angle and hduname != "OI_T3":
            u = data[idx[idata]].data['UCOORD']
            v = data[idx[idata]].data['VCOORD']
            B = np.sqrt(u**2+v**2)
            PA = np.rad2deg(np.arctan2(u, v))
        for i in range(shape[0]):
            namej = ""
            for j in range(shape[1]):
                namej += stanames[np.where(staindexes == staidx[i, j])[0]][0]
                if j < shape[1]-1:
                    namej += "-"
            if hduname != "OI_T3":
                if length:
                    namej += " {0:.0f}m".format(B[i])
                if angle:
                    namej += " {0:.0f}$^o$".format(PA[i])
            namei.append(namej)

        name.append(namei)

    if squeeze == True and len(name) == 1:
        name = name[0]
    return name


def getConfigName(oifits, hduname="OI_VIS2", extver=None, squeeze=True):
    # TODO : Add support for multiple extensions
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header['ARRNAME'] for i in idx_arr])

    if extver != None:
        data_arrnames = [data[hduname, extver].header['ARRNAME']]
    else:
        idx = np.where(extnames == hduname)[0]
        data_arrnames = [data[i].header['ARRNAME'] for i in idx]

    name = []
    for idata, data_arrname in enumerate(data_arrnames):

        iarr = idx_arr[np.where(arrnames == data_arrname)[0][0]]

        stanames = data[iarr].data['STA_NAME']
        staindexes = data[iarr].data['STA_INDEX']

        staidx = np.unique(data[idx[idata]].data['STA_INDEX'].flatten())
        s = staidx.size
        namei = ""
        for i in range(s):
            namei += stanames[np.where(staindexes == staidx[i])[0]][0]
            if i < s-1:
                namei += "-"
        name.append(namei)

    if squeeze == True and len(name) == 1:
        name = name[0]
    return name


def getBaselineLengthAndPA(oifits, arr="OI_VIS2", extver=None, squeeze=True):
    """Return a tuple (B, PA) of the baseline lengths and orientation
    (position angles) from a fits extension within an opened oifits file.

    By default it is reading the OI_VIS extension.

    Parameters
    ----------
    oifits: astropy.io.fits.hdu.hdulist.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr: str, optional
        The fits extension name. The default is "OI_VIS2".

    Returns
    -------
    B: numpy.ndarray
        the array containing the baselines length.
    PA: numpy.ndarray
        the array containing the baselines orientation (in deg).
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    if extver != None:
        data = [data[arr, extver]]
    else:
        extnames = np.array([di.name for di in data])
        idx = np.where(extnames == arr)[0]
        data = [data[i] for i in idx]

    B = []
    PA = []
    for idata, datai in enumerate(data):

        if arr != "OI_T3":
            u = datai.data["UCOORD"]
            v = datai.data["VCOORD"]

            B.append(np.sqrt(u**2+v**2))
            PA.append(np.rad2deg(np.arctan2(u, v)))
        else:
            u1 = datai.data["U1COORD"]
            v1 = datai.data["V1COORD"]
            u2 = datai.data["U2COORD"]
            v2 = datai.data["V2COORD"]
            u3 = u1+u2
            v3 = v1+v2
            B1 = np.sqrt(u1**2+v1**2)
            B2 = np.sqrt(u2**2+v2**2)
            B3 = np.sqrt(u3**2+v3**2)
            B.append(np.array([B1, B2, B3]))
            PA1 = np.rad2deg(np.arctan2(u1, v1))
            PA2 = np.rad2deg(np.arctan2(u2, v2))
            PA3 = np.rad2deg(np.arctan2(u3, v3))
            PA.append(np.array([PA1, PA2, PA3]))
    if squeeze == True and len(B) == 1:
        B = B[0]
        PA = PA[0]
    return B, PA


def getSpaFreq(oifits, arr="OI_VIS2", unit=None, extver=None, squeeze=True):
    """

    Parameters
    ----------
    oifits: TYPE
        DESCRIPTION.
    arr: TYPE, optional
        DESCRIPTION. The default is "OI_VIS2".
    unit: TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    spaFreq: TYPE
        DESCRIPTION.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    B, _ = getBaselineLengthAndPA(data, arr, extver, squeeze=False)

    if arr == "OI_T3":
        B = [np.max(Bi, axis=0) for Bi in B]

    extnames = np.array([di.name for di in data])

    if extver != None:
        arrays = [data[arr, extver]]
        insnames = np.array([arrays.header['INSNAME']])
    else:
        idx = np.where(extnames == arr)[0]
        insnames = [data[i].header['INSNAME'] for i in idx]
        arrays = [data[i] for i in idx]

    idx_wlarr = np.where(extnames == "OI_WAVELENGTH")[0]
    wl_insnames = np.array([data[i].header['INSNAME'] for i in idx_wlarr])

    if unit == "cycles/mas":
        mult = units.mas.to(units.rad)
    elif unit == "cycles/arcsec":
        mult = units.arcsec.to(units.rad)
    elif unit == "Mlam":
        mult = 1/(1e6)
    else:
        mult = 1

    spaFreq = []

    for iarr, arri in enumerate(arrays):

        iwlarr = idx_wlarr[np.where(wl_insnames == insnames[iarr])[0][0]]

        lam = data[iwlarr].data["EFF_WAVE"]
        nlam = np.size(lam)
        nB = np.size(B[iarr])

        spaFreqi = np.ndarray([nB, nlam])
        for iB in range(nB):
            spaFreqi[iB, :] = B[iarr][iB]/lam*mult

        spaFreq.append(spaFreqi)

    if squeeze == True and len(spaFreq) == 1:
        spaFreq = spaFreq[0]

    return spaFreq


def hdulistDeepCopy(hdulist):
    res = hdulist.copy()
    res._file = hdulist._file
    for iext, exti in enumerate(res):
        res[iext] = exti.copy()
        res[iext].header = exti.header.copy()

    return res


_cutArr = ['EFF_WAVE', 'EFF_BAND', 'VIS2DATA', 'VIS2ERR', 'FLAG', 'VISAMP', 'VISAMPERR',
           'VISPHI', 'VISPHIERR', 'T3AMP', 'T3AMPERR', 'T3PHI', 'T3PHIERR',
           'FLUXDATA', 'FLUXERR', 'FLAG']


def cutWavelengthRange(oifits, wlRange=None, addCut=[]):
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits
    extnames = np.array([data[i].name for i in range(len(data))])

    oiwl_idx = np.where(extnames == "OI_WAVELENGTH")[0]

    wlRange = np.array(wlRange)

    if wlRange.ndim == 1:
        wlRange = wlRange.reshape((1, len(wlRange)))

    for i in oiwl_idx:
        insname = data[i].header['INSNAME']

        idx_wl_cut = []

        for wlRangei in wlRange:
            idx_wl_cut.extend(np.where((data[i].data['EFF_WAVE'] >= wlRangei[0]) &
                                       (data[i].data['EFF_WAVE'] <= wlRangei[1]))[0])

        idx_wl_cut = np.sort(idx_wl_cut)
        nwl_cut = len(idx_wl_cut)
        for idata, datai in enumerate(data):
            if "INSNAME" in datai.header:
                if datai.header['INSNAME'] == insname:
                    colDefs = []
                    for col in datai.columns:
                        if np.isin(col.name, _cutArr) or np.isin(col.name, addCut):
                            format0 = col.format[-1]
                            shape = datai.data[col.name].shape
                            if len(shape) == 2:
                                arr = np.take(
                                    datai.data[col.name], idx_wl_cut, axis=1)
                                colDefs.append(fits.Column(name=col.name, format="{0}{1}".
                                                           format(nwl_cut, format0), array=arr, unit=col.unit))
                            else:
                                arr = np.take(datai.data[col.name], idx_wl_cut)
                                colDefs.append(fits.Column(name=col.name,
                                                           format=col.format, array=arr,
                                                           unit=col.unit))

                        else:
                            colDefs.append(fits.Column(name=col.name,
                                                       format=col.format, array=datai.data[col.name],
                                                       unit=col.unit))

                    cols = fits.ColDefs(colDefs)
                    hdu = fits.BinTableHDU.from_columns(cols)
                    hdu.header = datai.header
                    hdu.update()
                    data[idata] = hdu
    return data



def getWlFromFitsImageCube(header, outputUnit=None):
    """Returns the wl law from a chromatic cube image in the fits format

    Parameters
    ----------
    header: astropy.io.fits.header
        DESCRIPTION.
    outputUnit: astropy.unit, optional
        If set convert the result to the proper unit. The default is None.

    Returns
    -------
    wl: float
        The wavelength in the given unit of the fits cube or the units specified
        by the user if outputUnit is set

    """
    dwl = header['CDELT3']
    nwl = header['NAXIS3']
    wl0 = header['CRVAL3']
    try:
        x0 = header['CRPIX3']
    except:
        x0 = 0
    wl = wl0+(np.arange(nwl)-x0)*dwl

    if outputUnit:
        if "CUNIT3" in header:
            try:
                unit0 = units.Unit(header["CUNIT3"])
            except:
                unit0 = units.m
        else:
            unit0 = units.m
        wl*unit0.to(outputUnit)

    return wl


def createOiArray(arrname, arrx, arry, arrz, sta_name, tel_name, diameter, staxyz):
    """

    Parameters
    ----------
    arrname: TYPE
        DESCRIPTION
    arrx: TYPE
        DESCRIPTION
    arry: TYPE
        DESCRIPTION
    arrz: TYPE
        DESCRIPTION
    sta_name: TYPE
        DESCRIPTION
    tel_name: TYPE
        DESCRIPTION
    diameter: TYPE
        DESCRIPTION
    staxyz: TYPE
        DESCRIPTION

    Returns
    -------
    arr: TYPE
        DESCRIPTION
    """
    nstation = np.size(sta_name)
    tel_name = fits.Column(name="TEL_NAME", format="A16",
                           array=np.array(tel_name))
    sta_name = fits.Column(name="STA_NAME", format="A16",
                           array=np.array(sta_name))
    sta_index = fits.Column(name="STA_INDEX", format="I2",
                            array=np.arange(1, nstation+1))
    diameter = fits.Column(name="DIAMETER", format="E",
                           array=np.array(diameter), unit='m')
    staxyz = fits.Column(name="STAXYZ", format="3D",
                         array=np.array(staxyz), unit='m')

    cols = [tel_name, sta_name, sta_index, diameter, staxyz]
    arr = fits.BinTableHDU.from_columns(cols)

    arr.header['EXTVER'] = (1, 'ID number of this OI_ARRAY')
    arr.header['ARRAYX'] = float(arrx)  # ,'[m] Array center X coordinate')
    arr.header['ARRAYY'] = float(arry)  # ,'[m] Array center Y coordinate')
    arr.header['ARRAYZ'] = float(arrz)  # ,'[m] Array center Z coordinate')
    arr.header['FRAME'] = ('GEOCENTRIC', 'Coordinate frame')
    arr.header['EXTNAME'] = 'OI_ARRAY'
    arr.header['OI_REVN'] = (1, 'Revision number of the table definition')
    arr.header['ARRNAME'] = (arrname, 'Array name')

    return arr


def createOiTargetFromSimbad(names):
    """

    Parameters
    ----------
    names : TYPE
        DESCRIPTION.

    Returns
    -------
    tar : TYPE
        DESCRIPTION.
    """
    customSimbad = Simbad()
    customSimbad.add_votable_fields(
        'plx', 'plx_error', 'propermotions', 'sptype', 'velocity')
    if type(names) == type(""):
        names = [names]
    data = customSimbad.query_objects(names)

    ntargets = len(names)
    rad = Angle(data['RA'], unit="hourangle").deg
    dec = Angle(data['DEC'], unit="deg").deg
    ra_err = (data['COO_ERR_MAJA'].data*units.mas).to_value(unit='deg')
    dec_err = (data['COO_ERR_MINA'].data*units.mas).to_value(unit='deg')
    pmra = (data['PMRA'].data*units.mas).to_value(unit='deg')
    pmdec = (data['PMDEC'].data*units.mas).to_value(unit='deg')
    pmra_err = (data['PM_ERR_MAJA'].data*units.mas).to_value(unit='deg')
    pmdec_err = (data['PM_ERR_MINA'].data*units.mas).to_value(unit='deg')
    plx_value = (data['PLX_VALUE'].data*units.mas).to_value(unit='deg')
    plx_error = (data['PLX_ERROR'].data*units.mas).to_value(unit='deg')

    target_id = fits.Column(name="TARGET_ID", format="I",
                            array=np.arange(1, ntargets+1))
    target = fits.Column(name="TARGET", format="16A", array=names)
    raep0 = fits.Column(name="RAEP0", format="D", array=rad, unit="deg")
    decep0 = fits.Column(name="DECEP0", format="D", array=dec, unit="deg")
    equinox = fits.Column(name="EQUINOX", format="E",
                          array=np.repeat(2000, ntargets), unit="yr")
    ra_err = fits.Column(name="RA_ERR", format="D", array=ra_err, unit="deg")
    dec_err = fits.Column(name="DEC_ERR", format="D",
                          array=dec_err, unit="deg")
    sysvel = fits.Column(name="SYSVEL", format="D",
                         array=np.zeros([ntargets]), unit="m/s")  # TODO
    veltyp = fits.Column(name="VELTYP", format="8A",
                         array=np.repeat("UNKNOWN", ntargets))  # TODO
    veldef = fits.Column(name="VELDEF", format="8A",
                         array=np.repeat("OPTICAL", ntargets))  # TODO
    pmra = fits.Column(name="PMRA", format="D", array=pmra, unit="deg/yr")
    pmdec = fits.Column(name="PMDEC", format="D", array=pmdec, unit="deg/yr")
    pmra_err = fits.Column(name="PMRA_ERR", format="D",
                           array=pmra_err, unit="deg/yr")
    pmdec_err = fits.Column(name="PMDEC_ERR", format="D",
                            array=pmdec_err, unit="deg/yr")
    parallax = fits.Column(name="PARALLAX", format="E",
                           array=plx_value, unit="deg")
    para_err = fits.Column(name="PARA_ERR", format="E",
                           array=plx_error, unit="deg")
    spectyp = fits.Column(name="SPECTYP", format="16A", array=data['SP_TYPE'])

    cols = [target_id, target, raep0, decep0, equinox, ra_err, dec_err, sysvel, veltyp,
            veldef, pmra, pmdec, pmra_err, pmdec_err, parallax, para_err, spectyp]

    tar = fits.BinTableHDU.from_columns(cols)

    tar.header['OI_REVN'] = (1, 'Revision number of the table definition')
    tar.header['EXTVER'] = (1, 'ID number of this OI_TARGET')
    tar.header['EXTNAME'] = 'OI_TARGET'
    return tar


def createOiWavelength(insname, eff_wave, eff_band):
    """

    Parameters
    ----------
    insname : TYPE
        DESCRIPTION
    eff_wave : TYPE
        DESCRIPTION
    eff_band : TYPE
        DESCRIPTION

    Returns
    -------
    wave : TYPE
        DESCRIPTION

    """
    eff_wave = fits.Column(name="EFF_WAVE", format="E",
                           array=np.array(eff_wave), unit='m')
    eff_band = fits.Column(name="EFF_BAND", format="E",
                           array=np.array(eff_band), unit='m')

    cols = [eff_wave, eff_band]
    wave = fits.BinTableHDU.from_columns(cols)
    wave.header['EXTNAME'] = 'OI_WAVELENGTH'
    wave.header['EXTVER'] = (1, 'ID number of this OI_WAVELENGTH')
    wave.header['OI_REVN'] = (1, 'Revision number of the table definition')
    wave.header['INSNAME'] = (
        insname, 'Identifies corresponding OI_WAVELENGTH')

    return wave

###############################################################################


def createOiVis2(arrname, insname, target_id, time, mjd, int_time, vis2data, vis2err,
                 ucoord, vcoord, sta_index, flag, dateobs):
    """

    Parameters
    ----------
    arrname: TYPE
        DESCRIPTION
    insname: TYPE
        DESCRIPTION
    target_id: TYPE
        DESCRIPTION
    time: TYPE
        DESCRIPTION
    mjd: TYPE
        DESCRIPTION
    int_time: TYPE
        DESCRIPTION
    vis2data: TYPE
        DESCRIPTION
    vis2err: TYPE
        DESCRIPTION
    ucoord: TYPE
        DESCRIPTION
    vcoord: TYPE
        DESCRIPTION
    sta_index: TYPE
        DESCRIPTION
    flag: TYPE
        DESCRIPTION
    dateobs: TYPE
        DESCRIPTION

    Returns
    -------
    oivis2: TYPE
        DESCRIPTION
    """
    nb = len(target_id)

    v2shape = np.shape(vis2data)
    v2dim = len(v2shape)

    if (v2dim == 2):
        nlam = v2shape[1]
    elif (v2dim == 1):
        if v2shape[0] == nb:
            nlam = 1
        else:
            nlam = v2shape[0]

    target_id = fits.Column(name="TARGET_ID", format="I",
                            array=np.array(target_id))
    time = fits.Column(name="TIME", format="D",
                       array=np.array(time), unit="sec")
    mjd = fits.Column(name="MJD", format="D", array=np.array(mjd), unit="day")
    int_time = fits.Column(name="INT_TIME", format="D",
                           array=np.array(int_time), unit="sec")
    vis2data = fits.Column(name="VIS2DATA", format="{0}D".format(
        nlam), array=np.array(vis2data))
    vis2err = fits.Column(name="VIS2ERR", format="{0}D".format(
        nlam), array=np.array(vis2err))
    ucoord = fits.Column(name="UCOORD", format="1D",
                         array=np.array(ucoord), unit="m")
    vcoord = fits.Column(name="VCOORD", format="1D",
                         array=np.array(vcoord), unit="m")
    sta_index = fits.Column(name="STA_INDEX", format="2I",
                            array=np.array(sta_index))
    flag = fits.Column(name="FLAG", format="{0}L".format(
        nlam), array=np.array(flag))

    cols = [target_id, time, mjd, int_time, vis2data,
            vis2err, ucoord, vcoord, sta_index, flag]

    oivis2 = fits.BinTableHDU.from_columns(cols)
    oivis2.header['EXTNAME'] = 'OI_VIS2'
    oivis2.header['EXTVER'] = (1, 'ID number of this OI_VIS2')
    oivis2.header['OI_REVN'] = (1, 'Revision number of the table definition')
    oivis2.header['INSNAME'] = (
        insname, 'Identifies corresponding OI_WAVELENGTH')
    oivis2.header['DATE-OBS'] = dateobs
    oivis2.header['ARRNAME'] = arrname

    return oivis2


def createOiVis(arrname, insname, target_id, time, mjd, int_time, visamp, visamperr, visphi, visphierr,
                ucoord, vcoord, sta_index, flag, dateobs, amptyp="absolute", phityp="absolute"):
    nb = len(target_id)

    vshape = np.shape(visamp)
    vdim = len(vshape)

    if (vdim == 2):
        nlam = vshape[1]
    elif (vdim == 1):
        if vshape[0] == nb:
            nlam = 1
        else:
            nlam = vshape[0]

    target_id = fits.Column(name="TARGET_ID", format="I",
                            array=np.array(target_id))
    time = fits.Column(name="TIME", format="D",
                       array=np.array(time), unit="sec")
    mjd = fits.Column(name="MJD", format="D", array=np.array(mjd), unit="day")
    int_time = fits.Column(name="INT_TIME", format="D",
                           array=np.array(int_time), unit="sec")
    visamp = fits.Column(name="VISAMP", format="{0}D".format(
        nlam), array=np.array(visamp))
    visamperr = fits.Column(name="VISAMPERR", format="{0}D".format(
        nlam), array=np.array(visamperr))
    visphi = fits.Column(name="VISPHI", format="{0}D".format(
        nlam), array=np.array(visphi))
    visphierr = fits.Column(name="VISPHIERR", format="{0}D".format(
        nlam), array=np.array(visphierr))
    ucoord = fits.Column(name="UCOORD", format="1D",
                         array=np.array(ucoord), unit="m")
    vcoord = fits.Column(name="VCOORD", format="1D",
                         array=np.array(vcoord), unit="m")
    sta_index = fits.Column(name="STA_INDEX", format="2I",
                            array=np.array(sta_index))
    flag = fits.Column(name="FLAG", format="{0}L".format(
        nlam), array=np.array(flag))

    cols = [target_id, time, mjd, int_time, visamp, visamperr,
            visphi, visphierr, ucoord, vcoord, sta_index, flag]

    oivis = fits.BinTableHDU.from_columns(cols)
    oivis.header['EXTNAME'] = 'OI_VIS'
    oivis.header['EXTVER'] = (1, 'ID number of this OI_VIS')
    oivis.header['OI_REVN'] = (1, 'Revision number of the table definition')
    oivis.header['INSNAME'] = (
        insname, 'Identifies corresponding OI_WAVELENGTH')
    oivis.header['DATE-OBS'] = dateobs
    oivis.header['ARRNAME'] = arrname
    oivis.header['AMPTYP'] = amptyp
    oivis.header['PHITYP'] = phityp

    return oivis
    int_time = fits.Column(name="INT_TIME", format="D",
                           array=np.array(int_time), unit="sec")
    visamp = fits.Column(name="VISAMP", format="{0}D".format(
        nlam), array=np.array(visamp))
    visamperr = fits.Column(name="VISAMPERR", format="{0}D".format(
        nlam), array=np.array(visamperr))
    visphi = fits.Column(name="VISPHI", format="{0}D".format(
        nlam), array=np.array(visphi))
    visphierr = fits.Column(name="VISPHIERR", format="{0}D".format(
        nlam), array=np.array(visphierr))
    ucoord = fits.Column(name="UCOORD", format="1D",
                         array=np.array(ucoord), unit="m")
    vcoord = fits.Column(name="VCOORD", format="1D",
                         array=np.array(vcoord), unit="m")
    sta_index = fits.Column(name="STA_INDEX", format="2I",
                            array=np.array(sta_index))
    flag = fits.Column(name="FLAG", format="{0}L".format(
        nlam), array=np.array(flag))

    cols = [target_id, time, mjd, int_time, visamp, visamperr,
            visphi, visphierr, ucoord, vcoord, sta_index, flag]

    oivis = fits.BinTableHDU.from_columns(cols)
    oivis.header['EXTNAME'] = 'OI_VIS'
    oivis.header['EXTVER'] = (1, 'ID number of this OI_VIS')
    oivis.header['OI_REVN'] = (1, 'Revision number of the table definition')
    oivis.header['INSNAME'] = (
        insname, 'Identifies corresponding OI_WAVELENGTH')
    oivis.header['DATE-OBS'] = dateobs
    oivis.header['ARRNAME'] = arrname
    oivis.header['AMPTYP'] = amptyp
    oivis.header['PHITYP'] = phityp

    return oivis
