# -*- coding: utf-8 -*-
"""Various utilities for optical interferometry"""
import astropy.units as units
import numpy as np
from scipy.stats import circstd, circvar
from astropy.coordinates import Angle
from astropy.io import fits
from astroquery.simbad import Simbad
import oimodeler as oim
from os import PathLike


_oimDataType = ["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI", "FLUXDATA"]
_oimDataTypeErr = ["VIS2ERR", "VISAMPERR",
                   "VISPHIERR", "T3AMPERR", "T3PHIERR", "FLUXERR"]
_oimDataTypeArr = ["OI_VIS2", "OI_VIS", "OI_VIS", "OI_T3", "OI_T3", "OI_FLUX"]

_oimDataAnalysisInComplex = [False,False,True,False,True,False]

def getDataTypeIsAnalysisComplex(dataType):
    try:
        return _oimDataAnalysisInComplex[_oimDataType.index(dataType)]
    except:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")
        
def getDataArrname(dataType):
    try:
        return _oimDataTypeArr[_oimDataType.index(dataType)]
    except:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")
   
def getDataType(dataArrname):
    return[ datatypei for datatypei,arrnamei 
           in zip(_oimDataType,_oimDataTypeArr) if arrnamei==dataArrname]

def getDataTypeError(dataArrname):
    return[ datatypei for datatypei,arrnamei 
           in zip(_oimDataTypeErr,_oimDataTypeArr) if arrnamei==dataArrname]




def loadOifitsData(something, mode="listOfHdlulist"):
    """ 
    return the oifits data from either filenames, already opened oifts or a 
    oimData boject as either a list of hdlulist (default) or as a oimData 
    object using the option mode="oimData". This Function is used in oimData 
    and multiple plotting functions
    
    Parameters
    ----------
    something : astropy.io.fits.hdu.hdulist.HDUList, oimodeler.oimData, string, list
        The data to deal with. Can be a oimData object, a hdulist, a string 
        representing a filename or a list of these kind of object
    mode : str, optional
        The type of the return data, either "listOfHdlulist" or "oimData"
        The default is "listOfHdlulist"
    """
    
    if isinstance(something, oim.oimData):
        if mode =="oimData":
            data = something
        else:    
            data = something.data    
    else: 
               
        if isinstance(something,fits.hdu.hdulist.HDUList) or \
           isinstance(something,PathLike) or \
           isinstance(something,str):
            something=[something]
            
        if isinstance(something,list):
            data=[] 
            
            for el in something:
                if  isinstance(el,fits.hdu.hdulist.HDUList):
                    data.append(el)
                else:
                    try:
                        data.append(fits.open(el))
                    except:
                        raise ValueError("the path does not exist or is not a"\
                                         " valid fits files")
        else:
            raise TypeError("only oimData,hdulist,Path or string, or list of"\
                            " these kind of objects allowed ") 
            
        if mode =="oimData":
             data = oim.oimData(data)
    
    return data
    

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
        if hduname == "OI_FLUX":
            staidx = staidx[:,None]
        shape = np.shape(staidx)
        namei = []
        if length or angle and (hduname != "OI_T3" or hduname != "OI_FLUX"):
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


def getBaselineLengthAndPA(oifits, arr="OI_VIS2", extver=None, squeeze=True,
                           returnUV=False):
    """Return a tuple (B, PA) of the baseline lengths and orientation
    (position angles) from a fits extension within an opened oifits file.

    By default it is reading the OI_VIS extension.

    Parameters
    ----------
    oifits: astropy.io.fits.hdu.hdulist.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr: str, optional
        The fits extension name. The default is "OI_VIS2".
    returnUV: bool
        if True also return the u,v coordinates in m 
        the default is False

    Returns
    -------
    B: numpy.ndarray
        the array containing the baselines length.
    PA: numpy.ndarray
        the array containing the baselines orientation (in deg).
    ucoord: numpy.ndarray
        the array containing the u coordinate (in m)(optional)
    ucoord: numpy.ndarray
        the array containing the u coordinate (in m)(optional)        
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
    ucoord=[]
    vcoord=[]
    for idata, datai in enumerate(data):

        if arr != "OI_T3":
            u = datai.data["UCOORD"]
            v = datai.data["VCOORD"]
            
            ucoord.append(u)
            vcoord.append(v)
            
            B.append(np.sqrt(u**2+v**2))
            PA.append(np.rad2deg(np.arctan2(u, v)))
        else:
            u1 = datai.data["U1COORD"]
            v1 = datai.data["V1COORD"]
            u2 = datai.data["U2COORD"]
            v2 = datai.data["V2COORD"]
            u3 = u1+u2
            v3 = v1+v2
            
            ucoord.append([u1,u2,u3])
            vcoord.append([v1,v2,v3])  
            
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
        ucoord=ucoord[0]
        vcoord=vcoord[0]
        
    if returnUV:
        return B,PA,ucoord,vcoord
    else:
        return B, PA

def get2DSpaFreq(oifits, arr="OI_VIS2", unit=None, extver=None, squeeze=True):
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

    _ , _, ucoord, vcoord = getBaselineLengthAndPA(data, arr, extver, 
                                                   squeeze=False, returnUV=True)

    if arr == "OI_T3":
        raise TypeError("get2DSpaFreq does not accept OI_T3 extension")

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

    spaFreqU = []
    spaFreqV = []
    for iarr, arri in enumerate(arrays):

        iwlarr = idx_wlarr[np.where(wl_insnames == insnames[iarr])[0][0]]

        lam = data[iwlarr].data["EFF_WAVE"]
        nlam = np.size(lam)
        nB = np.size(ucoord[iarr])

        spaFreqUi = np.ndarray([nB, nlam])
        spaFreqVi = np.ndarray([nB, nlam])        
        for iB in range(nB):
            spaFreqUi[iB, :] = ucoord[iarr][iB]/lam*mult
            spaFreqVi[iB, :] = vcoord[iarr][iB]/lam*mult
        spaFreqU.append(spaFreqUi)
        spaFreqV.append(spaFreqVi)
    if squeeze == True and len(spaFreqU) == 1:
        spaFreqU = spaFreqU[0]
        spaFreqV = spaFreqV[0]
        
    return spaFreqU,spaFreqV



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


def getWlFromOifits(oifits, arr="OI_VIS2",extver=None,returnBand=False):
    
    if isinstance(arr,str):
        arr=oifits[arr,extver]
        
    insname=arr.header['INSNAME']
    
    oiwls = np.array([di for di in oifits if di.name == "OI_WAVELENGTH"])
    oiwls_insname = np.array([oiwli.header['INSNAME'] for oiwli in oiwls])
    
    iwl=np.where(oiwls_insname == insname)[0][0]
    oiwl = oiwls[iwl]
    
    wl  = oiwl.data['EFF_WAVE']
    
    if returnBand:
        dwl = oiwl.data['EFF_BAND']
        return wl,dwl
    else:
        return wl
    

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

        idx_wl_cut = np.sort(idx_wl_cut).astype("int64")
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


###############################################################################
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


###############################################################################
oi_target_keywords=[
    ("OI_REVN", False, "Revision number")]

oi_target_columns=[
    ("TARGET_ID", "I",   False, "Index number. Must be >=1", None),
    ("TARGET",    "16A", False, "Target name",None),
    ("RAEP0",     "D",   False, "RA at mean EQUINOX ","deg"),
    ("DECEP0",    "D",   False, "Dec at mean EQUINOX","deg"),
    ("EQUINOX",   "E",   False, "Equinox",None),
    ("RA_ERR",    "D",   False, "Error in RA","deg"),
    ("DEC_ERR",   "D",   False, "Error in Dec","deg"),
    ("SYSVEL",    "D",   False, "Systemic radial velocity","m/s"),
    ("VELTYP",   "8A",   False, "Reference for radial velocity:LSR, GEOCENTR...",None),
    ("VELDEF",   "8A",   False, "Definition of radial velocity:(OPTICAL,RADIO)",None),
    ("PMRA",     "D",    False, "Proper motion in RA","deg/yr"),
    ("PMDEC",    "D",    False, "Proper motion in Dec","deg/yr"),
    ("PMRA_ERR", "D",    False, "Error of proper motion in RA","deg/yr"),
    ("PMDEC_ERR","D",    False, "Error of proper motion in Dec","deg/yr"),
    ("PARALLAX", "E",    False, "Parallax","deg"),
    ("PARA_ERR", "E",    False, "Error in parallax ","deg"),
    ("SPECTYP",  "16A",  False, "Spectral type","deg"),    
    ("CATEGORY", "3A",   True,  "CAL or SCI",None)]    

oi_array_keywords=[
    ("OI_REVN", False, "Revision number"),
    ("ARRNAME", False, "A Array name, for cross-referencing"),
    ("FRAME",   False, "A Coordinate frame"),    
    ("ARRAYX",  False, "Array center x coordinates (m)"),
    ("ARRAYY",  False, "Array center y coordinates (m)"),
    ("ARRAYZ",  False, "Array center z coordinates (m)")]

oi_array_columns=[
    ("TEL_NAME","16A",False," Telescope name", None),
    ("STA_NAME","16A",False,"Station name",None),
    ("STA_INDEX","I",False,"Station number. Must be >=1",None),
    ("DIAMETER","E",False,"Element diameter","m"),
    ("STAXYZ","3D",False,"Station coordinates w.r.t. array center","m"),
    ("FOV","D",False," Photometric field of view","arcsec"),
    ("FOVTYPE","6A",False,"Model for FOV: FWHM or RADIUS",None)]

oi_wl_keywords=[
    ("OI_REVN",  False, "Revision number"),
    ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table")]

oi_wl_columns=[("EFF_WAVE","E",False,"Effective wavelength of channel", "m"),
               ("EFF_BAND","E",False,"Effective bandpass of channel","m")]

oi_vis2_keywords=[
    ("OI_REVN",  False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True,  "Identifies corresponding OI_CORR table")]

oi_vis2_columns=[
  ("TARGET_ID",         "I",    False, "Target number as index into OI_TARGET table","m"),
  ("TIME",              "D",    False, "Zero. For backwards compatibility",None),
  ("MJD",               "D",    False, "Modified Julian day","day"),
  ("INT_TIME",          "D",    False, "Integration time",None),
  ("VIS2DATA",          "NWLD", False, "Squared Visibility",None),
  ("VIS2ERR",           "NWLD", False, "Error in Squared Visibility",None),
  ("CORRINDX_VIS2DATA", "J",    True,  "Index into correlation matrix for 1st VIS2DATA element",None),
  ("UCOORD",            "D",    False, "U coordinate of the data","m"),
  ("VCOORD",            "D",    False, "V coordinate of the data","m"),
  ("STA_INDEX",         "2I",   False, "Station numbers contributing to the data",None),
  ("FLAG",              "NWLL", False, "Flag",None)]

oi_vis_keywords=[
    ("OI_REVN",  False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True,  "Identifies corresponding OI_CORR table"),
    ("AMPTYP",   True,  "absolute, differential, correlated flux"),
    ("PHITYP",   True,  "absolute, differential"),
    ("AMPORDER", True,  "Polynomial fit order for differential chromatic amplitudes"),
    ("PHIORDER", True,  "Polynomial fit order for differential chromatic phases")]

#TODO implement RIV VISREFMAP ...
oi_vis_columns=[
  ("TARGET_ID",       "I",    False, "Target number as index into OI_TARGET table","m"),
  ("TIME",            "D",    False, "Zero. For backwards compatibility",None),
  ("MJD",             "D",    False, "Modified Julian day","day"),
  ("INT_TIME",        "D",    False, "Integration time","s"),
  ("VISAMP",          "NWLD", False, "Visibility amplitude",None),
  ("VISAMPERR",       "NWLD", False, "Error in visibility amplitude",None),
  ("CORRINDX_VISAMP", "J",    True,  "Index into correlation matrix for 1st VISAMP element",None),
  ("VISPHI",          "NWLD", False, "Visibility phase","deg"),
  ("VISPHIERR",       "NWLD", False, "Error in visibility Phase","deg"),
  ("CORRINDX_VISPHI", "J",    True,  "Index into correlation matrix for 1st VISPHI element",None),
  ("UCOORD",          "D",    False, "U coordinate of the data","m"),
  ("VCOORD",          "D",    False, "V coordinate of the data","m"),
  ("STA_INDEX",       "2I",   False, "Station numbers contributing to the data",None),
  ("FLAG",            "NWLL", False, "Flag",None)]

oi_t3_keywords=[
    ("OI_REVN",  False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True,  "Identifies corresponding OI_CORR table")]

oi_t3_columns=[
  ("TARGET_ID",      "I",    False, "Target number as index into OI_TARGET table","m"),
  ("TIME",           "D",    False, "Zero. For backwards compatibility",None),
  ("MJD",            "D",    False, "Modified Julian day","day"),
  ("INT_TIME",       "D",    False, "Integration time","s"),
  ("T3AMP",          "NWLD", False, "Triple product amplitude",None),
  ("T3AMPERR",       "NWLD", False, "Error in triple product amplitude",None),
  ("CORRINDX_T3AMP", "J",    True,  "Index into correlation matrix for 1st T3AMP element",None),
  ("T3PHI",          "NWLD", False, "Triple Product Phase","deg"),
  ("T3PHIERR",       "NWLD", False, "Error in Triple Product Phase","deg"),
  ("CORRINDX_T3PHI", "J",    True,  "Index into correlation matrix for 1st T3PHI element",None),
  ("U1COORD",        "D",    False, "U coordinate of baseline AB of the triangle","m"),
  ("V1COORD",        "D",    False, "V coordinate of baseline AB of the triangle","m"),
  ("U2COORD",        "D",    False, "U coordinate of baseline BC of the triangle","m"),
  ("V2COORD",        "D",    False, "V coordinate of baseline BC of the triangle","m"),
  ("STA_INDEX",      "3I",   False, "Station numbers contributing to the data",None),
  ("FLAG",           "NWLL", False, "Flag",None)]

oi_flux_keywords=[
    ("OI_REVN",  False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True,  "Identifies corresponding OI_CORR table"),
    ("FOV",      True,  "Area on sky over which flux is integrated (arcsec)"),
    ("FOVTYPE",  True,  "Model for FOV: FWHM or RADIUS"),
    ("CALSTAT",  True,  "C: Spectrum is calibrated, U: uncalibrated")]

oi_flux_columns=[
  ("TARGET_ID",         "I",    False, "Target number as index into OI_TARGET table","m"),
  ("MJD",               "D",    False, "Modified Julian day","day"),
  ("INT_TIME",          "D",    False, "Integration time","s"),
  ("FLUXDATA",          "NWLD", False, "Flux in units of TUNITn","external"),
  ("FLUXERR",           "NWLD", False, "Corresponding flux error","external"),
  ("CORRINDX_FLUXDATA", "J",    True,  "Index into correlation matrix for 1st FLUXDATA element",None),
  ("STA_INDEX",         "I",   True, "Station number contributing to the data",None),
  ("FLAG",              "NWLL", False, "Flag",None)]

###############################################################################

def _createOiTab(extname,keywords_def,colums_def,dataTypeFromShape,**kwargs):
    keys={}

    for el in kwargs.keys():
        upper=el.upper()
        keys[upper]=kwargs[el]
        
    if "DATE_OBS"  in keys:
        keys["DATE-OBS"] = keys.pop("DATE_OBS")
        
    if "DATEOBS"  in keys:
        keys["DATE-OBS"] = keys.pop("DATEOBS")        

    try:
        nb = len(keys["TARGET_ID"])
    except:
        nb=1
        
    shape = np.shape(keys[dataTypeFromShape])
    dim = len(shape)

    if not("TIME"  in keys) and "TIME" :
        keys["TIME"] = np.zeros(nb)
    if not("OI_REVN" in keys):
        keys["OI_REVN"] = 2
         
         
    if (dim == 2):
        nwl = shape[1]
    elif (dim == 1):
        if shape[0] == nb:
            nwl = 1
        else:
            nwl = shape[0]
            
    cols = []

    for colname,form,optional,comment,unit in colums_def:

        if not(colname in keys) and optional==False:
            raise TypeError(f"Missing {colname} column")  
            
        elif colname in keys:
            if "NWL" in form:
                form= f"{nwl}{form[3:]}"
        
            if unit=="external":
                try:
                    unit=keys["UNIT"]
                except:
                    raise TypeError(f"missing unit for {colname}")   
            cols.append(fits.Column(name=colname, format=form,
                                    array=np.array(keys[colname]), unit=unit))
        
    hdu = fits.BinTableHDU.from_columns(cols)
    
    for keyword,optional,comment in keywords_def:
        
        if not(keyword in keys) and optional==False:
            raise TypeError(f"Missing {keyword} keyword")
            
        elif keyword in keys:
            hdu.header[keyword]=(keys[keyword], comment)
    
    hdu.header['EXTNAME'] = extname
    if "EXTVER" in keys:
        hdu.header['EXTVER'] = (keys["EXTVER"], f'ID number of this {extname}')
    return hdu,keys

###############################################################################
def createOiTarget(**kwargs):
    hdu,keys = _createOiTab("OI_TARGET",oi_target_keywords,oi_target_columns,
                            "TARGET_ID",**kwargs)
    return hdu

###############################################################################
def createOiArray(**kwargs):
    hdu,keys = _createOiTab("OI_ARRAY",oi_array_keywords,oi_array_columns,
                            "STA_INDEX",**kwargs)
    return hdu

###############################################################################
def createOiWavelength(**kwargs):
    hdu,keys = _createOiTab("OI_WAVELENGTH",oi_wl_keywords,oi_wl_columns,
                            "EFF_WAVE",**kwargs)
    return hdu

###############################################################################
def createOiVis2(**kwargs):
    hdu,_ = _createOiTab("OI_VIS2",oi_vis2_keywords,oi_vis2_columns,
                         "VIS2DATA",**kwargs)
    return hdu

###############################################################################
def createOiVis(**kwargs):
    hdu,_ = _createOiTab("OI_VIS",oi_vis_keywords,oi_vis_columns,
                         "VISAMP",**kwargs)
    return hdu

###############################################################################
def createOiT3(**kwargs):
    hdu,_ = _createOiTab("OI_T3",oi_t3_keywords,oi_t3_columns,
                         "T3AMP",**kwargs)
    return hdu

###############################################################################
def createOiFlux(**kwargs):
    hdu,keys = _createOiTab("OI_FLUX",oi_flux_keywords,oi_flux_columns,
                            "FLUXDATA",**kwargs)
    return hdu

###############################################################################
def createOiTargetFromSimbad(names):

    customSimbad = Simbad()
    customSimbad.add_votable_fields('plx', 'plx_error', 
                                    'propermotions', 'sptype', 'velocity')
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

    hdu = createOiTarget(target_id=np.arange(1, ntargets+1), target=names,
        raep0=rad, decep0=dec, equinox=np.repeat(2000, ntargets),
        ra_err = ra_err, dec_err = dec_err, sysvel=np.zeros([ntargets]),
        veltyp = np.repeat("UNKNOWN", ntargets), 
        veldef = np.repeat("OPTICAL", ntargets), pmra = pmra,pmdec = pmdec,
        pmra_err = pmra_err, pmdec_err = pmdec_err, parallax = plx_value,
        para_err = plx_error,spectyp = data['SP_TYPE'])

    return hdu


###############################################################################

#TODO this current method does not allow to shift baseline of the same oifits
# individually. 

def shiftWavelength(oifits,shift,verbose=False):
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits
    extnames=np.array([data[i].name for i in range(len(data))])

    wl_idx=np.where(extnames=="OI_WAVELENGTH")[0]
    for i in wl_idx:
        if verbose :
            print("OI_WAVELENGTH table at index {0}".format(i))
        insname=data[i].header['INSNAME']
        if verbose :
            print("INSNAME = {0}".format(insname))
        data[i].data['EFF_WAVE']+=shift

###############################################################################
#TODO for phases smoothing should be done in complex plan
def spectralSmoothing(oifits,kernsize,cols2Smooth="all",normalizeError=True):
    
    tableToSmooth=["OI_VIS","OI_VIS2","OI_T3","OI_FLUX"]
    
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits

    kernel=np.ones(kernsize)/kernsize
    
    if not( isinstance(cols2Smooth,list)):
        cols2Smooth=[cols2Smooth]
        
    if  "all" in cols2Smooth:
        cols2Smooth = ["VIS2DATA","VIS2ERR","VISAMP","VISAMPERR","VISPHI","VISPHIERR",
                     "T3AMP","T3AMPERR","T3PHI","T3PHIERR","FLUXDATA","FLUXDATAERR"]
        
    for i in range(1,len(data)):
        try:
            if data[i].name in tableToSmooth:
                cols = data[i].data.columns
                for coli in cols:
                    if coli.name in cols2Smooth:
                        dims=np.shape(data[i].data[coli.name])
    
                        if len(dims)==2:
                            nB=dims[0]
                            for iB in range(nB):
                                data[i].data[coli.name][iB,:]=np.convolve(
                                    data[i].data[coli.name][iB,:],kernel,"same")
                        else:
                              data[i].data[coli.name]=np.convolve(
                                  data[i].data[coli.name],kernel,"same")
                              
                        if normalizeError == True and coli.name in _oimDataTypeErr:
                            data[i].data[coli.name]/=np.sqrt(kernsize)
        except:
            pass

###############################################################################

def _rebin(arr, binsize,median=True):
    newsize = (arr.shape[0] // int(binsize)) * binsize
    arr = arr [:newsize]
    shape = (arr.shape[0]// binsize, binsize)
    if median:
        res=np.median(arr.reshape(shape),axis=-1)
    else:
        res=arr.reshape(shape).mean(-1)
    return res


def _rebinHdu(hdu,binsize,exception=[]):
    cols = hdu.data.columns
    newcols=[]

    dim2 = 2 in [len(np.shape(hdu.data[coli.name])) for coli in cols]

    if dim2==True:
        for coli in cols:
            #print(coli)
            newformat=coli.format
            shape = np.shape(hdu.data[coli.name])
            if len(shape) == 2 and not(coli.name in exception):
                bini =[]
                for jB in range(shape[0]):
                    binij = _rebin(hdu.data[coli.name][jB,:],binsize)
                    bini.append(binij)
                bini = np.array(bini)
                newformat = "{0}{1}".format(shape[1]//binsize,coli.format[-1])
            else:
                bini = hdu.data[coli.name]

            newcoli=fits.Column(name=coli.name,array=bini,unit=coli.unit,format=newformat)
            newcols.append(newcoli)
    else:
        for coli in cols:
            bini = _rebin(hdu.data[coli.name],binsize)
            newcoli=fits.Column(name=coli.name,array=bini,unit=coli.unit,format=coli.format)
            newcols.append(newcoli)

    newhdu=fits.BinTableHDU.from_columns(fits.ColDefs(newcols))
    newhdu.header=hdu.header
    newhdu.update()
    return newhdu


#TODO for phases bining should be done in complex plan
def binWavelength(oifits,binsize,normalizeError=True):
    

    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits

    tobin=["OI_WAVELENGTH","OI_VIS","OI_VIS2","OI_T3","OI_FLUX"]

    for i in range(1,len(data)):
        if data[i].name in tobin:
            data[i]=_rebinHdu(data[i],binsize,exception=["STA_INDEX"])
            if normalizeError:
                errname=getDataTypeError(data[i].name)
                for errnamei in errname:
                    data[i].data[errnamei]/=np.sqrt(binsize)
                

###############################################################################

def oifitsFlagWithExpression(data,arr,extver,expr,keepOldFlag = False):
    
    if not(isinstance(arr,list)):
        arr=[arr]
    
    if arr==['all']:
        arr=["OI_VIS","OI_VIS2","OI_T3","OI_FLUX"]
        
    for arri in arr:
        try:
            EFF_WAVE, EFF_BAND = oim.getWlFromOifits(data,arr=arri,
                                                     extver=extver,returnBand=True)
            nwl = np.size(EFF_WAVE)
            LENGTH, PA = oim.getBaselineLengthAndPA(data,arr=arri,extver=extver)
            nB = np.size(LENGTH)
            
            EFF_WAVE = np.tile(EFF_WAVE[None,:],(nB,1))
            EFF_BAND = np.tile(EFF_BAND[None,:],(nB,1))
            LENGTH   = np.tile(LENGTH[:,None],(1,nwl))
            PA       = np.tile(PA[:,None],(1,nwl))    
            
        
            for colname in data[arri].columns:
                coldata = data[arri].data[colname.name]
                s = coldata.shape
                
                if len(s) == 1 and s[0] == nB:
                    coldata = np.tile(LENGTH[:,None],(1,nwl))
                    
                exec(f"{colname.name}=coldata")
        except:
            pass
    
    f =eval(expr)
    for arri in arr:
        if keepOldFlag ==  True:
            data[arri].data["FLAG"] = np.logical_or(f, data[arri].data["FLAG"])
        else:
            data[arri].data["FLAG"] = f
        
    return f

###############################################################################

def computeDifferentialError(oifits,ranges=[[0,5]],excludeRange=False,
                             rangeType="index",dataType="VISPHI",extver= [None]):
    
    if rangeType == "index":
        dtype = 'int64'
    else:
        dtype = 'float64'
        
    ranges=np.array(ranges,dtype=dtype)
    
    if ranges.ndim == 1:
        ranges= ranges.reshape(1,ranges.size)
    
    if not(isinstance(dataType,list)):
        dataType=[dataType]
    
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits
        
    extnames=np.unique([getDataArrname(dti) for dti in dataType])
    
    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header["EXTVER"] in extver or extver == [None]:
                
                wl = getWlFromOifits(oifits, arr=datai.name,
                                     extver=datai.header["EXTVER"])
                nwl=wl.size
                nB=np.size(datai.data['TARGET_ID'])
                
               
                idx=np.array([],dtype="int64")              
                for ri in ranges:                    
                    if rangeType == "index":
                        idx = np.append(idx,np.arange(ri[0],ri[1]))
                    else:
                        idxi = np.where( (wl>= ri[0]) & (wl <= ri[1]) )[0]
                        idx = np.append(idx,idxi)
                        
                    
                    if excludeRange == True:
                        idx = np.delete(np.arange(nwl),idx)
                
                for dataTypei in getDataType(datai.name ):
                    if dataTypei in dataType:
                        #print(dataTypei)
                        dataTypeiErr = _oimDataTypeErr[_oimDataType.index(dataTypei)]
                        if getDataTypeIsAnalysisComplex(dataTypei) == True:
                            #print("complex analysis")
                            for iB in range(nB):
                                #print( np.mean(datai.data[dataTypeiErr][iB,:]))
                                err = circstd(datai.data[dataTypei][iB,idx],
                                                low=-180,high=180)
                                datai.data[dataTypeiErr][iB,:] = err
                                #print(err)

                        else:
                            for iB in range(nB):
                                #print( np.mean(datai.data[dataTypeiErr][iB,:]))
                                err = np.std(datai.data[dataTypei][iB,idx])
                                datai.data[dataTypeiErr][iB,:] = err
                                #print(err)

                    
###############################################################################   
   
def setMinimumError(oifits,dataTypes,values,extver=None):
                 
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits

    if not(isinstance(dataTypes,list)):
        dataTypes=[dataTypes]
                
    if not(isinstance(values,list)):
        values=[values] 
    
    if not(isinstance(extver,list)):
        extver=[extver]   
    
    extnames=np.unique([getDataArrname(dti) for dti in dataTypes])
    
    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header["EXTVER"] in extver or extver == [None]:
                
                for dataTypei in getDataType(datai.name ):
                    if dataTypei in dataTypes:
                        dataTypeiErr = _oimDataTypeErr[_oimDataType.index(dataTypei)]
                        vali = values[dataTypes.index(dataTypei)]
                        if getDataTypeIsAnalysisComplex(dataTypei) == True:
                            #print(f"Setting a minimum of {vali} deg to {dataTypeiErr} ")
                            mask = (datai.data[dataTypeiErr]<vali) \
                                    .astype(datai.data[dataTypeiErr].dtype)

                            datai.data[dataTypeiErr] = datai.data[dataTypeiErr] * \
                                         (1 - mask) +  mask * vali
                            #print(datai.data[dataTypeiErr])                        
                        else:
                            #print(f"Setting a minimum of {vali} to {dataTypeiErr} ")
                            vali = vali / 100
                            mask = ((datai.data[dataTypeiErr]/datai.data[dataTypei]) \
                                     <vali).astype(datai.data[dataTypeiErr].dtype)
                            #print(mask)
                            datai.data[dataTypeiErr] = datai.data[dataTypeiErr] * \
                                (1 - mask) + mask * vali * \
                                    datai.data[dataTypei] 
                            #print(datai.data[dataTypeiErr])     
                                
    
                        

