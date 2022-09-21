# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 08:37:23 2022

@author: Ame
"""
import numpy as np
from astropy.io import fits
import astropy.units as units
from astroquery.simbad import Simbad
from  astropy.coordinates import Angle

###############################################################################

def getBaselineName(oifits,hduname="OI_VIS2",length=False,angle=False):
    """
    Return the baseline names, i.e. telescopes names 
    separated by minus sign, in an extension of a oifits file. 
    By default it is reading the OI_VIS extension
    

    Parameters
    ----------
    oifits : astropy.io.fits.hdu.hdulist.HDUList
        An oifits file structure already opened with astropy.io.fits
    hduname : str, optional
        The fits extension name. The default is "OI_VIS2".
    length : bool, optional
        Add baseline length to the returned result. The default is False.
    angle : bool, optional
        Add baseline position angle ((in deg)à=)to the returned result. 
        The default is False.

    Returns
    -------
    name :  python list of str
        The array containing the baseline names (or triplet) and optionally 
        the baseline length and orientation.
    """
    stanames=oifits['OI_ARRAY'].data['STA_NAME']
    staindexes=oifits['OI_ARRAY'].data['STA_INDEX']

    staidx=oifits[hduname].data['STA_INDEX']

    shape=np.shape(staidx)
    name=[]
    if length or angle and hduname!="OI_T3":
        u=oifits[hduname].data['UCOORD']
        v=oifits[hduname].data['VCOORD']
        B=np.sqrt(u**2+v**2)
        PA=np.rad2deg(np.arctan2(u,v))
    for i in range(shape[0]):
        namei=""
        for j in range(shape[1]):
            namei+=stanames[np.where(staindexes==staidx[i,j])[0]][0]
            if j<shape[1]-1:
                namei+="-"
        if hduname!="OI_T3":
            if length:
                namei+=" {0:.0f}m".format(B[i])
            if angle:
                namei+=" {0:.0f}$^o$".format(PA[i])
        name.append(namei)

    return name

###############################################################################

def getBaselineLengthAndPA(oifits,arr="OI_VIS2"):
    """
    
    Return a tuple (B, PA) of the baseline lengths and orientation
    (position angles) from a fits extension within an opened oifits file.
    By default it is reading the OI_VIS extension.
    

    Parameters
    ----------
    oifits : astropy.io.fits.hdu.hdulist.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr : str, optional
        The fits extension name. The default is "OI_VIS2".

    Returns
    -------
    B : numpy.ndarray
        the array containing the baselines length.
    PA : numpy.ndarray
        the array containing the baselines orientation (in deg).

    """
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits
    u=data[arr].data["UCOORD"]
    v=data[arr].data["VCOORD"]
    B=np.sqrt(u**2+v**2)
    PA=np.rad2deg(np.arctan2(u,v))

    # TODO OI_T3
    return B,PA

###############################################################################

def getSpaFreq(oifits,arr="OI_VIS2",unit=None):
    """
    

    Parameters
    ----------
    oifits : TYPE
        DESCRIPTION.
    arr : TYPE, optional
        DESCRIPTION. The default is "OI_VIS2".
    unit : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    spaFreq : TYPE
        DESCRIPTION.

    """
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits


    if arr!="OI_T3":
        B,PA=getBaselineLengthAndPA(data,arr)

    else:
        u1=data[arr].data["U1COORD"]
        v1=data[arr].data["V1COORD"]
        u2=data[arr].data["U2COORD"]
        v2=data[arr].data["V2COORD"]
        u3=u1+u2
        v3=v1+v2
        B1=np.sqrt(u1**2+v1**2)
        B2=np.sqrt(u2**2+v2**2)
        B3=np.sqrt(u3**2+v3**2)
        Bs=[B1,B2,B3]
        B=np.max(Bs,axis=0)

    lam=data["OI_WAVELENGTH"].data["EFF_WAVE"]
    nlam=np.size(lam)
    nB=np.size(B)
    if unit=="cycles/mas":
        mult=units.mas.to(units.rad)
    elif unit=="cycles/arcsec":
        mult=units.arcsec.to(units.rad)
    elif unit == "Mlam":
        mult = 1/(1e6)
    else:
        mult=1

    spaFreq=np.ndarray([nB,nlam])
    for iB in range(nB):
        spaFreq[iB,:]=B[iB]/lam*mult

    return spaFreq

###############################################################################

def createOiArray(arrname,arrx,arry,arrz,sta_name,tel_name,diameter,staxyz):
    """
    

    Parameters
    ----------
    arrname : TYPE
        DESCRIPTION.
    arrx : TYPE
        DESCRIPTION.
    arry : TYPE
        DESCRIPTION.
    arrz : TYPE
        DESCRIPTION.
    sta_name : TYPE
        DESCRIPTION.
    tel_name : TYPE
        DESCRIPTION.
    diameter : TYPE
        DESCRIPTION.
    staxyz : TYPE
        DESCRIPTION.

    Returns
    -------
    arr : TYPE
        DESCRIPTION.

    """
    nstation=np.size(sta_name)
    tel_name=fits.Column(name="TEL_NAME",format="A16",array=np.array(tel_name))
    sta_name=fits.Column(name="STA_NAME",format="A16",array=np.array(sta_name))
    sta_index=fits.Column(name="STA_INDEX",format="I2",array=np.arange(1,nstation+1))
    diameter=fits.Column(name="DIAMETER",format="E",array=np.array(diameter),unit='m')
    staxyz=fits.Column(name="STAXYZ",format="3D",array=np.array(staxyz),unit='m')

    cols=[tel_name,sta_name,sta_index,diameter,staxyz]
    arr=fits.BinTableHDU.from_columns(cols)

    arr.header['EXTVER']=(1,'ID number of this OI_ARRAY')
    arr.header['ARRAYX']=float(arrx)#,'[m] Array center X coordinate')
    arr.header['ARRAYY']=float(arry)#,'[m] Array center Y coordinate')
    arr.header['ARRAYZ']=float(arrz)#,'[m] Array center Z coordinate')
    arr.header['FRAME']=('GEOCENTRIC','Coordinate frame')
    arr.header['EXTNAME']='OI_ARRAY'
    arr.header['OI_REVN']=(1,'Revision number of the table definition')
    arr.header['ARRNAME']=(arrname,'Array name')

    return arr
###############################################################################


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
    customSimbad.add_votable_fields('plx','plx_error','propermotions','sptype','velocity')
    if type(names)==type(""):
        names=[names]
    data=customSimbad.query_objects(names)

    ntargets=len(names)
    rad=Angle(data['RA'],unit="hourangle").deg
    dec=Angle(data['DEC'],unit="deg").deg
    ra_err=(data['COO_ERR_MAJA'].data*units.mas).to_value(unit='deg')
    dec_err=(data['COO_ERR_MINA'].data*units.mas).to_value(unit='deg')
    pmra=(data['PMRA'].data*units.mas).to_value(unit='deg')
    pmdec=(data['PMDEC'].data*units.mas).to_value(unit='deg')
    pmra_err=(data['PM_ERR_MAJA'].data*units.mas).to_value(unit='deg')
    pmdec_err=(data['PM_ERR_MINA'].data*units.mas).to_value(unit='deg')
    plx_value=(data['PLX_VALUE'].data*units.mas).to_value(unit='deg')
    plx_error=(data['PLX_ERROR'].data*units.mas).to_value(unit='deg')

    target_id=fits.Column(name="TARGET_ID",format="I",array=np.arange(1,ntargets+1))
    target=fits.Column(name="TARGET",format="16A",array=names)
    raep0=fits.Column(name="RAEP0",format="D",array=rad,unit="deg")
    decep0=fits.Column(name="DECEP0",format="D",array=dec,unit="deg")
    equinox=fits.Column(name="EQUINOX",format="E",array=np.repeat(2000,ntargets),unit="yr")
    ra_err=fits.Column(name="RA_ERR",format="D",array=ra_err,unit="deg")
    dec_err=fits.Column(name="DEC_ERR",format="D",array=dec_err,unit="deg")
    sysvel=fits.Column(name="SYSVEL",format="D",array=np.zeros([ntargets]),unit="m/s")   #TODO
    veltyp=fits.Column(name="VELTYP",format="8A",array=np.repeat("UNKNOWN",ntargets))  #TODO
    veldef=fits.Column(name="VELDEF",format="8A",array=np.repeat("OPTICAL",ntargets)) #TODO
    pmra=fits.Column(name="PMRA",format="D",array=pmra,unit="deg/yr")
    pmdec=fits.Column(name="PMDEC",format="D",array=pmdec,unit="deg/yr")
    pmra_err=fits.Column(name="PMRA_ERR",format="D",array=pmra_err,unit="deg/yr")
    pmdec_err=fits.Column(name="PMDEC_ERR",format="D",array=pmdec_err,unit="deg/yr")
    parallax=fits.Column(name="PARALLAX",format="E",array=plx_value,unit="deg")
    para_err=fits.Column(name="PARA_ERR",format="E",array=plx_error,unit="deg")
    spectyp=fits.Column(name="SPECTYP",format="16A",array=data['SP_TYPE'])


    cols=[target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,
          veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp]

    tar=fits.BinTableHDU.from_columns(cols)

    tar.header['OI_REVN']=(1,'Revision number of the table definition')
    tar.header['EXTVER']=(1,'ID number of this OI_TARGET')
    tar.header['EXTNAME']='OI_TARGET'
    return tar

###############################################################################

def createOiWavelength(insname,eff_wave,eff_band):
    """
    

    Parameters
    ----------
    insname : TYPE
        DESCRIPTION.
    eff_wave : TYPE
        DESCRIPTION.
    eff_band : TYPE
        DESCRIPTION.

    Returns
    -------
    wave : TYPE
        DESCRIPTION.

    """
    eff_wave=fits.Column(name="EFF_WAVE",format="E",array=np.array(eff_wave),unit='m')
    eff_band=fits.Column(name="EFF_BAND",format="E",array=np.array(eff_band),unit='m')

    cols=[eff_wave,eff_band]
    wave=fits.BinTableHDU.from_columns(cols)
    wave.header['EXTNAME']='OI_WAVELENGTH'
    wave.header['EXTVER']=(1,'ID number of this OI_WAVELENGTH')
    wave.header['OI_REVN']=(1,'Revision number of the table definition')
    wave.header['INSNAME']=(insname,'Identifies corresponding OI_WAVELENGTH')

    return wave

###############################################################################

def createOiVis2(arrname,insname,target_id,time,mjd,int_time,vis2data,vis2err,
                 ucoord,vcoord,sta_index,flag,dateobs):
    """
    

    Parameters
    ----------
    arrname : TYPE
        DESCRIPTION.
    insname : TYPE
        DESCRIPTION.
    target_id : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.
    mjd : TYPE
        DESCRIPTION.
    int_time : TYPE
        DESCRIPTION.
    vis2data : TYPE
        DESCRIPTION.
    vis2err : TYPE
        DESCRIPTION.
    ucoord : TYPE
        DESCRIPTION.
    vcoord : TYPE
        DESCRIPTION.
    sta_index : TYPE
        DESCRIPTION.
    flag : TYPE
        DESCRIPTION.
    dateobs : TYPE
        DESCRIPTION.

    Returns
    -------
    oivis2 : TYPE
        DESCRIPTION.

    """

    nb=len(target_id)
    
    v2shape=np.shape(vis2data)
    v2dim=len(v2shape)
    
    if (v2dim==2):
        nlam=v2dim[1]
    elif (v2dim==1):
        if v2shape[0]==nb:
            nlam=1
        else:
            nlam=v2shape[0]

    target_id=fits.Column(name="TARGET_ID",format="I",array=np.array(target_id))
    time=fits.Column(name="TIME",format="D",array=np.array(time),unit="sec")
    mjd=fits.Column(name="MJD",format="D",array=np.array(mjd),unit="day")
    int_time=fits.Column(name="INT_TIME",format="D",array=np.array(int_time),unit="sec")
    vis2data=fits.Column(name="VIS2DATA",format="{0}D".format(nlam),array=np.array(vis2data))
    vis2err=fits.Column(name="VIS2ERR",format="{0}D".format(nlam),array=np.array(vis2err))
    ucoord=fits.Column(name="UCOORD",format="1D",array=np.array(ucoord),unit="m")
    vcoord=fits.Column(name="VCOORD",format="1D",array=np.array(vcoord),unit="m")
    sta_index=fits.Column(name="STA_INDEX",format="2I",array=np.array(sta_index))
    flag=fits.Column(name="FLAG",format="L",array=np.array(flag))
    
    cols=[target_id,time,mjd,int_time,vis2data,vis2err,ucoord,vcoord,sta_index,flag]
      
    oivis2=fits.BinTableHDU.from_columns(cols)
    oivis2.header['EXTNAME']='OI_VIS2'
    oivis2.header['EXTVER']=(1,'ID number of this OI_VIS2')
    oivis2.header['OI_REVN']=(1,'Revision number of the table definition')
    oivis2.header['INSNAME']=(insname,'Identifies corresponding OI_WAVELENGTH')
    oivis2.header['DATE-OBS']=dateobs
    oivis2.header['ARRNAME']=arrname

    return oivis2

