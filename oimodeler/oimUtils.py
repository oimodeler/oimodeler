# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 08:37:23 2022

@author: Ame
"""
import numpy as np
from astropy.io import fits
import astropy.units as units


#--------------------------------------------baselinename-------------------------

def baselinename(oifits,hduname="OI_VIS2",length=False,angle=False):
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

#----------------------------------------------getBaselineLengthAndPA------------------------------------------

def getBaselineLengthAndPA(oifits,arr="OI_VIS2"):
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

#-------------------------------------------------getSpaFreq--------------------------------------------------

def getSpaFreq(oifits,arr="OI_VIS2",unit=None):
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


