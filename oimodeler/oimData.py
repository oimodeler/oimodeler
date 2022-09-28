# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""

import numpy as np
from astropy.io import fits
import os
from enum import IntFlag
import oimodeler as oim

_oimDataType=["VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI","FLUXDATA"]
_oimDataTypeErr=["VIS2ERR","VISAMPERR","VISPHIERR","T3AMPERR","T3PHIERR","FLUXERR"]
_oimDataTypeArr=["OI_VIS2","OI_VIS","OI_VIS","OI_T3","OI_T3","OI_FLUX"]


###############################################################################

###############################################################################    
    
def oimDataGetWl(data,arr,dwl=True):
        insname=arr.header['INSNAME']
        oiWlArr=[arri for arri in data if (arri.name=="OI_WAVELENGTH" 
                                           and arri.header['INSNAME']==insname)][0]
        if dwl==False:
            return oiWlArr.data["EFF_WAVE"]
        else:
            return oiWlArr.data["EFF_WAVE"],oiWlArr.data["EFF_BAND"]

###############################################################################

class oimDataType(IntFlag):
    NONE        = 0
    VIS2DATA    = 1
    VISAMP_ABS  = 2    
    VISAMP_DIF  = 4  
    VISAMP_COR  = 8      
    VISPHI_ABS  = 16
    VISPHI_DIF  = 32    
    T3AMP       = 64
    T3PHI       = 128
    FLUXDATA    = 256


###############################################################################

def oimGetDataValErrAndTypeFlag(arr):
    
    dtype=oimDataType(0)
    val=[]
    err=[]
    flag=[]
    if arr.name=="OI_VIS2":
        nv2=np.size(np.where(arr.data["VIS2DATA"]!=0))
        if nv2!=0:
            val.append(arr.data["VIS2DATA"]) 
            err.append(arr.data["VIS2ERR"])
            flag.append(arr.data["FLAG"])
            dtype|=oimDataType.VIS2DATA
    if arr.name=="OI_VIS":
        nvamp=np.size(np.where(arr.data["VISAMP"]!=0))
        if nvamp!=0:
            val.append(arr.data["VISAMP"]) 
            err.append(arr.data["VISAMPERR"])
            flag.append(arr.data["FLAG"])
            try:
                if arr.header["AMPTYP"].lower()=="absolute":
                    dtype|=oimDataType.VISAMP_ABS
                elif arr.header["AMPTYP"].lower()=="differential":
                    dtype|=oimDataType.VISAMP_DIF
                else:
                    dtype|=oimDataType.VISAMP_COR
            except:
                dtype|=oimDataType.VISAMP_ABS
        nvphi=np.size(np.where(arr.data["VISPHI"]!=0))
        if nvphi!=0:
            val.append(arr.data["VISPHI"]) 
            err.append(arr.data["VISPHIERR"])
            flag.append(arr.data["FLAG"])
            try:
                if arr.header["PHITYP"].lower()=="absolute":
                    dtype|=oimDataType.VISPHI_ABS
                else:
                    dtype|=oimDataType.VISPHI_DIF
            except:
                dtype|=oimDataType.VISPHI_ABS           
    if arr.name=="OI_T3": 
        t3amp=np.size(np.where(arr.data["T3AMP"]!=0))
        if t3amp!=0:
            val.append(arr.data["T3AMP"]) 
            err.append(arr.data["T3AMPERR"])
            flag.append(arr.data["FLAG"])
            dtype|=oimDataType.T3AMP
        t3phi=np.size(np.where(arr.data["T3PHI"]!=0))
        if t3phi!=0:
            val.append(arr.data["T3PHI"]) 
            err.append(arr.data["T3PHIERR"])
            flag.append(arr.data["FLAG"])
            dtype|=oimDataType.T3PHI
    if arr.name=="OI_FLUX":     
        nflx=np.size(np.where(arr.data["FLUXDATA"]!=0))
        if nflx!=0:
            val.append(arr.data["FLUXDATA"]) 
            err.append(arr.data["FLUXERR"])
            flag.append(arr.data["FLAG"])
            dtype|=oimDataType.FLUXDATA
    return dtype,val,err,flag

###############################################################################

def oimDataCheckData(arr):
    cdata=[]
    if np.size(np.where(arr.data["FLAG"]==False))!=0:
        if arr.name=="OI_VIS2":
            nv2=np.size(np.where(arr.data["VIS2DATA"]!=0))
            if nv2!=0:
                cdata.append("VIS2DATA")    
        if arr.name=="OI_VIS":
            nvamp=np.size(np.where(arr.data["VISAMP"]!=0))
            if nvamp!=0:
                cdata.append("VISAMP")   
            nvphi=np.size(np.where(arr.data["VISPHI"]!=0))
            if nvphi!=0:
                cdata.append("VISPHI")             
        if arr.name=="OI_T3": 
            t3amp=np.size(np.where(arr.data["T3AMP"]!=0))
            if t3amp!=0:
                cdata.append("T3AMP")   
            t3phi=np.size(np.where(arr.data["T3PHI"]!=0))
            if t3phi!=0:
                cdata.append("T3PHI")
        if arr.name=="OI_FLUX":     
            nflx=np.size(np.where(arr.data["FLUXDATA"]!=0))
            if nflx!=0:
                cdata.append("FLUXDATA")   
    return cdata

###############################################################################

def oimDataGetVectCoord(data,arr):
    wl,dwl=oimDataGetWl(data,arr)
    nwl=np.size(wl)    
    mjd=arr.data["TIME"]
    mjd=np.outer(mjd,np.ones(nwl)).flatten()
    
    #zero freq vector for vis normalization
    uv0=wl*0
    
    if arr.name=="OI_T3":
        nB=np.shape(arr.data["T3PHI"])[0]*3+1
        um1=arr.data["U1COORD"]
        vm1=arr.data["V1COORD"]
        um2=arr.data["U2COORD"]
        vm2=arr.data["V2COORD"]        
        um3=arr.data["U1COORD"]+arr.data["U2COORD"]
        vm3=arr.data["V1COORD"]+arr.data["V2COORD"]
        
        
        u1=np.outer(um1,1./wl).flatten()
        v1=np.outer(vm1,1./wl).flatten()
        u2=np.outer(um2,1./wl).flatten()
        v2=np.outer(vm2,1./wl).flatten()
        u3=np.outer(um3,1./wl).flatten()
        
        v3=np.outer(vm3,1./wl).flatten()
        
        u=np.concatenate((u1,u2,u3))
        v=np.concatenate((v1,v2,v3)) 
        
        mjd=np.outer(np.ones(3),mjd).flatten()
        
    elif arr.name=="OI_FLUX":
        nB=np.shape(arr.data["FLUXDATA"])[0]
        u=np.zeros(nB*nwl)
        v=np.zeros(nB*nwl)
        
    else:
        um=arr.data["UCOORD"]
        vm=arr.data["VCOORD"]
        nB=np.size(um)+1
       
        u=np.outer(um,1./wl).flatten()
        v=np.outer(vm,1./wl).flatten()
        
    if  arr.name!="OI_FLUX":
        u=np.concatenate((uv0,u))
        v=np.concatenate((uv0,v))
    
    wl=np.outer(np.ones(nB),wl).flatten()
    dwl=np.outer(np.ones(nB),dwl).flatten()
        
    return u,v,wl,dwl,mjd,nB,nwl
   
###############################################################################

class oimData(object):
    """
    A class to hold and manipulate data
    """
    def __init__(self,dataOrFilename=None,filt=None):
        self._data=[]
        self.dataInfo=[]
        self.vect_u=None
        self.vect_v=None
        self.vect_wl=None
        self.vect_dwl=None
        self.vect_mjd=None
        
        self._prepared = False
        
        self._filter=filt
        self._useFilter=False
        self._filteredData = None
        self._filteredDataReady = False
        
 
        if dataOrFilename:
           self.addData(dataOrFilename) 
        
        self.prepareData()
            
    @property
    def data(self):
        if self._useFilter == False or self._filter==None:
            return self._data
        else:
            if self._filteredDataReady == False:
                self.applyFilter()
            return self._filteredData
   

    def addData(self,dataOrFilename):
        if type(dataOrFilename)==type([]):
            for el in dataOrFilename:
                self.addData(el)
        else:        
            if type(dataOrFilename)==str:
                self._data.append(fits.open(dataOrFilename))
            else:
                self._data.append(dataOrFilename)
            self._analyzeOIFitFile(self.data[-1])
            self.prepared = False
        
        self._filteredDataReady = False
        self.prepareData()
        
    def removeData(self,dataOrIndex):
        self._prepared = False
        self._filteredDataReady = False
        self.prepareData()
        #TODO
      
        
    def setFilter(self,filt=None,useFilter=True):
        self._filter=filt
        self._filteredDataReady = False
        self.useFilter=useFilter
        
        
    def applyFilter(self):
        
        self._filteredData=[]
        for data in self._data:
            self._filteredData.append(oim.hdulistDeepCopy(data))
            
        if self._filter!=None:
            self._filter.applyFilter(self._filteredData)
            
        self._filteredDataReady =True
        self.prepareData()
        
  
    @property
    def useFilter(self):
        return self._useFilter
        
    @useFilter.setter
    def useFilter(self,val):
        self._useFilter=val
        if val==True:
            if self._filteredDataReady == False:
                self.applyFilter()
       
            
    
                
    def _analyzeOIFitFile(self,data):
        dataInfo=[]
        
        for iarr,arri in enumerate(data):
            info=None
            if arri.name in _oimDataTypeArr:
                info={'arr':arri.name, 'idx':iarr}
                cdata=oimDataCheckData(arri)
                if arri.name=="OI_VIS2":
                    nB=np.shape(arri.data["VIS2DATA"])
                if arri.name=="OI_VIS":
                    nB=np.shape(arri.data["VISAMP"])
                if arri.name=="OI_T3": 
                    nB=np.shape(arri.data["T3AMP"])
                if arri.name=="OI_FLUX":     
                    nB=np.shape(arri.data["FLUXDATA"])
                    
                info["nB"]=nB
                cdata=oimDataCheckData(arri)
                info["data"]=cdata
            if info:
                dataInfo.append(info)
        
        self.dataInfo.append(dataInfo)
        
    def prepareData(self):
        self.vect_u=np.array([])
        self.vect_v=np.array([])
        self.vect_wl=np.array([])
        self.vect_dwl=np.array([])
        self.vect_mjd=np.array([])
                   
        self.struct_u=[]
        self.struct_v=[]
        self.struct_wl=[]
        self.struct_dwl=[]
        self.struct_mjd=[]
        self.struct_nB=[]
        self.struct_nwl=[]
        self.struct_val=[]
        self.struct_err=[]  
        self.struct_flag=[]
        self.struct_arrNum=[]
        self.struct_arrType=[]
        self.struct_dataType=[]
        
        for idata,datai in enumerate(self.data):
            #print("File {}".format(idata))
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
            for iarr,arri in enumerate(datai):
                if arri.name in _oimDataTypeArr:
                    
                    #print("arr {} : type={}".format(iarr,arri.name))
                    dataTypeFlag,val,err,flag=oimGetDataValErrAndTypeFlag(arri)
                    
                    if dataTypeFlag != oimDataType.NONE:
                        u,v,wl,dwl,mjd,nB,nwl=oimDataGetVectCoord(datai,arri)
                        
                        #print(np.shape(u))
                        self.vect_u=np.concatenate((self.vect_u,u))
                        self.vect_v=np.concatenate((self.vect_v,v))
                        self.vect_wl=np.concatenate((self.vect_wl,wl))
                        self.vect_dwl=np.concatenate((self.vect_dwl,dwl))
                        self.vect_mjd=np.concatenate((self.vect_mjd,mjd))
                        
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
        
        
    def __str__(self):
       nfiles=np.size(self.data)
       txt="oimData containing {} file(s)\n".format(nfiles)
       for ifile,fi in enumerate(self.dataInfo):
           fname=os.path.basename(self.data[ifile].filename())
           txt+="{}\n".format(fname)
           for di in fi:
               shapetxt="({},{})".format(di["nB"][0],di["nB"][1])
               txt+="\t"+di["arr"]+shapetxt+"\n"
               
               
               for ddi in di["data"]:
                   txt+="\t\t"+ddi+"\n"
       return txt
            
    
###############################################################################    


