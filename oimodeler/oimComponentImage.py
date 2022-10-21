# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:19:06 2022

@author: Ame
"""

import numpy as np
from astropy import units as units
from scipy import interpolate
import oimodeler as oim
from oimodeler import oimParam,oimComponent,_standardParameters
from datetime import datetime


class oimComponentImage(oimComponent):    
    """
    Base class for components define in 2D : x,y (regular grid) in the image plan.
    """
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        
        self._wl=[0]
        self._t=[0]
        
        self._pixSize=0 #in rad
        
        self._allowExternalRotation=True
               
        self.params["dim"]=oimParam(name="dim",value=256,
                         description="Image dimension",unit=1,free=False)
        self.params["pa"]=oimParam(**_standardParameters["pa"])
        
        if 'FTBackend' in kwargs:
             self.FTBackend=kwargs['FTBackend']
        else: 
            self.FTBackend=oim.oimOptions['FTBackend']
            
        self.FTBackendData=None
        
        self._eval(**kwargs)   
    
    
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
        
        if wl is None:
            wl=ucoord*0
        if t is None:
            t=ucoord*0
        
        im0=self.getInternalImage()
        im=self._padImage(im0)
        pix=self._pixSize

        if self._allowExternalRotation==True:   
            pa_rad=(self.params["pa"](wl,t)+90)* \
                        self.params["pa"].unit.to(units.rad)      
            co=np.cos(pa_rad)
            si=np.sin(pa_rad)
            fxp=ucoord*co-vcoord*si
            fyp=ucoord*si+vcoord*co
            ucoord=fxp
            vcoord=fyp
        if self.FTBackend.check(self.FTBackendData,im,pix,self._wl,
                                                   self._t,ucoord,vcoord,wl,t)==False:

                   self.FTBackendData=self.FTBackend.prepare(im,pix,self._wl,
                                                   self._t,ucoord,vcoord,wl,t)
                   
        
        vc=self.FTBackend.compute(self.FTBackendData,im,pix,self._wl,
                                                   self._t,ucoord,vcoord,wl,t)
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    
    
    def getImage(self,dim,pixSize,wl=None,t=None):
  
        

        if wl is None:
            wl=0
        if t is None:
            t=0
         
        
        t=np.array(t).flatten()
        nt=t.size
        wl=np.array(wl).flatten()
        nwl=wl.size       
        dims=(nt,nwl,dim,dim)

        v=np.linspace(-0.5,0.5,dim)       
        vx,vy=np.meshgrid(v,v)
        
        vx_arr=np.tile(vx[None,None,:,:], (nt,nwl, 1, 1))
        vy_arr=np.tile(vy[None,None,:,:], (nt,nwl, 1, 1))
        wl_arr=np.tile(wl[None,:,None,None], (nt,1, dim, dim))
        t_arr=np.tile(t[:,None,None,None], (1,nwl, dim, dim))
        
        x_arr=(vx_arr*pixSize*dim).flatten()   
        y_arr=(vy_arr*pixSize*dim).flatten()   
        wl_arr=wl_arr.flatten()
        t_arr=t_arr.flatten()
    
        x_arr,y_arr=self._directTranslate(x_arr,y_arr,wl_arr,t_arr)
        
        
        if self._allowExternalRotation==True :
            pa_rad=(self.params["pa"](wl_arr,t_arr)+90)* \
                               self.params["pa"].unit.to(units.rad)
                               
            xp=x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp=x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            x_arr=xp
            y_arr=yp
        
        im0 = self.getInternalImage()
        grid=self._getInternalGrid()
        
        coord=np.transpose(np.array([t_arr,wl_arr,x_arr,y_arr]))

        im=interpolate.interpn(grid,im0,coord,bounds_error=False,fill_value=0)
        f0=np.sum(im0)
        f=np.sum(im)
        im=im/f*f0
        return im.reshape(dims)
    
    def getInternalImage(self):
        res=self._internalImage()
        
        if res is None:
            t_arr,wl_arr,x_arr,y_arrx_arr=self._getInternalGrid(simple=False)
            res = self._imageFunction(x_arr,y_arrx_arr,t_arr,wl_arr)
        
        return res
            
    def _internalImage(self):
        return None
    
    def _padImage(self,im):
        im0=np.sum(im,axis=(0,1))
        dimy=im0.shape[0]
        dimx=im0.shape[1]
        
        im0x=np.sum(im0,axis=1)
        im0y=np.sum(im0,axis=1)
        
        s0x=np.trim_zeros(im0x).size
        s0y=np.trim_zeros(im0y).size
        
        min_sizex=s0x*oim.oimOptions["FTpaddingFactor"]
        min_sizey=s0y*oim.oimOptions["FTpaddingFactor"]
        
        min_pow2x=2**(min_sizex - 1).bit_length()
        min_pow2y=2**(min_sizey - 1).bit_length()
        
        padx=(min_pow2x-dimx)//2
        pady=(min_pow2y-dimy)//2        
        
        return np.pad(im, ((0,0),(0,0),(padx, padx),(pady, pady)),'constant',constant_values=0)
      
    
    def _imageFunction(self,xx,yy,wl,t):
        image=xx*0+1
        return image    
    
    
    def _directTranslate(self,x,y,wl,t):
        x=x-self.params["x"](wl,t)
        y=y-self.params["y"](wl,t)
        return x,y
    
    def _getInternalGrid(self,simple=True,flatten=False):
        
        pix=self._pixSize*units.rad.to(units.mas)
        v=np.linspace(-0.5,0.5,self.params["dim"].value)        
        xy=v*pix*self.params["dim"].value
        
        if simple==True:
            return self._t,self._wl,xy,xy
        
        else:
            t=np.array(self._t).flatten()
            nt=t.size
            wl=np.array(self._wl).flatten()
            nwl=wl.size
            
            xx,yy=np.meshgrid(xy,xy)
            x_arr=np.tile(xx[None,None,:,:], (nt,nwl, 1, 1))
            y_arr=np.tile(yy[None,None,:,:], (nt,nwl, 1, 1))
            wl_arr=np.tile(wl[None,:,None,None], (nt,1, self.params["dim"].value, self.params["dim"].value))
            t_arr=np.tile(t[:,None,None,None], (1,nwl, self.params["dim"].value, self.params["dim"].value))
            
            if flatten==True:
                return t_arr.flaten(),wl_arr.flaten(), x_arr.flaten(),y_arr.flaten()
            else:
                return t_arr,wl_arr,x_arr,y_arr
            
       
                
        
    
    