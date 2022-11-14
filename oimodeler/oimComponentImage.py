# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:19:06 2022

@author: Ame
"""

import numpy as np
from astropy import units as units
from scipy import interpolate,integrate
from scipy.special import j0
import oimodeler as oim
from oimodeler import oimParam,oimComponent,_standardParameters


class oimComponentImage(oimComponent):    
    """
    Base class for components define in 2D : x,y (regular grid) in the image plan.
    """
    
    elliptic=False
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        
        self._wl=None
        self._t=None
        self._pixSize=0 #in rad
        
        self._allowExternalRotation=True
        self.normalizeImage=True
               
        self.params["dim"]=oimParam(name="dim",value=256,
                         description="Image dimension",unit=1,free=False)
        
        
        self.params["dim"]=oimParam(**_standardParameters["pixSize"])
        
        self.params["pa"]=oimParam(**_standardParameters["pa"])
        
        
       
        #Add ellipticity 
        if self.elliptic==True:
            self.params["elong"]=oimParam(**_standardParameters["elong"])
            
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
        
        im0=self.getInternalImage(wl,t)
        im=self._padImage(im0)
        pix=self._pixSize
        
        tr=self._ftTranslateFactor(ucoord,vcoord,wl,t)*self.params["f"](wl,t)
        
        
        if self._allowExternalRotation==True:   
            pa_rad=(self.params["pa"](wl,t))* \
                        self.params["pa"].unit.to(units.rad)      
            co=np.cos(pa_rad)
            si=np.sin(pa_rad)
            fxp=ucoord*co-vcoord*si
            fyp=ucoord*si+vcoord*co
            vcoord=fyp
            if self.elliptic==True:
                ucoord=fxp/self.params["elong"](wl,t)
            else:
                ucoord=fxp
                
        if self._wl is None:
            wl0=np.sort(np.unique(wl))
        else:
            wl0= self._wl
            
        if self._t is None:
            t0=np.sort(np.unique(t))
        else:
            t0= self._t         
                
        if self.FTBackend.check(self.FTBackendData,im,pix,wl0,
                                     t0,ucoord,vcoord,wl,t)==False:

                   self.FTBackendData=self.FTBackend.prepare(im,pix,wl0,
                                                   t0,ucoord,vcoord,wl,t)
                   
        vc=self.FTBackend.compute(self.FTBackendData,im,pix,wl0,
                                                  t0,ucoord,vcoord,wl,t)
        
        return vc*tr*self.params["f"](wl,t)
    
    
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
            pa_rad=(self.params["pa"](wl_arr,t_arr))* \
                               self.params["pa"].unit.to(units.rad)
                               
            xp=x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp=x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            y_arr=yp
            
            if self.elliptic==True:
                x_arr=xp*self.params["elong"](wl_arr,t_arr)
            else:
                x_arr=xp
            

        im0=self._internalImage()
        
        if im0 is None:
            im = self._imageFunction(x_arr,y_arr,wl_arr,t_arr)
        else:      
            im0=np.swapaxes(im0,-2,-1)
            grid=self._getInternalGrid()
            coord=np.transpose(np.array([t_arr,wl_arr,x_arr,y_arr]))

            im=interpolate.interpn(grid,im0,coord,bounds_error=False,fill_value=0)
            f0=np.sum(im0)
            f=np.sum(im)
            im=im/f*f0
            
        im=im.reshape(dims)
        
        if self.normalizeImage==True:
            #TODO no loop for normalization
            tot=np.sum(im,axis=(2,3))
            for it,ti in enumerate(t):
                for iwl,wli in enumerate(wl):
                    if tot[it,iwl]!=0:  
                        im[it,iwl,:,:] = im[it,iwl,:,:]  \
                                  / tot[it,iwl]*self.params["f"](wli,ti)    
            
            
        
        return im
    
    def getInternalImage(self,wl,t):
        res=self._internalImage()
        
        if res is None:
            t_arr,wl_arr,x_arr,y_arr=self._getInternalGrid(simple=False,wl=wl,t=t)
            res = self._imageFunction(x_arr,y_arr,wl_arr,t_arr)
            
        
        for it in range(res.shape[0]):
            for iwl in range(res.shape[1]):
                res[it,iwl,:,:]=res[it,iwl,:,:]/np.sum(res[it,iwl,:,:])
        
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
        
        return np.pad(im, ((0,0),(0,0),(padx, padx),(pady, pady)),
                      'constant',constant_values=0)
      
    
    def _imageFunction(self,xx,yy,wl,t):
        image=xx*0+1
        return image    
    
   
    def _getInternalGrid(self,simple=True,flatten=False,wl=None,t=None):
        
        
        if self._wl is None:
            wl0=np.sort(np.unique(wl))
        else:
            wl0= self._wl
            
        if self._t is None:
            t0=np.sort(np.unique(t))
        else:
            t0= self._t            
        
        pix=self._pixSize*units.rad.to(units.mas)
        v=np.linspace(-0.5,0.5,self.params["dim"].value)        
        xy=v*pix*self.params["dim"].value
        
        if simple==True:
            return t0,wl0,xy,xy
        
        else:
            t=np.array(t0).flatten()
            nt=t.size
            wl=np.array(wl0).flatten()
            nwl=wl.size
            
            xx,yy=np.meshgrid(xy,xy)
            x_arr=np.tile(xx[None,None,:,:], (nt,nwl, 1, 1))
            y_arr=np.tile(yy[None,None,:,:], (nt,nwl, 1, 1))
            wl_arr=np.tile(wl[None,:,None,None], (nt,1, 
                        self.params["dim"].value, self.params["dim"].value))
            t_arr=np.tile(t[:,None,None,None],(1,nwl,
                        self.params["dim"].value, self.params["dim"].value))
            
            if flatten==True:
                return t_arr.flatten(),wl_arr.flatten(), \
                        x_arr.flatten(),y_arr.flatten()
            else:
                return t_arr,wl_arr,x_arr,y_arr
            
       
###############################################################################


class oimComponentRadialProfile(oimComponent):    
    """
    Base class for components define by their radial profile
    """
    
    elliptic=False
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.normalizeImage=True
        self.hankel=self.fht     

        self._t=None
        self._wl=None
        self._r=None
        
        #Add ellipticity 
        if self.elliptic==True:
            self.params["pa"]=oimParam(**_standardParameters["pa"])
            self.params["elong"]=oimParam(**_standardParameters["elong"])

        self._eval(**kwargs)      
                

    @staticmethod
    def fht(Ir,r,wlin,tin,sfreq,wl,t):
        
        pad=oim.oimOptions['FTpaddingFactor']
        nr=r.size
        ntin=tin.size
        nwlin=wlin.size
        
        fov=r[-1]-r[0]
        dsfreq0=1/(fov*pad)
        sfreq0=np.linspace(0,pad*nr-1,pad*nr)*dsfreq0      

        
        r2D=np.tile(r[None,None,:,None], (ntin,nwlin,1,pad*nr))
        Ir2D=np.tile(Ir[:,:,:,None], (1,1,1,pad*nr))
        sf2D=np.tile(sfreq0[None,None,None,:], ( ntin,nwlin,nr,1)) 
      

        res0=integrate.trapezoid(2.*np.pi*r2D*Ir2D*j0(2.*np.pi*r2D*sf2D), r2D,axis=2) 

      
        for it in range(ntin):
            for iwl in range(nwlin):
                res0[it,iwl,:]/= integrate.trapezoid(2.*np.pi*r*Ir[it,iwl,:], r)
                
                #res0[it,iwl,:]/=np.sum(r*Ir[it,iwl,:]*dr)
        
        grid=(tin,wlin,sfreq0)  
        
        coord=np.transpose([t,wl,sfreq])  
        real=interpolate.interpn(grid,np.real(res0),coord,bounds_error=False,fill_value=None)
        imag=interpolate.interpn(grid,np.imag(res0),coord,bounds_error=False,fill_value=None)
        vc=real+imag*1j
        
        #print("interp {:.2f}ms".format((time.time() - s0)*1000))    
         
        
        return vc
        
    
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
        
        if wl is None:
            wl=ucoord*0
        if t is None:
            t=ucoord*0
        
        if self.elliptic==True:   
            pa_rad=(self.params["pa"](wl,t))* \
                        self.params["pa"].unit.to(units.rad)      
            co=np.cos(pa_rad)
            si=np.sin(pa_rad)
            fxp=ucoord*co-vcoord*si
            fyp=ucoord*si+vcoord*co
            vcoord=fyp
            ucoord=fxp/self.params["elong"](wl,t)      
        
        spf=np.sqrt(ucoord**2+vcoord**2)
        
               
        if self._wl is None:
            wl0=np.sort(np.unique(wl))
        else:
            wl0= self._wl
            
        if self._t is None:
            t0=np.sort(np.unique(t))
        else:
            t0= self._t    
        
        
        Ir=self.getInternalRadialProfile(wl0,t0)

        vc=self.hankel(Ir, self._r*units.mas.to(units.rad),wl0,t0,spf,wl,t)
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    
    
    def _radialProfileFunction(self,r=None,wl=None,t=None):
        return 0
    
    def getInternalRadialProfile(self,wl,t):
        res=self._internalRadialProfile()
        
        
        if res is None:
            t_arr,wl_arr,r_arr=self._getInternalGrid(simple=False,wl=wl,t=t)
            res = self._radialProfileFunction(r_arr,wl_arr,t_arr)
            
        
        for it in range(res.shape[0]):
            for iwl in range(res.shape[1]):
                res[it,iwl,:]=res[it,iwl,:]/np.sum(res[it,iwl,:])
        
        return res
    
    def _getInternalGrid(self,simple=True,flatten=False,wl=None,t=None):
        
        if self._wl is None:
            wl0=np.sort(np.unique(wl))
        else:
            wl0= self._wl
            
        if self._t is None:
            t0=np.sort(np.unique(t))
        else:
            t0= self._t  
     
        if self._r is None:
            pix=self._pixSize*units.rad.to(units.mas)
            r=np.linspace(0,(self.params["dim"].value-1)*pix,self.params["dim"].value)   
        else:
            r=self._r
            
        if simple==True:
            return r,wl,t
        else:
            nt=np.array(t0).flatten().size
            nwl=np.array(wl0).flatten().size
            nr=r.flatten().size

            r_arr=np.tile(r[None,None,:], (nt,nwl, 1))
            wl_arr=np.tile(wl[None,:,None], (nt,1, nr))
            t_arr=np.tile(t[:,None,None],  (1,nwl, nr))
            
            if flatten==True:
                return t_arr.flatten(),wl_arr.flatten(), r_arr.flatten()
            else:
                return t_arr,wl_arr,r_arr
    
    def _internalRadialProfile(self):
        return None
    
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
        
        if self.elliptic==True :
            pa_rad=(self.params["pa"](wl_arr,t_arr))* \
                               self.params["pa"].unit.to(units.rad)
                               
            xp=x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp=x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            y_arr=yp
            x_arr=xp*self.params["elong"](wl_arr,t_arr)

        r_arr=np.sqrt(x_arr**2+y_arr**2)
        im=self._radialProfileFunction(r_arr,wl_arr,t_arr)
       
        im=im.reshape(dims)
        
        if self.normalizeImage==True:
            #TODO no loop for normalization
            tot=np.sum(im,axis=(2,3))
            for it,ti in enumerate(t):
                for iwl,wli in enumerate(wl):
                    if tot[it,iwl]!=0:  
                        im[it,iwl,:,:] = im[it,iwl,:,:]  \
                                  / tot[it,iwl]*self.params["f"](wli,ti)    

        return im

