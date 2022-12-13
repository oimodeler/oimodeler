# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""

import numpy as np
from astropy import units as units
from astropy.io import fits
from scipy.special import j0,j1,jv
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import oimodeler as oim
from oimodeler import oimParam,_standardParameters,oimInterpWl, \
                oimParamInterpWl,oimInterpTime,oimParamInterpTime, \
                oimParamLinker,oimParamInterpolator,oimInterp


###############################################################################
"""
Useful definitions.
Maybeto be  put in a separate file later
"""
mas2rad=units.mas.to(units.rad)
I=complex(0,1)
###############################################################################

def getFourierComponents():
    """
    A function to get the list of all available components deriving from the 
    oimComponentFourier class.

    Returns
    -------
    res : list
        list of all available components deriving from the oimComponentFourier class.
    """
    fnames=dir(oim)
    res=[]
    for f in fnames:
        try:
            if issubclass(oim.__dict__[f], oim.oimComponentFourier):
                res.append(f)
        except:
            pass
    return res

###############################################################################

class oimComponent(object):
    """
    The OImComponent class is the parent abstract class for all types of 
    components that can be added to a OImModel. It has a similar interface than
    the oimModel and allow to compute images (or image cubes fore wavelength 
    dependent models or time dependent models) and complex coherentFluxes for a
    vector of u,v,wl, and t coordinates. 
    
    variables:
        name: the name of the component
        shortname: a short name for the component
        description: a detail description of the component
        params: the dictionary of the component parameters
        
    """  
    
  
    name = "Generic component"
    shortname = "Gen comp"
    description = "This is the class from which all components derived"
    
    def __init__(self,**kwargs):
        """
        Create and initiliaze a new instance of the oimComponent class
        All components have at least three parameters the position 
        x and y and their flux f.

        Returns
        -------
        None.

        """
        self.params={}
        self.params["x"]=oimParam(**_standardParameters["x"])
        self.params["y"]=oimParam(**_standardParameters["y"])
        self.params["f"]=oimParam(**_standardParameters["f"])
       
        self._eval(**kwargs)
       
    def _eval(self,**kwargs):
        for key, value in kwargs.items():
            if key in self.params.keys():
                if isinstance(value,oimInterpWl):
                    if not(isinstance(self.params[key],oimParamInterpWl)):
                       self.params[key]=oimParamInterpWl(self.params[key],value)
                elif isinstance(value,oimInterpTime):
                    if not(isinstance(self.params[key],oimParamInterpTime)):
                        self.params[key]=oimParamInterpTime(self.params[key],value)  
                elif isinstance(value,oimInterp):

                    if not(isinstance(self.params[key],oimParamInterpolator)):
                        self.params[key]=value.type(self.params[key],**value.kwargs)
                else:
                    self.params[key].value=value
   
                            
    def getComplexCoherentFlux(self,u,v,wl=None,t=None):
        """
        Compute and return the complex coherent flux for an array of u,v 
        (and optionally wavelength and time ) coordinates.

        Parameters
        ----------
        u : list or numpy array
            spatial coordinate u (in cycles/rad) 
        v : list or numpy array
            spatial coordinate vu (in cycles/rad) .
        wl : list or numpy array, optional
            wavelength(s) in meter. The default is None.
        t :  list or numpy array, optional
            time in s (mjd). The default is None.

        Returns
        -------
        A numpy array of  the same size as u & v
            The complex coherent flux.

        """
        return np.array(u)*0;
    
    def getImage(self,dim,pixSize,wl=None,t=None):
        """
        Compute and return an image or and image cube (if wavelength and time 
        are given). The returned image as the x,y dimension dim in pixel with
        an angular pixel size pixSize in rad. Image is returned as a numpy 
        array unless the keyword fits is set to True. In that case the image is
        returned as an astropy.io.fits hdu.

        Parameters
        ----------
        dim : integer
            image x & y dimension in pixels.
        pixSize : float
            pixel angular size.in rad
        wl : integer, list or numpy array, optional
             wavelength(s) in meter. The default is None.
        t :  integer, list or numpy array, optional
            time in s (mjd). The default is None.
        fits : bool, optional
            if True returns result as a fits hdu. The default is False.

        Returns
        -------
        a numpy 2D array (or 3 or 4D array if wl, t or both are given) or an
        astropy.io.fits hdu. image hdu if fits=True.
            The image of the component with given size in pixels and rad.

        """          
        return np.zeros(dim,dim);    

        
    def _ftTranslateFactor(self,ucoord,vcoord,wl,t): 
        x=self.params["x"](wl,t)*self.params["x"].unit.to(units.rad)
        y=self.params["y"](wl,t)*self.params["y"].unit.to(units.rad)
        return np.exp(-2*I*np.pi*(ucoord*x+vcoord*y))
    
    def _directTranslate(self,x,y,wl,t):
        x=x-self.params["x"](wl,t)
        y=y-self.params["y"](wl,t)
        return x,y
    
    
    def __str__(self):
        txt=self.name
        for name,param in self.params.items():  
            txt+=" {0}={1:.2f}".format(param.name,param.value)
        return txt
    
#


class oimComponentFourier(oimComponent):
    """
    Class  for all component analytically defined in the Fourier plan.
    Inherit from the oimComponent.
    Implements translation in direct and Fourier space, getImage from the
    Fourier definition of the object, ellipticity (i.e. flatening)
    Children classes should only implement the _visFunction and _imageFunction
    functions. 
    """
    elliptic=False
    def __init__(self,**kwargs):      

        super().__init__(**kwargs)
        
        #Add ellipticity if either elong or pa is specified in kwargs
        if ("elong" in kwargs) or ("pa" in kwargs) or self.elliptic==True:
            self.params["elong"]=oimParam(**_standardParameters["elong"])
            self.params["pa"]=oimParam(**_standardParameters["pa"])
            self.elliptic=True
        
        self._eval(**kwargs)

        
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
        
        if self.elliptic==True:   
            pa_rad=(self.params["pa"](wl,t))* \
                        self.params["pa"].unit.to(units.rad)      
            co=np.cos(pa_rad)
            si=np.sin(pa_rad)
            fxp=ucoord*co-vcoord*si
            fyp=ucoord*si+vcoord*co
            rho=np.sqrt(fxp**2/self.params["elong"](wl,t)**2+fyp**2) 
        else:
            fxp=ucoord
            fyp=vcoord
            rho=np.sqrt(fxp**2+fyp**2)            
               
        vc=self._visFunction(fxp,fyp,rho,wl,t)
              
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    
          
    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        return ucoord*0
    
    def getImage(self,dim,pixSize,wl=None,t=None):  

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
        
        if self.elliptic:
       
            pa_rad=(self.params["pa"](wl_arr,t_arr))* \
                               self.params["pa"].unit.to(units.rad)
                               
            xp=x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp=x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            
            x_arr=xp*self.params["elong"](wl_arr,t_arr)
            y_arr=yp
            
        image = self._imageFunction(x_arr.reshape(dims),y_arr.reshape(dims),
                                    wl_arr.reshape(dims),t_arr.reshape(dims))
        
        tot=np.sum(image,axis=(2,3))
        
        #TODO no loop for normalization
        for it,ti in enumerate(t):
            for iwl,wli in enumerate(wl):
                if tot[it,iwl]!=0:  
                    image[it,iwl,:,:] = image[it,iwl,:,:]  \
                              / tot[it,iwl]*self.params["f"](wli,ti)    

        return image    
    
    
    def _imageFunction(self,xx,yy,wl,t):
        image=xx*0+1
        return image
    
    
      
###############################################################################
class oimPt(oimComponentFourier):
    name="Point source"
    shortname = "Pt"   
    def __init__(self,**kwargs):        
         super().__init__(**kwargs) 
         self._eval(**kwargs)

    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        return 1
    
    def _imageFunction(self,xx,yy,wl,t):
        image=xx*0
        val=np.abs(xx)+np.abs(yy)
        idx=np.unravel_index(np.argmin(val),np.shape(val))
        image[idx]=1
        return image
    
###############################################################################
class oimBackground(oimComponentFourier):
    name="Background"
    shortname = "Bckg"     
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)

    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        vc=rho*0
        idx=np.where(rho==0)[0]
        if np.size(idx)!=0:
            vc[idx]=1
        return vc

    def _imageFunction(self,xx,yy,wl,t):
        return xx*0+1
    
###############################################################################
class oimUD(oimComponentFourier):
    name="Uniform Disk"
    shortname = "UD"
    def __init__(self,**kwargs):   
        super().__init__(**kwargs)         
        self.params["d"]=oimParam(**_standardParameters["d"])
        self._eval(**kwargs)

    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        xx=np.pi*self.params["d"](wl,t)*self.params["d"].unit.to(units.rad)*rho
        return np.nan_to_num(np.divide(2*j1(xx),xx),nan=1)

    
    def _imageFunction(self,xx,yy,wl,t):
        return ((xx**2+yy**2)<=(self.params["d"](wl,t)/2)**2).astype(float)
    
###############################################################################
class oimEllipse(oimUD):
    name="Uniform Ellipse"
    shortname = "eUD"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
###############################################################################   

class oimGauss(oimComponentFourier):
    name="Gaussian Disk"
    shortname = "GD"
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self.params["fwhm"]=oimParam(**_standardParameters["fwhm"])
         self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):
        return np.exp(-1*(np.pi*self.params["fwhm"](wl,t)*
               self.params["fwhm"].unit.to(units.rad)*rho)**2/(4*np.log(2)))
    
    def _imageFunction(self,xx,yy,wl,t):        
        r2=(xx**2+yy**2)
        return np.sqrt(4*np.log(2*self.params["fwhm"](wl,t))/np.pi)* \
                   np.exp(-4*np.log(2)*r2/self.params["fwhm"](wl,t)**2)

###############################################################################   
class oimEGauss(oimGauss):
    name="Gaussian Ellipse"
    shortname = "EG"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
###############################################################################

class oimIRing(oimComponentFourier):
    name="Infinitesimal Ring"
    shortname = "IR"
    def __init__(self,**kwargs):    
        super().__init__(**kwargs)   
        self.params["d"]=oimParam(**_standardParameters["d"])    
        self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):     
        xx=np.pi*self.params["d"](wl,t)*self.params["d"].unit.to(units.rad)*rho
        return j0(xx)
    
    def _imageFunction(self,xx,yy,wl,t,minPixSize=None):
        r2=(xx**2+yy**2)
        dx=np.max([np.abs(1.*(xx[0,0,0,1]-xx[0,0,0,0])),
                   np.abs(1.*(yy[0,0,1,0]-yy[0,0,0,0]))])
        return ((r2<=(self.params["d"](wl,t)/2+dx)**2) & 
                (r2>=(self.params["d"](wl,t)/2)**2)).astype(float)

###############################################################################   
class oimEIRing(oimIRing):
    name="Ellitical infinitesimal Ring"
    shortname = "EIR"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
         
###############################################################################

class oimRing(oimComponentFourier):
    name="Ring"
    shortname = "R"
    def __init__(self,**kwargs):        
        super().__init__(**kwargs)
        self.params["din"]=oimParam(**_standardParameters["din"]) 
        self.params["dout"]=oimParam(**_standardParameters["dout"])           
        self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):     
        xxin=np.pi*self.params["din"](wl,t)* \
                         self.params["din"].unit.to(units.rad)*rho
        xxout=np.pi*self.params["dout"](wl,t)* \
                         self.params["dout"].unit.to(units.rad)*rho
    
        fin=(self.params["din"](wl,t))**2
        fout=(self.params["dout"](wl,t))**2
           
        return np.nan_to_num(2*(j1(xxout)/xxout*fout-j1(xxin)/xxin*fin)/(fout-fin),nan=1)
          
    
    def _imageFunction(self,xx,yy,wl,t):

        r2=(xx**2+yy**2)
        return ((r2<=(self.params["dout"](wl,t)/2)**2) & 
                (r2>=(self.params["din"](wl,t)/2)**2)).astype(float)
    
###############################################################################

class oimRing2(oimComponentFourier):
    name="Ring2"
    shortname = "R2"
    def __init__(self,**kwargs):        
        super().__init__(**kwargs)
        self.params["d"]=oimParam(**_standardParameters["d"]) 
        self.params["w"]=oimParam(**_standardParameters["d"])           
        self.params["w"].name="w"
        self.params["w"].description="width of the ring"
        self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):
        
        d=self.params["d"](wl,t)* self.params["d"].unit.to(units.rad)
        w=self.params["w"](wl,t)* self.params["w"].unit.to(units.rad)
        
        
        xx=np.pi*(d+w/2)*rho
        dxx=np.pi*w*rho
    
        return j0(xx)*np.nan_to_num(np.divide(2*j1(dxx),dxx),nan=1)

    
    def _imageFunction(self,xx,yy,wl,t):

        r2=(xx**2+yy**2)
        return ((r2<=(self.params["d"](wl,t)/2+self.params["w"](wl,t)/2)**2 & 
                (r2>=(self.params["d"](wl,t)/2-self.params["w"](wl,t)/2)**2))).astype(float)
    
###############################################################################   
class oimERing(oimRing):
    name="Ellitical  Ring"
    shortname = "ER"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
         
###############################################################################   
class oimERing2(oimRing2):
    name="Ellitical  Ring2"
    shortname = "ER2"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)

                  
###############################################################################
          
class oimESKIRing(oimComponentFourier):
    name="Skewed Ellitical Infinitesimal Ring"
    shortname = "SKEIR"
    elliptic=True   
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["d"]=oimParam(**_standardParameters["d"])   
        self.params["skw"]=oimParam(**_standardParameters["skw"])    
        self.params["skwPa"]=oimParam(**_standardParameters["skwPa"])      
        self._eval(**kwargs)
    
    def _visFunction(self,xp,yp,rho,wl,t):

        xx=np.pi*self.params["d"](wl,t)*self.params["d"].unit.to(units.rad)*rho    
        
        phi=  (self.params["skwPa"](wl,t)-self.params["pa"](wl,t))* \
            self.params["skwPa"].unit.to(units.rad) +  np.arctan2(yp, xp);
       
        return np.nan_to_num(j0(xx)-I*np.sin(phi)*j1(xx)*self.params["skw"](wl,t),nan=1)
          
    def _imageFunction(self,xx,yy,wl,t):
        r2=(xx**2+yy**2)  
        # dr=np.sqrt(np.abs(np.roll(r2,(-1,-1),(0,1))-np.roll(r2,(1,1),(0,1))))
        phi=np.arctan2(yy,xx)  +  self.params["skwPa"](wl,t)* \
                                  self.params["skwPa"].unit.to(units.rad)
        dx=np.abs(1*(xx[0,0,0,1]-xx[0,0,0,0]
                    +xx[0,0,1,0]-xx[0,0,0,0])*self.params["elong"](wl,t))
        #dx=np.abs(1*(xx[0,1]-xx[0,0]+xx[1,0]-xx[0,0]))*3
        F=1+self.params["skw"](wl,t)*np.cos(phi)           
        return  ((r2<=(self.params["d"](wl,t)/2+dx/2)**2) & 
                 (r2>=(self.params["d"](wl,t)/2-dx/2)**2)).astype(float)*F
        
###############################################################################
          
class oimESKRing(oimComponentFourier):
    name="Skewed Ellitical Ring"
    shortname = "SKER"
    elliptic=True   
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["din"]=oimParam(**_standardParameters["din"]) 
        self.params["dout"]=oimParam(**_standardParameters["dout"])    
        self.params["skw"]=oimParam(**_standardParameters["skw"])    
        self.params["skwPa"]=oimParam(**_standardParameters["skwPa"])      
        self._eval(**kwargs)
    
    def _visFunction(self,xp,yp,rho,wl,t):
        
        
        xxin=np.pi*self.params["din"](wl,t)* \
                         self.params["din"].unit.to(units.rad)*rho
        xxout=np.pi*self.params["dout"](wl,t)* \
                         self.params["dout"].unit.to(units.rad)*rho
    
        xx=(xxin+xxout)/2
        xx2=(xxout-xxin)/2
        
        phi =  (self.params["skwPa"](wl,t)-self.params["pa"](wl,t))* \
            self.params["skwPa"].unit.to(units.rad) +  np.arctan2(yp, xp)
            

        res=(j0(xx)-1j*np.sin(phi)*j1(xx)*self.params["skw"](wl,t))  \
                *np.nan_to_num(np.divide(2*j1(xx2),xx2),nan=1)
     
        return res
        
        
          
    def _imageFunction(self,xx,yy,wl,t):
        r2=(xx**2+yy**2)  

        phi =  (self.params["skwPa"](wl,t)-self.params["pa"](wl,t))* \
                  self.params["skwPa"].unit.to(units.rad) + np.arctan2(yy,xx) 

        F=(1+self.params["skw"](wl,t)*np.sin(phi)  )/(1+self.params["skw"](wl,t))   
        
        return  ((r2<=(self.params["dout"](wl,t)/2)**2) & 
                 (r2>=(self.params["din"](wl,t)/2)**2)).astype(float)*F
    
###############################################################################    

class oimLorentz(oimComponentFourier):
    #From Lazareff 2017 A&A 599, 85
    #TODO : Small difference between images using direct formula or inverse of vis function
    name="Pseudo-Lorentzian"
    shortname = "LO"
    elliptic=True   
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["fwhm"]=oimParam(**_standardParameters["fwhm"])  
        self._eval(**kwargs)
    
    def _visFunction(self,xp,yp,rho,wl,t):
               
        xx=np.pi*self.params["fwhm"](wl,t)* self.params["fwhm"].unit.to(units.rad)*rho
        return np.exp(-2*np.pi*xx/3**0.5)
            
        
    def _imageFunction(self,xx,yy,wl,t):
        r2=(xx**2+yy**2)  
        a=np.pi*self.params["fwhm"](wl,t)* self.params["fwhm"].unit.to(units.mas)
        return a/(2*np.pi*3**0.5)*(a**2/3+r2)**(-1.5)
    
###############################################################################    
    
class oimELorentz(oimLorentz):
    name="Ellitical  Pseudo-Lorentzian"
    shortname = "ELO"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
    
###############################################################################    
   
class oimLinearLDD(oimComponentFourier):
    name="Linear Limb Darkened Disk "
    shortname = "LLDD" 
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["d"]=oimParam(**_standardParameters["d"])  
        self.params["a"]=oimParam(name="a",value=0,description="Linear LDD coeff",
                                  unit=units.one,mini=-1,maxi=1)
               
        self._eval(**kwargs)
    
    def _visFunction(self,xp,yp,rho,wl,t):
               
        xx=np.pi*self.params["d"](wl,t)* self.params["d"].unit.to(units.rad)*rho
        
        a=self.params["a"](wl,t)
  
        c1=2*np.divide(j1(xx),xx)
        c2=1.5*(np.pi*2)**0.5*np.divide(jv(1.5,xx),xx**1.5)    
        return  np.nan_to_num((1-a)*c1+a*c2,nan=1)
            
###############################################################################    

class oimQuadLDD(oimComponentFourier):
    #From Domiciano de Souza 2003 (phd thesis)
    name="Quadratic Limb Darkened Disk "
    shortname = "QLDD"  
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["d"]=oimParam(**_standardParameters["d"])  
        self.params["a1"]=oimParam(name="a1",value=0,description="1st QLDD coeff",
                                   unit=units.one,mini=-1,maxi=1)
        self.params["a2"]=oimParam(name="a2",value=0,description="2nd QLDD coeff",
                                   unit=units.one,mini=-1,maxi=1)
                
        self._eval(**kwargs)
    
    def _visFunction(self,xp,yp,rho,wl,t):
               
        xx=np.pi*self.params["d"](wl,t)* self.params["d"].unit.to(units.rad)*rho
        
        a1=self.params["a1"](wl,t)
        a2=self.params["a2"](wl,t)

        
        
        c1=np.divide(j1(xx),xx)
        c2=(np.pi/2)**0.5*np.divide(jv(1.5,xx),xx**1.5)
        c3=2*np.divide(jv(2.,xx),xx**2.)
        s=(6-2*a1-a2)/12
        return np.nan_to_num(((1-a1-a2)*c1+(a1+2*a2)*c2-a2*c3)/s,nan=1)
    
###############################################################################    

class oimConvolutor(oimComponentFourier):
    def __init__(self,component1, component2,**kwargs):        
        super().__init__(**kwargs)
        
        self.component1=component1
        self.component2=component2
        self.name="Convonlution Component"
        self.shortname = "Conv"
         
        for key in component2.params:
            self.params["c2_"+key]=component2.params[key]
            
        for key in component1.params:
               self.params["c1_"+key]=component1.params[key]
         
        self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):     
    
        return  self.component1.getComplexCoherentFlux(xp,yp,wl,t)* \
                self.component2._visFunction(xp,yp,rho,wl,t)
     
    def _imageFunction(self,xx,yy,wl,t):
        return self.component1._imageFunction(xx,yy,wl,t)* \
               self.component2._imageFunction(xx,yy,wl,t)



###############################################################################    
class oimModel(object):
    """
    The oimModel class hold a model made of one or more components (derived 
    from the oimComponent class), and allow to compute images (or image cubes 
    for wavelength or time dependent models) and complex coherent fluxes for a 
    vector of u,v,wl, and t coordinates.
    """
    def __init__(self,*components): 
        """
        Parameters
        ----------
        *components : list or serie of oimComponents
           The components of the model

        Returns
        -------
        None.
        """
        
    
        if len(components)==1 and type(components[0])==list:
            self.components=components[0]
        else:    
            self.components=components
        
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
        """
        Compute and return the complex coherent flux for an array of u,v 
        (and optionally wavelength and time ) coordinates.

        Parameters
        ----------
        u : list or numpy array
            spatial coordinate u (in cycles/rad) 
        v : list or numpy array
            spatial coordinate vu (in cycles/rad) .
        wl : list or numpy array, optional
            wavelength(s) in meter. The default is None.
        t :  list or numpy array, optional
            time in s (mjd). The default is None.

        Returns
        -------
        A numpy array of  the same size as u & v
            The complex coherent flux.

        """
        
        res=complex(0,0)
        for c in self.components:
            res+=c.getComplexCoherentFlux(ucoord,vcoord,wl,t)
   
        return res
        
    def getParameters(self,free=False): 
        """
        
        Get the Model paramters (or free parameters)
        
        Parameters
        ----------
        free : Bool, optional
            If True retrieve the free parameters of the models only. 
            The default is False.

        Returns
        -------
        params : Dict of oimParam
            a Dictionnary of the model parameters (or free parameters).

        """

        params={}
        for i,c in enumerate(self.components):
            for name,param in c.params.items():
                if not(param in params.values()):
                    if     isinstance(param,oimParamInterpWl) \
                        or isinstance(param,oimParamInterpTime) \
                        or isinstance(param,oimParamInterpolator):
                        for iparam,parami in enumerate(param.params):
                            if not(parami in params.values()):
                                if (parami.free==True or free==False):
                                    params["c{0}_{1}_{2}_interp{3}".format(
                                        i+1, c.shortname.replace(" ", "_"), 
                                        name, iparam+1)]=parami
                    elif isinstance(param,oimParamLinker):
                        pass
                    else:
                        if (param.free==True or free==False):
                            
                             params["c{0}_{1}_{2}".format(i+1, 
                                  c.shortname.replace(" ", "_"), name)]=param                           
        return params

    def getFreeParameters(self):  
        """
        Get the Model free paramters 

        Returns
        -------
        Dict of oimParam
            A Dictionnary of the model free parameters.
        """
        return self.getParameters(free=True)    



    def getImage(self,dim,pixSize,wl=None,t=None,toFits=False, 
                 fromFT=False,squeeze=True,normalize=False):
        """
        Compute and return an image or and image cube (if wavelength and time 
        are given). The returned image as the x,y dimension dim in pixel with
        an angular pixel size pixSize in rad. Image is returned as a numpy 
        array unless the keyword fits is set to True. In that case the image is
        returned as an astropy.io.fits hdu.

        Parameters
        ----------
        dim : integer
            image x & y dimension in pixels..
        pixSize : float
            pixel angular size.in mas
        wl : integer, list or numpy array, optional
            wavelength(s) in meter. The default is None.
        t :  integer, list or numpy array, optional
            time in s (mjd). The default is None.
        fits : bool, optional
            if True returns result as a fits hdu. The default is False.
        fromFT : bool, optional
            If True compute the image using FT formula when available
            The default is False.
        squeeze : bool, optional
            If False returns a (nt,nwl,dim,dim) array even if nt and/or nwl equal 1
            The default is True
            
        Returns
        -------
            numpy.ndarray or astropy.io.fits.hdu
             a numpy 2D array (or 3 or 4D array if wl, t or both are given) or an
             astropy.io.fits hdu.imagehdu if fits=True.
             The image of the component with given size in pixels and mas or rasd
        """
        
        #TODO : maybe we should change all None to zero as default values
        if wl is None:
            wl=0
        if t is None:
            t=0
         
        
        t=np.array(t).flatten()
        nt=t.size
        wl=np.array(wl).flatten()
        nwl=wl.size
        dims=(nt,nwl,dim,dim)
        
        if fromFT==True:

            v=np.linspace(-0.5,0.5,dim)
            vx,vy=np.meshgrid(v,v)
            
            vx_arr=np.tile(vx[None,None,:,:], (nt,nwl, 1, 1))
            vy_arr=np.tile(vy[None,None,:,:], (nt,nwl, 1, 1))
            wl_arr=np.tile(wl[None,:,None,None], (nt,1, dim, dim))
            t_arr=np.tile(t[:,None,None,None], (1,nwl, dim, dim))
            
            spfx_arr=(vx_arr/pixSize/mas2rad).flatten()   
            spfy_arr=(vy_arr/pixSize/mas2rad).flatten()   
            wl_arr=wl_arr.flatten()
            t_arr=t_arr.flatten()

            ft=self.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr,t_arr).reshape(dims)
            
            image=np.abs(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ft,axes=[-2,-1]),axes=[-2,-1]),axes=[-2,-1]))
           
        else:
            image=np.zeros(dims)
            for c in self.components:
                image+=c.getImage(dim,pixSize,wl,t)
                
        if normalize==True:
            for it in range(nt):
                for iwl in range(nwl):
                    image[it,iwl,:,:]/=np.max(image[it,iwl,:,:])
            
        #Always squeeze dim which are equal to one if exported to fits format
        if squeeze==True or toFits==True:
            image= np.squeeze(image)
            

        if toFits==True:
            
            hdu = fits.PrimaryHDU(image)
            hdu.header['CDELT1']=pixSize*mas2rad
            hdu.header['CDELT2']=pixSize*mas2rad
            hdu.header['CRVAL1']=0
            hdu.header['CRVAL2']=0
            hdu.header['CRPIX1']=dim/2
            hdu.header['CRPIX2']=dim/2
            hdu.header['CUNIT1']="rad"
            hdu.header['CUNIT2']="rad"
            hdu.header['CROTA1']=0
            hdu.header['CROTA2']=0 
            
            naxis=3
            if nwl!=1:
                dwl=(np.roll(wl,-1)-wl)[:-1]
            
                if np.all(np.abs(dwl-dwl[0])<1e-12):
                    dwl=dwl[0]
                    
                    hdu.header['CDELT{}'.format(naxis)]=dwl
                    hdu.header['CRPIX{}'.format(naxis)]=1
                    hdu.header['CRVAL{}'.format(naxis)]=wl[0]
                    hdu.header['CUNIT{}'.format(naxis)]="m" 
                    naxis+=1
                    
                else:
                    raise TypeError("Wavelength vector is not regular. Fit image" \
                                    " with irregular grid not yet implemented")
            
            if nt!=1:
                dt=(np.roll(t,-1)-t)[:-1]
            
                if np.all(np.abs(dt-dt[0])<1e-12):
                    dt=dt[0]
                    
                    hdu.header['CDELT{}'.format(naxis)]=dt
                    hdu.header['CRPIX{}'.format(naxis)]=1
                    hdu.header['CRVAL{}'.format(naxis)]=t[0]
                    hdu.header['CUNIT{}'.format(naxis)]="day"   
                    
                else:
                    raise TypeError("Time vector is not regular. Fit image" \
                                    " with irregular grid not yet implemented")  
            return hdu
        else:
            return image


    def saveImage(self,filename,dim,pixSize,wl=None,t=None,fromFT=False,normalize=False):
        im=self.getImage(dim,pixSize,wl=wl,t=t,toFits=True,
                         fromFT=fromFT,normalize=normalize)
        
        im.writeto(filename,overwrite=True)
        return im
    

    def showModel(self,dim,pixSize,wl=None,t=None, 
        fromFT=False,axe=None,normPow=0.5,figsize=(3.5,2.5),savefig=None,
        colorbar=True,legend=False,swapAxes=True,kwargs_legend={},
        normalize=False,**kwargs):
        """
        
        Show the mode Image or image-Cube

        Parameters
        ----------
        dim : integer
            image x & y dimension in pixels..
        pixSize : float
            pixel angular size.in mas
        wl : integer, list or numpy array, optional
            wavelength(s) in meter. The default is None.
        t :  integer, list or numpy array, optional
            time in s (mjd). The default is None.
        fits : bool, optional
            if True returns result as a fits hdu. The default is False.
        fromFT : bool, optional
            If True compute the image using FT formula when available
            The default is False.
        axe : matplotlib.axes.Axes, optional
            If provided the image will be shown in this axe. If not a new figure 
            will be created. The default is None.
        normPow : float, optional
            Exponent for the Image colorcale powerLaw normalisation.
            The default is 0.5.
        figsize : tuple of int, optional
            The Figure size in inches. The default is (8,6).
        savefig : str, optional
            Name of the files for saving the figure If None the figure is not saved.
            The default is None.
        colorbar : bool, optional
            Add a colobar to the Axe. The default is True.
        **kwargs : dict
            Arguments to be passed to the plt.imshow function

        Returns
        -------
        fig : matplotlib.figure.Figure
            The Figure created if needed
        axe : matplotlib.axes.Axes
            The Axes instances, created if needed.
        im  : the image(s) as a numpy array

        """

        im=self.getImage(dim,pixSize,wl,t,fromFT=fromFT,squeeze=False,normalize=normalize)
              
        t=np.array(t).flatten()      
        wl=np.array(wl).flatten()
              
        if swapAxes:
            t,wl=wl,t
            
        nt=t.size    
        nwl=wl.size
    
        if axe is None:
            fig,axe=plt.subplots(nwl,nt,figsize=(figsize[0]*nt,figsize[1]*nwl)
                ,sharex=True,sharey=True,subplot_kw=dict(projection='oimAxes'))    
        else:
            try:
                fig=axe.get_figure()
            except:
                fig=axe.flatten()[0].get_figure()
         

        axe=np.array(axe).flatten().reshape((nwl,nt))
        
        
        if not('norm' in kwargs):
            kwargs['norm']=colors.PowerNorm(gamma=normPow)
            
        
        for iwl,wli in enumerate(wl):
            for it,ti in enumerate(t):
                if swapAxes==False:
                    cb=axe[iwl,it].imshow(im[it,iwl,:,:],
                        extent=[-dim/2*pixSize,dim/2*pixSize,
                                -dim/2*pixSize,dim/2*pixSize],
                        origin='lower',**kwargs)
                else:
                    cb=axe[iwl,it].imshow(im[iwl,it,:,:],
                        extent=[-dim/2*pixSize,dim/2*pixSize,
                                -dim/2*pixSize,dim/2*pixSize],
                        origin='lower',**kwargs)
                
                axe[iwl,it].set_xlim(dim/2*pixSize,-dim/2*pixSize)
                
                
                if iwl==nwl-1:
                    axe[iwl,it].set_xlabel("$\\alpha$(mas)")
                if it==0:     
                    axe[iwl,it].set_ylabel("$\\delta$(mas)")
            
                if legend==True:
                    txt=""
                    if swapAxes==False:

                        if wl[0]!=None:
                            txt+="wl={:.4f}$\mu$m\n".format(wli*1e6)
                        if t[0]!=None:
                            txt+="Time={}".format(ti)
                        if not('color' in kwargs_legend):
                            kwargs_legend['color']="w"
                    else: 
                        if t[0]!=None:
                            txt+="wl={:.4f}$\mu$m\n".format(ti*1e6)
                        if wl[0]!=None:
                            txt+="Time={}".format(wli)
                        if not('color' in kwargs_legend):
                            kwargs_legend['color']="w"
                    axe[iwl,it].text(0,0.95*dim/2*pixSize,txt,
                            va='top',ha='center',**kwargs_legend)
                    
            
        if colorbar!=False:
            fig.colorbar(cb, ax=axe,label="Normalized Intensity")
            
            
        if savefig!=None:
            plt.savefig(savefig)
            
        return fig,axe,im

        







