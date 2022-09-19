# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""

import numpy as np
from astropy import units as units
from scipy.special import j0,j1
from scipy import integrate
import numbers
from scipy import interpolate

np.seterr(invalid='ignore')
###############################################################################
"""
Useful definitions.
Maybeto be  put in a separate file later
"""
mas2rad=units.mas.to(units.rad)
I=complex(0,1)
###############################################################################

class oimParam(object):
    """
    Class of model parameters
    """
    def __init__(self,name=None,value=None,mini=-1*np.inf,maxi=np.inf,
                 description="",unit=1,free=True,error=0):
        """
        Create and initiliaze a new instance of the oimParam class
        Parameters
        ----------
        name : string, optional
            name of the Parameter. The default is None.
        value : float, optional
            value of the parameter. The default is None.
        mini : float, optional
            mininum value allowed for the parameter. The default is -1*np.inf.
        maxi : float, optional
            maximum value allowed for the parameter. The default is np.inf.
        description : string, optional
            A description of the parameter. The default is "".
        unit : 1 or astropy.unit, optional
            unit of the parameter. The default is 1.

        Returns
        -------
        None.

        """
        self.name=name
        self.value=value 
        self.error=error
        self.min=mini
        self.max=maxi
        self.free=free              
        self.description=description
        self.unit=unit

    
    def __call__(self,wl=None,t=None):
        """ The call function will be useful for wavelength or time dependent
        parameters. In a simple oimParam it only return the parameter value
        """
        return self.value
    
    
    def __str__(self):
        return "oimParam={}={}+-{} range=[{},{}] free={} ".format(self.name,
                self.value,self.error,self.min,self.max,self.free)

    def __repr__(self):
        return "oimParam at {} : {}={}+-{} range=[{},{}] free={} ".format(hex(id(self)),self.name,
                self.value,self.error,self.min,self.max,self.free)


###############################################################################

class oimInterpWl(object):
    """
    Structure for creating oimParamInterpWl directly in oimParam defintion
    """
    def __init__(self,wl=[],value=None):
        self.wl=wl
        self.value=value

###############################################################################
        
class oimParamInterpWl(oimParam):
    def __init__(self,param,interpWl):
        
        self.name=param.name
        self.description=param.description
        self.unit=param.unit
        
        value= interpWl.value
        wl=interpWl.wl
        nwl=len(wl)
      
        if value==None:
            value=[self.value]*nwl
        elif isinstance(value,numbers.Number):
            value=[value]*nwl
        else:
            if len(value)!=nwl:
                raise TypeError("wl and val should have the same length :"  
                            "len(x)={}, len(y)={}".format(len(wl), len(value)))
                
        self.params=[]
        self._wl=wl
        self._nwl=nwl
        
        for i in range(nwl):
           
            pi=oimParam(name=param.name,value=value[i],mini=param.min,
                        maxi=param.max,description=param.description,
                        unit=param.unit,free=param.free,error=param.error)
            self.params.append(pi)
            


                    

    def __call__(self,wl=None,t=None):
        values=np.array([pi.value for pi in self.params])
        return np.interp(wl,self.wl,values,left=values[0], right=values[-1])
  
             
        
    @property
    def wl(self):
        return self._wl
    
    @wl.setter
    def wl(self,_wl):
        nwl=len(_wl)
        if nwl== self._nwl:
            self._wl=np.array(_wl)
        else:
             raise TypeError("Can't modify number of key wls in oimParamInterpWl "
                             "after creation. Consider creating a new parameter")    

###############################################################################
class oimParamLinker(object):
    def __init__(self,param,operator="add",fact=0):
        self.param=param    
        self.fact=fact        
        
        self.op=None
        self._setOperator(operator)
       
          
    @property
    def unit(self):
        return self.param.unit
    
       

    def _setOperator(self,operator):
         if operator=="add":
             self.op=self._add
         elif operator=="mult":
             self.op=self._mult            
             
    def _add(self,val):
        return val+ self.fact
    
    def _mult(self,val):
        return val+ self.fact
    
        
    def __call__(self,wl=None,t=None):  
        return self.op(self.param.__call__(wl,t))
        
        

###############################################################################


    
#Here is a list of standard parameters to be used when defining new components
_standardParameters={
    "x":{"name":"x","value":0,"description":"x position","unit":units.mas,"free":False},
    "y":{"name":"x","value":0,"description":"y position","unit":units.mas,"free":False},
    "f":{"name":"f","value":1,"description":"flux","unit":1},
    "fwhm":{"name":"fwhm","value":0,"description":"FWHM","unit":units.mas},
    "d":{"name":"d","value":0,"description":"Diameter","unit":units.mas},
    "din":{"name":"din","value":0,"description":"Inner Diameter","unit":units.mas},
    "dout":{"name":"dout","value":0,"description":"Outer Diameter","unit":units.mas},    
    "elong":{"name":"elong","value":1,"description":"Elongation Ratio","unit":1},
    "pa":{"name":"pa","value":0,"description":"Major-axis Position angle","unit":units.deg},
    "skw":{"name":"skw","value":0,"description":"Skewedness","unit":"1"},
    "skwPa":{"name":"skwPa","value":0,"description":"Skewedness Position angle","unit":units.deg},
    "pixSize":{"name":"pixSize","value":0.1,"description":"Pixel Size","unit":units.mas}
    }
    

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
                if isinstance(value,numbers.Number):
                    self.params[key].value=value
                elif isinstance(value,oimInterpWl):
                    if not(isinstance(self.params[key],oimParamInterpWl)):
                        self.params[key]=oimParamInterpWl(self.params[key],value)
                    
                    
                            
    def getComplexComplexCoherentFlux(self,u,v,wl=None,t=None):
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
    

    
    def __str__(self):
        txt=self.name
        for name,param in self.params.items():  
            txt+=" {0}={1:.2f}".format(param.name,param.value)
        return txt
    

###############################################################################
 
class oimComponentImage(oimComponent):    
    """
    Base class for components define in 2D : x,y (regular grid) in the image plan.
    """
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)

        self.dim=0
        self.pixSize=0
        
        self.params["pixSize"]=oimParam(**_standardParameters["pixSize"])  
        self.params["dim"]=oimParam(name="dim",value=256,
                         description="Image dimension",unit=1,free=False)
        
        self._eval(**kwargs)   
        self._compute2DCoordGrid(force=True)
        
      
    def _compute2DCoordGrid(self,force=False):
        dim=self.params["dim"].value
        pixSize=self.params["pixSize"].value 
        
        #recompute grid only if pixSize or dim have changed
        
        if (dim!=self.dim or pixSize!=self.pixSize ):
           self.xx=np.tile((np.arange(dim)-dim/2)*pixSize,(dim,1))
           self.yy=np.transpose(self.xx)
           self.dim=dim
           self.pixSize=pixSize        

    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
        
        self._compute2DCoordGrid()        
        image= self._imageFunction(self.xx,self.yy,wl,t)
        s=np.shape(image) 
        
        fft2D=np.fft.fftshift(np.fft.fft2(np.fft.fftshift(image),
                                s=[s[0],s[1]]))

        pixSize=self.params["pixSize"]()*self.params["pixSize"].unit.to(units.rad)
        
        freqVectX=np.fft.fftshift(np.fft.fftfreq(s[0],pixSize))
        freqVectY=np.fft.fftshift(np.fft.fftfreq(s[1],pixSize))
        
        grid=(freqVectX,freqVectY)
        coord=np.transpose([ucoord,vcoord])
               
        real=interpolate.interpn(grid,np.real(fft2D),coord,method='linear',bounds_error=False,fill_value=None)
        imag=interpolate.interpn(grid,np.imag(fft2D),coord,method='linear',bounds_error=False,fill_value=None)
        
        
        vc=real+imag*np.complex(0,1)
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
          
    def _imageFunction(self,xx,yy,wl,t):
        dim=self.params["dim"].value
        return np.zeros([dim,dim])


    def _directTranslate(self,x,y,wl,t):
         x=x-self.params["x"](wl,t)
         y=y-self.params["y"](wl,t)
         return x,y

    def getImage(self,dim,pixSize,wl=None,t=None):  
           x=np.tile((np.arange(dim)-dim/2)*pixSize,(dim,1))
           y=np.transpose(x)
           x,y=self._directTranslate(x,y,wl,t)
           
           image = self._imageFunction(x,y,wl,t)
           tot=np.sum(image)
           if tot!=0:  
               image = image / np.sum(image,axis=(0,1))*self.params["f"](wl,t)    

           return image
                                 
    def getInternalImage(self,wl=None,t=None):  
           self._compute2DCoordGrid()        
           image= self._imageFunction(self.xx,self.yy,wl,t)
           tot=np.sum(image)
           if tot!=0:  
               image = image / np.sum(image,axis=(0,1))*self.params["f"](wl,t)    

           return image
                                                
###############################################################################
    

class oimComponentChromaticCube(oimComponent):    
    """
    Base class for components define in 3D : x,y (regular) and wl (non-regular) 
    in Image plan.
    """
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.params["pixSize"]=oimParam(**_standardParameters["pixSize"])
        self.image=None
        self.dim=0
        self.wl_vector=None
        self._eval(**kwargs)      
      
    def setImage(self,image,wl):
        self.image=image
        s=np.shape(image)
        
        self.dim=s[-1]
        self.wl_vector=wl
        self.nwl=s[0]
      
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
              
        fft2D=np.fft.fftshift(np.fft.fft2(np.fft.fftshift(self.image,axes=[2,1]),
                                s=[self.dim,self.dim]),axes=(2,1))

        pixSize=self.params["pixSize"]()*self.params["pixSize"].unit.to(units.rad)
        freqVect=np.fft.fftshift(np.fft.fftfreq(self.dim,pixSize))
        
        grid=(self.wl_vector,freqVect,freqVect)
        coord=np.transpose([wl,ucoord,vcoord])
               
        real=interpolate.interpn(grid,np.real(fft2D),coord,method='linear',bounds_error=False,fill_value=None)
        imag=interpolate.interpn(grid,np.imag(fft2D),coord,method='linear',bounds_error=False,fill_value=None)
        
        
        vc=real+imag*np.complex(0,1)
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    
        

###############################################################################
    

class oimComponentRadialProfile(oimComponent):    
    """
    Base class for components define by a non-regular radial profile in Image plan.
    """
    def __init__(self,**kwargs): 
        super().__init__(**kwargs)
        self.hankel=self._trapeze     
        self._eval(**kwargs)      
      

    @staticmethod
    def _trapeze(Ir,r,sfreq):
        res = np.zeros_like(sfreq)    
        for i, sf in enumerate(sfreq):       
            res[i] = integrate.trapezoid(2.*np.pi*r*Ir*j0(2.*np.pi*r*sf), r)     
            
        flux = integrate.trapezoid(2.*np.pi*r*Ir, r)
        res = res/flux
        return res
    
    
    def getComplexCoherentFlux(self,ucoord,vcoord,wl=None,t=None):
             
        spf=np.sqrt(ucoord**2+vcoord**2)
        radius,intensity=self._radialProfileFunction()
        vc=self.hankel(intensity, radius,spf)
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    
    def _radialProfileFunction(self):
        return 0,0



###############################################################################


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
            pa_rad=(self.params["pa"](wl,t)+90)* \
                        self.params["pa"].unit.to(units.rad)      
            co=np.cos(pa_rad)
            si=np.sin(pa_rad)
            fxp=ucoord*co-vcoord*si
            fyp=ucoord*si+vcoord*co
            rho=np.sqrt(fxp**2+fyp**2/self.params["elong"](wl,t)**2) 
        else:
            fxp=ucoord
            fyp=vcoord
            rho=np.sqrt(fxp**2+fyp**2)            
               
        vc=self._visFunction(fxp,fyp,rho,wl,t)
               
        return vc*self._ftTranslateFactor(ucoord,vcoord,wl,t)* \
                                                     self.params["f"](wl,t)
    

          
    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        return ucoord*0
    
    def _directTranslate(self,x,y,wl,t):
        x=x-self.params["x"](wl,t)
        y=y-self.params["y"](wl,t)
        return x,y
    
    def getImage(self,dim,pixSize,wl=None,t=None):  
        x=np.tile((np.arange(dim)-dim/2)*pixSize,(dim,1))
        y=np.transpose(x)
        x,y=self._directTranslate(x,y,wl,t)
        
        if self.elliptic:
            pa_rad=(self.params["pa"](wl)+90)* \
                               self.params["pa"].unit.to(units.rad)
            xp=x*np.cos(pa_rad)-y*np.sin(pa_rad)
            yp=x*np.sin(pa_rad)+y*np.cos(pa_rad)
            x=xp
            y=yp*self.params["elong"].value
        image = self._imageFunction(x,y,wl,t)
        tot=np.sum(image)
        if tot!=0:  
            image = image / np.sum(image,axis=(0,1))*self.params["f"](wl,t)    

        return image
    
    def _imageFunction(self,xx,yy,wl=None):
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
         super().__init__(self,**kwargs)
         self._eval(**kwargs)

    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        vc=rho*0
        idx=np.where(rho==0)[0]
        if np.size(idx)!=0:
            vc[idx]=1
        return vc

    def _imageFunction(self,xx,yy,wl):
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
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self.name="Infinitesimal Ring"
         self.shortname = "IR"
         self.params["d"]=oimParam(**_standardParameters["d"])    
         self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):     
        xx=np.pi*self.params["d"](wl,t)*self.params["d"].unit.to(units.rad)*rho
        return j0(xx)
    
    def _imageFunction(self,xx,yy,wl,t):
        r2=(xx**2+yy**2)   
        dx=np.max([np.abs(1.*(xx[0,1]-xx[0,0])),np.abs(1.*(yy[1,0]-yy[0,0]))])
        return ((r2<=(self.params["d"](wl,t)/2+dx)**2) & 
                (r2>=(self.params["d"](wl,t)/2)**2)).astype(float)

###############################################################################   
class oimEIRing(oimIRing):
    name="Ellitical infinitesimal ring"
    shortname = "EIR"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
         
###############################################################################

class oimRing(oimComponentFourier):
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self.name="Ring"
         self.shortname = "R"
         self.params["din"]=oimParam(**_standardParameters["din"]) 
         self.params["dout"]=oimParam(**_standardParameters["din"])           
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
        #TODO Fix the elong pb!
        if self.elliptic==True:
            r2=(xx**2+yy**2*self.params["elong"](wl,t)**2) 
        else:
            r2=(xx**2+yy**2)
        return ((r2<=(self.params["dout"](wl,t)/2)**2) & 
                (r2>=(self.params["din"](wl,t)/2)**2)).astype(float)

###############################################################################   
class oimERing(oimRing):
    name="Ellitical  ring"
    shortname = "ER"
    elliptic=True
    def __init__(self,**kwargs):        
         super().__init__(**kwargs)
         self._eval(**kwargs)
                  
###############################################################################
          
class ESKRing(oimComponentFourier):
    name="Skewed Ellitical ring"
    shortname = "SKER"
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
       
        return j0(xx)-I*np.sin(phi)*j1(xx)*self.params["skw"](wl,t)
          
    def _imageFunction(self,xx,yy,wl,t):
        r2=(xx**2+yy**2)  
        # dr=np.sqrt(np.abs(np.roll(r2,(-1,-1),(0,1))-np.roll(r2,(1,1),(0,1))))
        phi=np.arctan2(yy,xx)  +  self.params["skwPa"](wl,t)* \
                                  self.params["skwPa"].unit.to(units.rad)
        dx=np.abs(1*(xx[0,1]-xx[0,0]+xx[1,0]-xx[0,0])*self.params["elong"](wl,t))
        #dx=np.abs(1*(xx[0,1]-xx[0,0]+xx[1,0]-xx[0,0]))*3
        F=1+self.params["skw"](wl,t)*np.cos(phi)           
        return  ((r2<=(self.params["d"](wl,t)/2+dx/2)**2) & 
                 (r2>=(self.params["d"](wl,t)/2-dx/2)**2)).astype(float)*F
        

###############################################################################    
class oimModel(object):
    """
    The oimModel class hold a model made of one or more components (derived 
    from the oimComponent class), and allow to compute images (or image cubes 
    for wavelength or time dependent models) and complex coherent fluxes for a 
    vector of u,v,wl, and t coordinates.
    """
    def __init__(self,components=[]): 
        """
        Create and initialize a new instance of oimModel

        Returns
        -------
        None.

        """
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
        
        return None;
    def getImage(self,dim,pixSize,wl=None,t=None,fits=False, fromFT=False):
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
        
        if fromFT:
            #TODO multiple wavelengths and time
            if wl==None:
                wl=np.array([1e-6])
            wl=np.array(wl)
            try:
                nwl=len(wl)    
            except:
                wl=np.array([wl])
                nwl=1
            
            Bmax=wl/pixSize/mas2rad #meter
           
            # TODO remplace loop by better stuff
            ims=np.ndarray([nwl,dim,dim])
            for iwl,wli in enumerate(wl):
                B=np.linspace(-Bmax[iwl]/2,Bmax[iwl]/2,dim)
                spfy=np.outer(B,B*0+1)/wli
                spfx=np.transpose(spfy)
                spfy=spfy.flatten()
                spfx=spfx.flatten()        
                ft=self.getComplexCoherentFlux(spfx,spfy,np.ones(dim*dim)*wli).reshape(dim,dim)
                
                ims[iwl,:,:]=np.abs(np.fft.fftshift(np.fft.ifft2(ft))) 
                
                if nwl==1:
                    ims=ims[0,:,:]
            return ims
        
        else:
            #TODO add time (which can lead to 4D hypercube images!!)
            #Used if wl is a scalar
            wl=np.array(wl)
            nwl=np.size(wl)

            if np.shape(wl)==():
                image=np.zeros([dim,dim])
                for c in self.components:
                    image+=c.getImage(dim,pixSize,wl)
            else:
                #TODO : this is very slow!!!
                image=np.zeros([nwl,dim,dim])
                for iwl,wli in enumerate(wl):
                    for c in self.components:
                        image[iwl,:,:]+=c.getImage(dim,pixSize,wli)

            return image;    

    def getParameters(self,free=False):       
        params={}
        for i,c in enumerate(self.components):
            for name,param in c.params.items():
                if not(param in params.values()):
                    if isinstance(param,oimParamInterpWl):
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
        return self.getParameters(free=True)    











