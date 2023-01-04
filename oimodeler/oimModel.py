# -*- coding: utf-8 -*-
"""
creation of models

"""

import numpy as np
from astropy.io import fits
from astropy import units as units
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from oimodeler import oimParamLinker,oimParamInterpolator

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
                    if  isinstance(param,oimParamInterpolator):
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
            
            spfx_arr=(vx_arr/pixSize/units.mas.to(units.rad)).flatten()   
            spfy_arr=(vy_arr/pixSize/units.mas.to(units.rad)).flatten()   
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
            hdu.header['CDELT1']=pixSize*units.mas.to(units.rad)
            hdu.header['CDELT2']=pixSize*units.mas.to(units.rad)
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

        







