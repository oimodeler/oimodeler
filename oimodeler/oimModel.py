# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""


class oimParam(object):
    
    def __init__(self):
        """
        

        Returns
        -------
        None.

        """
        pass
    
    

class oimComponent(object):
    """
    Parent abstract class for all OIModeler models components
    """    
    def __init__(self):
        """
        

        Returns
        -------
        None.

        """
        pass
    def getComplexComplexCoherentFlux(u,v,wl=None,t=None):
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
        return None;
    def getImage(dim,pixSize,wl=None,t=None,fits=False):
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
        return None;    

    
class oimModel(object):
    
    def __init__(self):
        pass
    
    