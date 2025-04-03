# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import numpy as np
from astropy import units as units

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam


def fastRotator(dim0, size, incl, rot, Tpole, lam, beta=0.25,a1=0,a2=0,ldd=None):
    
    """
    Equations are taken from Domiciano+ 2018
    https://www.aanda.org/articles/aa/pdf/2018/11/aa33450-18.pdf
    """

    nlam = np.size(lam)
    incl = np.deg2rad(incl)

    x0 = np.linspace(-size, size, num=dim0)
    idx = np.where(np.abs(x0) <= 1.5)
    x = np.take(x0, idx)
    dim = np.size(x)
    unit = np.ones(dim)
    x = np.outer(x, unit)
    x = np.einsum('ij, k->ijk', x, unit)

    y = np.swapaxes(x, 0, 1)
    z = np.swapaxes(x, 0, 2)

    yp = y*np.cos(incl)+z*np.sin(incl)
    zp = y*np.sin(incl)-z*np.cos(incl)

    r = np.sqrt(x**2+yp**2+zp**2)

    theta = np.arccos(zp/r)

    eps = rot**2/3 # flatening parameter
    ome = 1.5*(1-eps)*np.sqrt(3*eps) # angular rate

    Rtheta = (1-eps)*np.sin(1/3*np.arcsin(ome*np.sin(theta)))/(1/3*ome*np.sin(theta))
    Rtheta = Rtheta/Rtheta.min()
    Req =  Rtheta.max()
    
    dr = (Rtheta-r) >= 0
    
    
    Fc = rot**2*1/Req**2*Rtheta*np.sin(theta)
    
    G = 3/(2*Rtheta**2)
    
    G_z = G*np.cos(theta)
    G_rho = G*np.sin(theta)
    
    
    geff=np.sqrt(G_z**2+(Fc-G_rho)**2)
    Teff = Tpole*geff**beta

    
    mu=np.rot90(np.sum(dr,axis=2))
    mu=mu/mu.max()
    
    
    h = 6.63e-34
    c = 3e8
    kb = 1.38e-23
    K1 = h*c/kb

    if ldd == "linear":
        ldd_im = (1-a1*(1-mu))
    elif ldd == "quadratic":
        ldd_im = (1-a1*(1-mu)-a2*(1-mu)**2)
    else:
        ldd_im = 1
    
    if nlam == 1:
        flx = 1./(np.exp(K1/(lam*Teff))-1)*2*h*c**2/lam**5

        im = np.zeros([dim, dim])

        for iz in range(dim):
            im = im*(im != 0)+(im == 0) * \
                dr[:, :, iz]*flx[:, :, iz] 

        im = np.rot90(im)
        im = im * ldd_im


        tot = np.sum(im)
        im = im/tot
        im0 = np.zeros([dim0, dim0,1])

        im0[dim0//2-dim//2:dim0//2+dim//2, dim0//2-dim//2:dim0//2+dim//2,0] = im

        
        return im0

    else:
        unit = np.zeros(nlam)+1
        dr = np.einsum('ijk, l->ijkl', dr, unit)
        Teff2=Teff[:,:,:,np.newaxis]
        lam2 = lam [np.newaxis,np.newaxis,np.newaxis,:]
        flx = 1./(np.exp(K1/(lam2*Teff2))-1)*2*h*c**2/lam2**5
        #flx = 1./(np.exp(K1/np.einsum('ijk, l->ijkl', Teff, lam))-1)

        im = np.zeros([dim, dim, nlam])

        for iz in range(dim):
            im = im*(im != 0)+dr[:, :, iz, :]*flx[:, :, iz, :]*(im == 0)

        im = np.rot90(im)
        if ldd :
            im = im * ldd_im[:,:,np.newaxis]
        
        tot = np.sum(im, axis=(0, 1))

        for ilam in range(nlam):
            im[:, :, ilam] = im[:, :, ilam]/tot[ilam]

        im0 = np.zeros([dim0, dim0, nlam])
        im0[dim0//2-dim//2:dim0//2+dim//2, dim0//2-dim//2:dim0//2+dim//2, :] = im
        return im0



###############################################################################

class oimFastRotator(oimComponentImage):
    name = "Fast Rotator"
    shortname = "FRot"

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)

        # Component parameters. Note that as it inherits from the oimComponentImage class it already has
        # x,y,f and dim as parameters
        self.params["incl"] = oimParam(name="incl", value=0, description="Inclination angle", unit=units.deg)
        self.params["rot"] = oimParam(name="rot", value=0, description="Rotation Rate", unit=units.one)
        self.params["Tpole"] = oimParam(name="Tpole", value=20000, description="Polar Temperature", unit=units.K)
        self.params["dpole"] = oimParam(name="dplot", value=1, description="Polar diameter", unit=units.mas)
        self.params["beta"] = oimParam(name="beta", value=0.25, description="Gravity Darkening Exponent", unit=units.one)

        # constant value <=> static model
        self._t = np.array([0])

        # The component is chromatic. Here we set a fixed array of reference wavelengths. This can be
        # modified later as, in our case the model is recomputed at each call to the fastRotator function
        self._wl = np.linspace(0.5e-6, 15e-6, num=10)

        # Finally evalutating paramters as for all other components
        self._eval(**kwargs)

    def _internalImage(self):
        dim = self.params["dim"].value
        incl = self.params["incl"].value
        rot = self.params["rot"].value
        Tpole = self.params["Tpole"].value
        dpole = self.params["dpole"].value
        beta = self.params["beta"].value

        im = fastRotator(dim, 1.5, incl, rot, Tpole, self._wl, beta=beta)

        # make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im = np.tile(np.moveaxis(im, -1, 0)[None, :, :, :], (1, 1, 1, 1))

        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize = 1.5*dpole/dim*units.mas.to(units.rad)

        return im

class oimFastRotatorLLDD(oimComponentImage):
    name = "Fast Rotator"
    shortname = "FRot"

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)

        # Component parameters. Note that as it inherits from the oimComponentImage class it already has
        # x,y,f and dim as parameters
        self.params["incl"] = oimParam(name="incl", value=0, description="Inclination angle", unit=units.deg)
        self.params["rot"] = oimParam(name="rot", value=0, description="Rotation Rate", unit=units.one)
        self.params["Tpole"] = oimParam(name="Tpole", value=20000, description="Polar Temperature", unit=units.K)
        self.params["dpole"] = oimParam(name="dplot", value=1, description="Polar diameter", unit=units.mas)
        self.params["beta"] = oimParam(name="beta", value=0.25, description="Gravity Darkening Exponent", unit=units.one)
        self.params["a"] = oimParam(name="a", value=0, description="Linear LDD coeff",unit=units.one, mini=-1, maxi=1)

        # constant value <=> static model
        self._t = np.array([0])

        # The component is chromatic. Here we set a fixed array of reference wavelengths. This can be
        # modified later as, in our case the model is recomputed at each call to the fastRotator function
        self._wl = np.linspace(0.5e-6, 15e-6, num=10)

        # Finally evalutating paramters as for all other components
        self._eval(**kwargs)

    def _internalImage(self):
        dim = self.params["dim"].value
        incl = self.params["incl"].value
        rot = self.params["rot"].value
        Tpole = self.params["Tpole"].value
        dpole = self.params["dpole"].value
        beta = self.params["beta"].value
        a = self.params["a"].value
        
        im = fastRotator(dim, 1.5, incl, rot, Tpole, self._wl, beta=beta,ldd="linear",a1=a)

        # make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im = np.tile(np.moveaxis(im, -1, 0)[None, :, :, :], (1, 1, 1, 1))

        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize = 1.5*dpole/dim*units.mas.to(units.rad)

        return im
    
    
class oimFastRotatorQuadLDD(oimComponentImage):
    name = "Fast Rotator"
    shortname = "FRot"

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)

        # Component parameters. Note that as it inherits from the oimComponentImage class it already has
        # x,y,f and dim as parameters
        self.params["incl"] = oimParam(name="incl", value=0, description="Inclination angle", unit=units.deg)
        self.params["rot"] = oimParam(name="rot", value=0, description="Rotation Rate", unit=units.one)
        self.params["Tpole"] = oimParam(name="Tpole", value=20000, description="Polar Temperature", unit=units.K)
        self.params["dpole"] = oimParam(name="dplot", value=1, description="Polar diameter", unit=units.mas)
        self.params["beta"] = oimParam(name="beta", value=0.25, description="Gravity Darkening Exponent", unit=units.one)
        self.params["a1"] = oimParam(name="a1", value=0, description="1st SLDD coeff",unit=units.one, mini=-1, maxi=1)
        self.params["a2"] = oimParam(name="a2", value=0, description="2nd SLDD coeff",unit=units.one, mini=-1, maxi=1)
        # constant value <=> static model
        self._t = np.array([0])

        # The component is chromatic. Here we set a fixed array of reference wavelengths. This can be
        # modified later as, in our case the model is recomputed at each call to the fastRotator function
        self._wl = np.linspace(0.5e-6, 15e-6, num=10)

        # Finally evalutating paramters as for all other components
        self._eval(**kwargs)

    def _internalImage(self):
        dim = self.params["dim"].value
        incl = self.params["incl"].value
        rot = self.params["rot"].value
        Tpole = self.params["Tpole"].value
        dpole = self.params["dpole"].value
        beta = self.params["beta"].value
        a1 = self.params["a1"].value
        a2 = self.params["a2"].value
        
        im = fastRotator(dim, 1.5, incl, rot, Tpole, self._wl, beta=beta,ldd="linear",a1=a1,a2=a2)

        # make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im = np.tile(np.moveaxis(im, -1, 0)[None, :, :, :], (1, 1, 1, 1))

        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize = 1.5*dpole/dim*units.mas.to(units.rad)

        return im