# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 15:45:17 2019

@author: ame
"""
from datetime import datetime

import numpy as np


def fastRotator  (dim0,size,incl,rot,Tpole,lam,beta=0.25):

    h=6.63e-34
    c=3e8
    kb=1.38e-23

    a=2./3*(rot)**0.4+1e-9
    K=np.sin(1./3.)*np.pi


    K1=h*c/kb
    nlam=np.size(lam)
    incl=np.deg2rad(incl)

    x0 = np.linspace(-size,size,num=dim0)
    idx=np.where(np.abs(x0)<=1.5)
    x=np.take(x0,idx)
    dim=np.size(x)
    unit=np.ones(dim)
    x = np.outer(x,unit)
    x = np.einsum('ij, k->ijk', x, unit)

    y=np.swapaxes(x,0,1)
    z=np.swapaxes(x,0,2)


    yp=y*np.cos(incl)+z*np.sin(incl)
    zp=y*np.sin(incl)-z*np.cos(incl)


    r=np.sqrt(x**2+yp**2+zp**2)


    theta=np.arccos(zp/r)

    x0=(1.5*a)**1.5*np.sin(1e-99)
    r0=a*np.sin(1/3.)*np.arcsin(x0)/(1.0/3.*x0)

    x2=(1.5*a)**1.5*np.sin(theta)
    rin=a*np.sin(1/3.)*np.arcsin(x2)/(1.0/3.*x2)


    rhoin=rin*np.sin(theta)/a/K

    dr=(rin/r0-r)>=0
    #dr=(rin/(r0*1.5)*1.5-r)>=0


    Teff=Tpole*(np.abs(1-rhoin*a)**beta)

    #TODO : implement a correct limb-darkening law
    #limb=np.abs((np.cos(np.arctan2(np.sqrt(x**2+y**2),-np.abs(z)))))*0+1


    if nlam==1:
        flx=1./(np.exp(K1/(lam*Teff))-1)

        im=np.zeros([dim,dim])

        for iz in range(dim):
            im=im*(im!=0)+(im==0)*dr[:,:,iz]*flx[:,:,iz]#*limb[:,:,iz]

        im=np.rot90(im)

        tot=np.sum(im)
        im=im/tot
        im0=np.zeros([dim0,dim0])

        im0[dim0//2-dim//2:dim0//2+dim//2,dim0//2-dim//2:dim0//2+dim//2]=im

        return im0

    else:
        unit=np.zeros(nlam)+1
        dr=np.einsum('ijk, l->ijkl', dr, unit)
        flx = 1./(np.exp(K1/np.einsum('ijk, l->ijkl',Teff,lam))-1)

        im=np.zeros([dim,dim,nlam])


        for iz in range(dim):
            im=im*(im!=0)+dr[:,:,iz,:]*flx[:,:,iz,:]*(im==0)


        im=np.rot90(im)

        tot=np.sum(im,axis=(0,1))


        for ilam in range(nlam):
            im[:,:,ilam]=im[:,:,ilam]/tot[ilam]

        im0=np.zeros([dim0,dim0,nlam])
        im0[dim0//2-dim//2:dim0//2+dim//2,dim0//2-dim//2:dim0//2+dim//2,:]=im

        return im0










