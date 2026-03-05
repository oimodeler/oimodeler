# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:29:45 2016

@author: domicia

DISCO: DIsc and Star COntinuum 

Model of star+gas disc (continuuum IR flux)
Based on:
Vieira et al. 2015, MNRAS, 454, 2107
Carciofi & Bjorkman 2006, ApJ, 639, 1081

NOTE: Computing time increases as N**2 and as Nwlen/2
"""

"""
Circumstellar disc:

Rd = disc outer radius (Rsun)
Bounds: Rd > Rstar

Td0 = temperature at the basis of the disc (K)
Bounds: Td0<Tstar
Suggested value for Be stars: Td0=0.6*Tstar

rho0 = density at disc basis (kg/m3)
Bounds: rho0 > 0

powrho = disc density power
Value for a steady state isothermal outflow: powrho=-3.5

powTd = disc-temperature power (T)
Bounds: powTd <= 0
Value for a reprocessing disk (optically thick, geometrically thin): powTd=-0.75 
Suggested value for inner regions of Be discs: powTd=-0.75

powHd = disc-flaring power
Value for a reprocessing disk (optically thick, geometrically thin): powHd=9/8 
Value for a thin isothermal disc: powHd=1.5
Bounds: powHd>0

ionfrac = ionization fraction
Value for a fully ionized disc : ionfrac = 1

Central star:

* Spherical star and blackbody model

Rstar = stellar radius (Rsun)
Bounds: Rstar > 0

Tstar = effective temperature of the central star (K)
Bounds: Tstar > 0

Mstar = stellar mass (Msun)
Bounds: Mstar > 0


Geometrical parameters:
i = inclination of star-disc (deg)

#Only necessary outside disco (to create a fits file for example)
d = distance (pc)
or
angsize = total angular size of the image (mas)

Spectral range
Wavelength values (different units accepted)
Bounds: Wavelength > 0
"""


"""
TODO:

*) Check all equations, write them down in a document

*) Check/resolve the question of mean molecular weight and ionization fraction

*) Criar opçoes para diferentes tipos de estrelas: uniform, limb darkenend,
elliptical, Roche, Kurucz

*) Introduce the free-free and bound-free opacities interpolated from 
my routines

*) Include option to force intensity to go to zero at the image borders

*) Output images as fits files (JMMC compatible, for uv computations)

*) Compare final version with Rodrigo's papers and FRACS II

*) Introduce Thomson (electronic) scattering ? See Lamers&Waters84 and 
email from Rodrigo (2017/26/01)

*) Complete list of input parameters and examples in the beginning of this routine
    
*) Add dust: (a) BB as in Beckwith 1999; (b) opacity form Mie (see email Gilles 2017-09-22)


"""

import numpy as np
#from astropy.analytic_functions import blackbody_lambda
#from astropy.constants import R_sun #http://docs.astropy.org/en/v0.2.1/constants/index.html
#from astropy.constants import M_sun

#CREATE DICTIONARIES WITH THE PARAMETERS: PHYSICAL, STAR, DISC

#Physical constants
cLight=299792458. #m/s (speed of light in vacuum); NIST source
mH=1.660539040E-27  #kg (atomic mass;~H or proton mass); NIST source
#kB=1.38064852E-23 #J/K (Boltzmann constant ); NIST source
#h=6.626070040E-34 #Js (Planck constant ); NIST source
#G=6.67408E-11  #kg (Newtonian constant of gravitation); NIST source   
hkB=4.79924466221135e-11 #=h/kB
pc=3.08567758E16  #m

#%%
def getfullimg_from_quad(img_quad, hN, N, input_Nwlen=0):
    """
    Compute full image NxN (4 quadrants) from the 1st quadrant image (hNxhN)
    Returns a NwlenxNxN image if input_Nwlen>0
    """

    if input_Nwlen != 0:
        img=np.empty((input_Nwlen, N, N), dtype=float)
        img[:, hN: , hN: ]=img_quad
        img[:, hN: , 0:hN]=img_quad[:,:,::-1] #left-right
        img[:, 0:hN, :   ]=img[:, N:hN-1:-1, :] #up-down
    else:
        img=np.empty((N, N), dtype=float)
        img[hN: , hN: ]=img_quad
        img[hN: , 0:hN]=img_quad[:,::-1] #left-right
        img[0:hN, :   ]=img[N:hN-1:-1, :] #up-down
        
    return img
#%%    

#%%
def BBsphere(wlen, Nwlen, Nstar, varpis, Teq, Req, **kwargs):
    """
    Compute intensity map (NstarXNstar) of spherical star with 
    radius Rs and Teff as a function of wavelength.
    The spectific intensity is given by a blackbody with T=Teff
    
    Input parameters
    Nstar: nb of pixels of intensity map (NstarXNstar)
    Nwlen: nb of wavelength values
    wlen : wavelength in meters    
    Req   : stellar radius in meters
    Teq : effective temperature in K
    """
    
    Istar=np.zeros((Nwlen, Nstar, Nstar), dtype=float)
    Istar[:, varpis <= Req]=BBwlen(wlen, Teq).reshape(Nwlen,1)

    return Istar
#%%


#%%
def gettauzdisc(wlen, Td, nu, rho0, ionfrac, powrho, powHd, \
             Xab, Req, Mstar, varpiReq, idxd, Nwlen, Np):
    """
    Compute vertical disc optical depth  excluding the 
    region occupied by the star.
    
    Vieira et al. 2015 (Eq. 4)
    
    kB=1.38064852E-23 #J/K (Boltzmann constant ); NIST source
    h=6.626070040E-34 #Js (Planck constant ); NIST source
    mH=1.660539040E-27  #kg (atomic mass;~H or proton mass); NIST source
    G=6.67408E-11  #kg (Newtonian constant of gravitation); NIST source
    """
    
    #mH=1.660539040E-27  #kg (atomic mass;~H or proton mass); NIST source    
    #hkB=4.79924466221135e-11 #=h/kB
    tauz_c1=0.03692  #=4(e2/4pieps0)^3/3mH2hc*sqrt(2pimH/3kB)
    mu=0.5 #Xab  #mean molecular weight => REPLACE BY CORRECT FORMULATION FOR SOUND SPEED
    tauz_c2=19783186.850612704*np.sqrt(Req**3/(mu*Mstar)) #=np.sqrt(np.pi*kB*Rstar**3/(mu*mH*G*Mstar))

    tauz=np.zeros((Nwlen, Np, Np), dtype=float)   
    """
    Better to write in the second form to ensure good broadcasting
    tauz[:, idxd[0], idxd[1]]=tauz_c1*(ionfrac*rho0/mH)**2 / (nu**3).reshape(Nwlen,1) * \
    (avggffdisc(wlen,Td[idxd[0],idxd[1]])+totgbfdisc(wlen,Td[idxd[0],idxd[1]]))* \
    (1.-np.exp(-hkB*nu.reshape(Nwlen,1)/Td[idxd[0],idxd[1]])) * \
    tauz_c2 * varpiReq[idxd[0], idxd[1]**(2*powrho+powHd)]
    """
    tauz[:, idxd[0], idxd[1]]= (tauz_c1*tauz_c2 * (ionfrac*rho0/mH)**2) *  \
    (avggffdisc(wlen,Td[idxd[0],idxd[1]])+totgbfdisc(wlen,Td[idxd[0],idxd[1]]))* \
    (1.-np.exp(-hkB*nu.reshape(Nwlen,1)/Td[idxd[0],idxd[1]])) * \
    (varpiReq[idxd[0], idxd[1]]**(2*powrho+powHd) / (nu**3).reshape(Nwlen,1))
    
    return tauz
#%%   

#%%
def avggffdisc(wlen, Td):
    """
    Free-free Gaunt factors
    Vieira2015 Table A1: Gaunt factor fitted parameters
    Coefficients Gn are computed for logT=[3.70, 3.82, 3.94, 4.06, 4.18, 4.30]
    """
    
    Nwlen=len(wlen)
    
    logTd=np.log10(Td)
    #Force values at the limits of temperature from Vieira15
    logTd[logTd < 3.70 ] = 3.70
    logTd[logTd > 4.30 ] = 4.30

    #Array of temperatures chosen by Vieira15     
    logT=np.array([3.70, 3.82, 3.94, 4.06, 4.18, 4.30])
    
    #Coefficients for free-free Gaunt factors
    G0=np.array([0.0952 ,  0.1001,  0.1097,  0.1250,  0.1470,  0.1761])
    G1=np.array([0.0215 ,  0.0421,  0.0639,  0.0858,  0.1071,  0.1269])
    G2=np.array([0.0145 ,  0.0130,  0.0111,  0.0090,  0.0068,  0.0046])

    #Convert wlen to microns and compute natural logarithm to match coefficients
    log_wlen_micron=np.log(wlen*1.E6)

    #Eq. A1 (Vieira2015)
    #Reshape for broadcasting correctly
    g_vieira=np.exp(G0.reshape(6,1)+
                    G1.reshape(6,1)*log_wlen_micron+
                    G2.reshape(6,1)*log_wlen_micron**2.)
      
    avggff=np.empty((Nwlen, len(logTd)), dtype=float)
    for iwlen in range(0, Nwlen):
        avggff[iwlen,:]=np.interp(logTd, logT, g_vieira[:,iwlen])
    
    #Simplified approximation
    #avggff=1.

    return avggff
#%% 

#%%
def totgbfdisc(wlen, Td):
    """
    Bound-free Gaunt factors
    Vieira2015 Table A1: Gaunt factor fitted parameters
    Coefficients Bn are computed for logT=[3.70, 3.82, 3.94, 4.06, 4.18, 4.30]
    """

    Nwlen=len(wlen)
    
    logTd=np.log10(Td)
    #Force values at the limits of temperature from Vieira15
    logTd[logTd < 3.70 ] = 3.70
    logTd[logTd > 4.30 ] = 4.30

    #Array of temperatures chosen by Vieira15     
    logT=np.array([3.70, 3.82, 3.94, 4.06, 4.18, 4.30])
   
    #Coefficients for bound-free Gaunt factors
    B0=np.array([ 2.2125,  1.6304,  1.1316,  0.6927,  0.2964, -0.0690])
    B1=np.array([-1.5290, -1.3884, -1.2866, -1.2128, -1.1585, -1.1185])
    B2=np.array([ 0.0563,  0.0413,  0.0305,  0.0226,  0.0169,  0.0126])

    #Convert wlen to microns and compute natural logarithm to match coefficients
    log_wlen_micron=np.log(wlen*1.E6)
    
    #Eq. A2 (Vieira2015)
    #Reshape for broadcasting correctly
    b_vieira=np.exp(B0.reshape(6,1)+
                    B1.reshape(6,1)*log_wlen_micron+
                    B2.reshape(6,1)*log_wlen_micron**2.)

    totgbf=np.empty((Nwlen, len(logTd)), dtype=float)
    for iwlen in range(0, Nwlen):
        totgbf[iwlen,:]=np.interp(logTd, logT, b_vieira[:,iwlen])
 
    #Simplified approximation    
    #totgbf=1.

    return totgbf
#%% 

#%%
def getIdisc(wlen, Bd, taui, idxd, Nwlen):
    """
    Compute disc specific intensity
    (Vieira15, Eq. 11 for A_disc, i.e. without star)
    """

    Id=np.zeros_like(taui) #faster then Bd=np.zeros((Nwlen, Np, Np), dtype=float)
    Id[:, idxd[0], idxd[1]]=Bd[:, idxd[0], idxd[1]] * \
        (1. - np.exp(-taui[:, idxd[0], idxd[1]]))

    return Id
#%%  

#%%
def getBdisc(wlen, Td, tauz, idxd, Nwlen): 
    """
    Compute disc blackbody emission excluding the region occupied by the star
    """
    
    Bd=np.zeros_like(tauz) #faster then Bd=np.zeros((Nwlen, Np, Np), dtype=float)
    Bd[:, idxd[0], idxd[1]]=BBwlen(wlen.reshape(Nwlen,1), Td[idxd[0], idxd[1]])
    
    return Bd
#%%  
 
 
#%%
def getTdisc(Td0, powTd, varpiReq, idxd): 
    """
    Compute disc temperature excluding the region occupied by the star
    """
    
    Td=np.zeros_like(varpiReq)
    Td[idxd[0], idxd[1]]=Td0 * varpiReq[idxd[0], idxd[1]]**powTd
    
    return Td
#%% 
 
#%%
def BBwlen(wlen, T):
    """
    Planck's law (blackbody radiation).
    Units of BBwlen are J/s/m/sr/m2 (SI)
    https://en.wikipedia.org/wiki/Planck's_law
    BBwlen=2 h c**2 / wlen**5  /(exp(h c/(wlen kT)-1)    
    Numerical constants used (from astropy.constants):
    2hc^2=2.*cte.h.si.value*cte.c.si.value**2=1.1910428681415875e-16
    hc/k=cte.h.value*cte.c.value/cte.k_B.value=0.014387769599838155
    
    Input parameters:
    wlen : wavelength in meters
    T    : temperature in K

    Test:
    All expressions below give 67204552 (J/s/m/sr/m2):
    blackbody_lambda(1.E-6 * u.m, 1000 * u.K).to(u.J/u.m/u.s/u.sr/u.m/u.m)
    blackbody_lambda(1.E-6 * u.m, 1000 * u.K)*1.E-7/1.E-10/1.E-4
    disco.BBwlen(1.E-6, 1000.)
    """

    return 1.1910428681415875e-16/wlen**5 / \
           (np.exp(0.014387769599838155/(wlen*T))-1.)
#%%     
    

    
    
#%%
# def kappa():
    
#     return klamb
#%%       
    
#%%
def H0(k_B, Td, mu, mH, G, Mstar, Req):
    """
    Compute the scale height
    """
    
    #Put physical constants k_B, mu, mH, G here ?
    
    return np.sqrt(k_B*Td*Req**3/(mu*mH*G*Mstar))
#%%

    
star_models={
    "BBsphere" : BBsphere
    #"BBellipse": BBellipse
}
    
#%%
def model(input_wlen=np.array([5.], dtype=float), wlen_unit='micron', 
          N=512,                  #Nb of pixels in image
          Rsun=695508000.0, #R_sun.si.value #m
          Msun=1.9891e+30,  #M_sun.si.value #kg
          Xab=1.,
          #Central star parameters
          star={"model":"BBsphere",
                "Req_sun"  :4.,  #Rsun
                "Teq"  :20000.0, #K
                "Mstar_sun":5., #Msun
                },
          #Disc parameters
          Rd=16.0, #Rsun
          Td0=12000., #K
          incl=0., #deg
          powTd=-0.75, #disc temperature power
          rho0 = 1.0E-8,  #kg/m3
          powrho=-3.5, #disc density power
          powHd=1.5,  #disc flaring power
          ionfrac=1.,
          maps=False
          ):

    """
    wlen=np.array([1., 5., 10], dtype=float)
    wlen_unit='micron' 
    N=512                  #Nb of pixels in image
    Rsun=695508000.0 #R_sun.si.value #m
    Msun=1.9891e+30  #M_sun.si.value #kg
    #Central star parameters
    Req=4.  #Rsun
    Tstar=20000.0 #K
    Mstar=5. #Msun
    #Disc parameters
    Rd=16.0 #Rsun
    Td0=5000. #K
    incl=60. #deg
    powTd=0.75
    """

    #Central star parameters (convert to MKS)
    #Rsun=R_sun.si.value  #m
    #Req*=Rsun #convert from Rsun to m
    #Mstar*=Msun #convert from Msun to kg
    star["Req"  ]=star["Req_sun"  ]*Rsun #convert from Rsun to m
    star["Mstar"]=star["Mstar_sun"]*Msun #convert from Msun to kg
    
    #Gas disc parameters (convert to MKS)
    Rd*=Rsun #convert from Rsun to m   
    
    #Numbers of pixels
    #N and Nstar have to be an even. 
    N = N if ((N % 2) == 0) else N+1 #force N to be even by adding 1 if necessary.
    #OR N+=(N % 2) #force N to be even    
    hN=int(N/2)
    Nstar=int(np.ceil(star["Req"]/Rd*N))  #highest integer (to span whole star)
    Nstar+=(Nstar % 2) #force N to be even
    hNstar=int(Nstar/2)

    #Wavelength
    wlen=np.copy(input_wlen)
    Nwlen=len(wlen)
    if (wlen_unit.lower() != 'm'):
        #Dictionary to convert to meters
        wlen_convert={ 
            #'m'       : 1.0,
            'mm'      : 1.0E-3,                  
            'micron'  : 1.0E-6,
            'nm'      : 1.0E-9,                   
            'angstrom': 1.0E-10
            }
        wlen*=wlen_convert[wlen_unit.lower()]  #force lower case to match dict
        
    #Compute frequency nu associated to wlen
    nu=cLight/wlen

    #Inclination
    cosi=np.cos(np.deg2rad(incl))

    #print('here')
    #print(Req/Rsun, Mstar/Msun, Rd/Rsun, rho0, Tstar, Td0, wlen)

    #Pole-on disc coordinates

    #(x',y') cartesian plane (pole-on disc)
    xp = np.linspace(-Rd, Rd, N, dtype=float)
    yp=xp/cosi
    
    #2D meshgrid for pole-on disc coordinates
    xxp, yyp = np.meshgrid(xp, yp)
    
    #2D meshgrid for pole-on disc in positive quadrant
    xxpquad, yypquad = np.meshgrid(xp[hN:], yp[hN:])

    #Distance to the central axis (pole-on disc)
    varpipquad=np.sqrt(xxpquad**2.+yypquad**2.) #circles for a given varpi constant
    #faster than np.hypot(xxpquad,yypquad)
    varpipquadReq=varpipquad/star["Req"]

    #indexes of pole-on disc excluding star (positive quadrant)
    idxd_quad=np.where(varpipquad > star["Req"])  
    
    
    #Coordinates of image plane xy (whole image)
    #Note that y=yp*cosi
    x = y = xp
    
    #2D meshgrid for pole-on disc coordinates (whole image)
    xx, yy = np.meshgrid(x, y)
    xxquad, yyquad = np.meshgrid(x[hN:], y[hN:]) #Positive quadrant    
    
    
    #---------------- Disc ----------------
    #Disc temperature in positive quadrant
    Td_quad=getTdisc(Td0, powTd, varpipquadReq, idxd_quad)
    

    #Disc vertical optical depth in positive quadrant 
    #(Vieira15, Eq. 4)
    """
    tauz_quad=0.03692(rho0/mH*ionfrac)**2 * nu**3* \
                (avggffdisc(wlen, Td_quad)+totgbfdisc(wlen, Td_quad))* \
                (1.-np.exp(h*nu/kB*))
    """
    
    tauz_quad=gettauzdisc(wlen, Td_quad, nu, rho0, ionfrac, powrho, powHd, \
            Xab, star["Req"], star["Mstar"], varpipquadReq, idxd_quad, Nwlen, hN)
    
    
    #Optical depth for disc seen at angle i in positive quadrant
    #(Vieira15, Eq. 6)
    #taui=tau0/cosi*(varpi/Req)**taupow
    taui_quad=tauz_quad/cosi
    
     
    #Blackbody specific intensity in the disc at temperature Td    
    Bdisc_quad=getBdisc(wlen, Td_quad, tauz_quad, idxd_quad, Nwlen)

  
    #Disc specific intensity
    #(Vieira15, Eq. 11 for A_disc, i.e. without star)
    #Positive quadrant
    Idisc_quad=getIdisc(wlen, Bdisc_quad, taui_quad, idxd_quad, Nwlen)
    #Whole image (4 quadrants)
    Idisc=getfullimg_from_quad(Idisc_quad, hN, N, input_Nwlen=Nwlen)
    #-----------------end disc------------------------------


    """
    CREATE ROUTINE TO COMPUTE ALL xs, xxs, xx, varpi, etc AND STORE IN DICT
    DOES NOT RECOMPUTE IF DICT IS GIVEN AS INPUT
    MAYBE NOT A GOOD IDEA SINCE IT DEPENDS ON i, Req, Rd, etc...
    """
    
    #---------------- Central star ----------------    
    #Can use symmetries if useful
    #check /Users/domicia/informatica/imaging_graphics/python/my_tests/analytical_models.py
    
    #Using slices to make sure that xs and ys match the x and y arrays for the whole image
    
    xs=ys=x[hN-hNstar:hN+hNstar]
    #or xs=ys=x[(x <= Req) & (x >= -Req)] #much slower
    
    xxs = xx[hN-hNstar:hN+hNstar, hN-hNstar:hN+hNstar]
    yys = yy[hN-hNstar:hN+hNstar, hN-hNstar:hN+hNstar]
    #or xxs, yys = np.meshgrid(xs, ys) #slower
    
    varpis=np.sqrt(xxs**2.+yys**2.)    
    
    #Istar=np.zeros((Nwlen, Nstar, Nstar), dtype=float)    
    Istar=star_models[star["model"]](wlen, Nwlen, Nstar, varpis, **star)
    #-----------------end central star----------------------


    #---------------- Disc+Central star ----------------    
    #Upper half part of disc+star and star regions
    Imap=np.copy(Idisc)
    upIdisc_star=Imap[:,hN:hN+hNstar,hN-hNstar:hN+hNstar] #already has the disc
    upIstar=Istar[:,hNstar:,:]
    
    #Upper inner region without disc contribution (star alone or in front of disc)
    #Vieira et al. 2015, Eq. 11 (A*^1/2)
    idx=upIstar > 0.0
    upIdisc_star[idx]=upIstar[idx]

    #Since upIdisc_star is the same as Imap, Imap already has the disc+star contribution
    #Imap[:,hN:hN+hNstar,hN-hNstar:hN+hNstar]=upIdisc_star
    
    #Lower part of disc+star and star regions
    lowIdisc_star=Imap[:,hN-hNstar:hN,hN-hNstar:hN+hNstar] #already has the disc contribution
    lowIstar=Istar[:,0:hNstar,:]  

    #Get taui from taui_quad for the lower disc+star part
    lowtauidisc_star=np.zeros((Nwlen, hNstar, Nstar), dtype=float)
    #Copy taui_quad using up-down symmetry
    lowtauidisc_star[:,:,hNstar:]=taui_quad[:,hNstar-1::-1,0:hNstar] #up-down
    #Copy the right half part of lowtauidisc_star to the left half part using symmetry
    lowtauidisc_star[:,:,0:hNstar]=lowtauidisc_star[:,:,Nstar-1:hNstar-1:-1] #left-right

    #Intersection disc-star (region defined by rays crossing the disc until the star)
    #Vieira et al. 2015, Eq. 11 (A*^-1/2)
    idx=((lowIdisc_star > 0.0) & (lowIstar > 0.0))
    lowIdisc_star[idx]+=lowIstar[idx]*np.exp(-lowtauidisc_star[idx])
    #same as lowIdisc_star[idx]=lowIdisc_star[idx]+lowIstar[idx]*np.exp(-lowtauidisc_star[idx])

    #Lower inner region without disc contribution (only the star)
    idx=(lowIdisc_star == 0.0)
    lowIdisc_star[idx]=lowIstar[idx]
 
    #Since lowIdisc_star is the same as Imap, Imap already has the disc+star contribution
    #Imap[:,hN-hNstar:hN,hN-hNstar:hN+hNstar]=lowIdisc_star

    #-----------------end Disc+Central star--------------------


    """
    If maps is given as a dict to disco then return taui, Idisc, Istar, etc in maps
    If necessary images are created from 1st quadrant using getfullimg_from_quad
    Also return the total flux (SED) of each image (see tests in disco_for_amhra_sv1.py).
    """    
    if type(maps) == dict:

        
        maps["wlen"]=wlen

        maps["N"]=N
        maps["x"]=x
        maps["y"]=y
    
        maps["Nstar"]=Nstar
        maps["xs"   ]=xs
        maps["ys"   ]=ys        

        maps["taui" ]=getfullimg_from_quad(taui_quad, hN, N, input_Nwlen=Nwlen)
        maps["Tdisc"]=getfullimg_from_quad(Td_quad, hN, N)
        maps["Idisc"]=Idisc
        maps["Istar"]=Istar
        maps["Imap" ]=Imap


        dOmega=(x[1]-x[0])*(y[1]-y[0])/(10.0*pc)**2.  #rad^2
        #maps["Flux_10pc"]=np.sum(maps['Imap']*dOmega, axis=(1,2))*1.0e-6 #(J/s/micron/m2)
        maps["Flux_10pc"]=np.sum(maps['Imap']*dOmega, axis=(1,2)) #J/s/m/m2 (SI)
    
    return Imap, xx, yy
    
#%%   


  