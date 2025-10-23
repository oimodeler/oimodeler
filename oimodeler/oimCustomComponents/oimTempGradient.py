import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import constants as const
from ..oimOptions import oimOptions
from ..oimParam import oimParam
from ..oimUtils import blackbody, linear_to_angular


class oimTempGrad(oimComponentRadialProfile):
    """A ring defined by a radial temperature profile in r^q and a radial dust
    surface density profile in r^p.

    Parameters
    ----------
    rin : float
        Inner radius of the disk [au].
    rout : float
        Outer radius of the disk [au].
    r0 : float
        Reference radius [au].        
    T0 : float
        Temperature at reference radius r0 [K].
    Sigma0 : float
         Dust surface density at reference radius r0 [M_sun].
    Mdust : float
         Mass of the dusty disk [M_sun].
    q : float
        Power-law exponent for the temperature profile.
    p : float
        Power-law exponent for the dust surface density profile.
    
    dist : float
        Distance of the star [pc].

    Attributes
    ----------
    params : dict with keys of str and values of oimParam
        Dictionary of parameters.
    _r : array_like
    _wl : array_like
        Wavelengths [micron].
    _t : array_like
        Times [second].

    Methods
    -------
    _radialProfileFunction(r, wl, t)
        Calculates a radial temperature gradient profile via a dust-surface
        density- and temperature profile.
    """

    name = "Temperature Gradient"
    shortname = "TempGrad"
    elliptic = True

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(
            name="rin",
            value=0,
            unit=u.au,
            description="Inner radius of the disk",
        )
        self.params["rout"] = oimParam(
            name="rout",
            value=0,
            unit=u.au,
            description="Outer radius of the disk",
        )
        self.params["r0"] = oimParam(
            name="r0",
            value=1,
            unit=u.au,
            free=False,
            description="Reference radius",
        )
        self.params["temp0"] = oimParam(
            name="temp0",
            value=300,
            unit=u.K,
            description="Temperature at reference radius",
        )
        self.params["sigma0"] = oimParam(
            name="sigma0",
            value=0,
            unit=u.g / u.cm**2,
            description="Dust surface density at reference radius",
        )
        self.params["dust_mass"] = oimParam(
            name="dust_mass",
            value=0,
            unit=u.M_sun,
            description="Mass of the dusty disk",
        )
        self.params["q"] = oimParam(
            name="q",
            value=-0.5,
            unit=u.one,
            description="Power-law exponent for the temperature",
        )
        self.params["p"] = oimParam(
            name="p",
            value=-1,
            unit=u.one,
            description="Power-law exponent for the dust surface density",
        )
        self.params["dist"] = oimParam(
            name="dist",
            value=0,
            unit=u.pc,
            free=False,
            description="Distance of the star",
        )

        self.params["f"].free = False #Flux contribution parameter of the TG disk (not fitted)
        self.params["f"].value = 1.0  #Set to 1.0 since the temperature gradient model returns itself a physical flux in Jy

        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static
        self._wl_kappa_abs = None #Wavelength array of the opacity file (in m)
        self._kappa_abs = None #opacity values of the opacity file (mass absorption coefficient in cm**2/g)
        self._eval(**kwargs)

    def _radialProfileFunction(
        self, r: np.ndarray, wl: np.ndarray, t: np.ndarray
    ) -> np.ndarray:
        """Calculates a radial temperature gradient profile via a dust-surface
        density- and temperature profile.

        Parameters
        ----------
        r : numpy.ndarray
            Radial grid [mas].
        wl : numpy.ndarray
            Wavelengths [m].
        t : numpy.ndarray
            Times [second].

        Results
        -------
        radial_profile : numpy.ndarray
        """
        # HACK: Sets the multi wavelength coordinates properly.
        # Does not account for time, improves computation time.
        # TODO: Check if this is still needed and if not remove it

        #Getting the wavelength array from the data
        wl = np.unique(wl)
        
        #Getting the distance to the star
        dist = self.params["dist"].value
        
        #Getting the opacity table provided by the user
        kappa_abs = self.get_kappa(wl)
        
        #Getting the disk elongation
        elong=self.params["elong"].value
        
        #Adding the correct dimensions to the radial grid, wavelength array and opacity table 
        if len(r.shape) == 3:
            r = r[0, 0][np.newaxis, np.newaxis, :]
            wl, kappa_abs = map(
                lambda x: x[np.newaxis, :, np.newaxis], [wl, kappa_abs]
            )
        else:
            wl, kappa_abs = map(lambda x: x[:, np.newaxis], [wl, kappa_abs])
            r = r[np.newaxis, :]
            
        #Getting and converting the different radii in cm for the computation of the surface density profile 
        rin, rout, r0 = map(lambda x: self.params[x](wl, t), ["rin", "rout", "r0"])
        rin_cm, rout_cm, r0_cm = map(
            lambda x: x * self.params["rin"].unit.to(u.cm), [rin, rout, r0]
        )
        
        #Getting the exponents of the temperature and surface density power laws
        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        
        #Converting the reference radius into mas
        r0 = linear_to_angular(r0, dist) * 1e3 
        
        # Radial temperature profile computation
        temp = self.params["temp0"](wl, t) * (r / r0) ** (q) 
        # Radial surface density profile computation
        sigma0 = self.params["sigma0"](wl, t)
        if sigma0 == 0: #Computation of sigma0 and then the surface density profile sigma using dust_mass
            dust_mass = self.params["dust_mass"](wl, t) * self.params[
                "dust_mass"
            ].unit.to(u.g)
            if p == -2:
                sigma0 = dust_mass / (
                    2.0 * np.pi * np.log(rout_cm / rin_cm) * r0_cm**2
                )
            else:
                f = (rout_cm**(2 + p) - rin_cm**(2 + p)) / (2 + p)
                sigma0 = dust_mass / (2.0 * np.pi * f * r0_cm**(-p)) 
            sigma = sigma0 * (r / r0) ** (p)
        else: #Computation of the surface density profile sigma directly from sigma0
            sigma = sigma0 * (r / r0) ** (p)

        #Computation of the emissivity factor from sigma and kappa_abs
        epsilon = 1 - np.exp(-sigma * kappa_abs * elong)
        
        #Computation of the radial intensity profile denoted as 'image' hereafter
        nu=const.c/wl
        spectral_density = blackbody(temp, nu) * epsilon * (1./elong)
        rin_mas, rout_mas = map(
            lambda x: linear_to_angular(x, dist) * 1e3, [rin, rout]
        )
        radial_profile = ((r > rin_mas) & (r < rout_mas)).astype(int)
        image = np.nan_to_num(radial_profile * spectral_density, nan=0)
        if len(r.shape) == 3:
            return image

        return image

    @property
    def _r(self):
        """Gets the radial profile (mas)."""
        rin = (
            linear_to_angular(
                self.params["rin"].value, self.params["dist"].value
            )
            * 1e3
        )
        rout = (
            linear_to_angular(
                self.params["rout"].value, self.params["dist"].value
            )
            * 1e3
        )
        if oimOptions.model.grid.type == "linear":
            return np.linspace(rin, rout, self.params["dim"].value)
        return np.logspace(
            0.0 if rin == 0 else np.log10(rin),
            np.log10(rout),
            self.params["dim"].value,
        )

    @_r.setter
    def _r(self, value):
        """Sets the radial profile [mas]."""
        return
    
    def get_kappa(self, wl: np.ndarray):
        wl_unique=np.unique(wl)
        if np.array_equal(wl_unique,self._wl_kappa_abs):
            return self._kappa_abs
        else:
            kappa_abs_interp=np.interp(wl_unique,self._wl_kappa_abs,self._kappa_abs)
            return kappa_abs_interp
            
    def setKappaAbs(self, kappa_file,um=1):
        kappa_data=np.loadtxt(kappa_file,usecols=(0,1))
        if um==1:
            self._wl_kappa_abs=kappa_data[:,0]*1e-6
        else:
            self._wl_kappa_abs=kappa_data[:,0]
        self._kappa_abs=kappa_data[:,1]


class oimAsymTempGrad(oimTempGrad):
    """A ring defined by a radial temperature profile in r^q and a radial dust
    surface density profile in r^p.

    Parameters
    ----------
    rin : float
        Inner radius of the disk [au].
    rout : float
        Outer radius of the disk [au].
    Tin : float
        Inner radius temperature [K].
    Mdust : float
        Mass of the dusty disk [M_sun].
    q : float
        Power-law exponent for the temperature profile.
    p : float
        Power-law exponent for the dust surface density profile.
    kappa_abs : float or oimInterp
        Dust mass absorption coefficient [cm2.g-1].
    dist : float
        Distance of the star [pc].

    Attributes
    ----------
    params : dict with keys of str and values of oimParam
        Dictionary of parameters.
    _r : array_like
    _wl : array_like
        Wavelengths [micron].
    _t : array_like
        Times [second].

    Methods
    -------
    _radialProfileFunction(r, wl, t)
        Calculates a radial temperature gradient profile via a dust-surface
        density- and temperature profile.
    """

    name = "Asymmetric Temperature Gradient"
    shortname = "AsymTempGrad"
    elliptic = True
    asymmetric = True
