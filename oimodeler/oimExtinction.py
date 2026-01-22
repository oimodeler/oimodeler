from pathlib import Path
import numpy as np
from scipy import interpolate

fitzindeb = np.genfromtxt(Path(__file__).parent / 'extlaws' / 'FitzIndeb_3.1_VOSA.dat', unpack=True)
fitzindebspline = interpolate.splrep(fitzindeb[0]/1e10, fitzindeb[1], s=1)

def extlaw_FitzIndeb(wavelength, A_V=10.0):

    kappa = interpolate.splev(wavelength, fitzindebspline, der=0)
    return A_V * (kappa / 211.4)