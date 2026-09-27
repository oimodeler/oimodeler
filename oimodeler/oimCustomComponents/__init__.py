# -*- coding: utf-8 -*-
"""Custom model components added by the community

If you add python files with custom components in the oimCustomComponents
directory add the import here so that all component will be available through
the customComponent module
"""

# TODO: Find a way to automatically add all class that derives from oimComponent
# from all files in the oimCustomComponents directory

from .oimAsymRing import oimAEIRing
from .oimBinaryOrbit import oimBinaryOrbit
from .oimBipolar import oimBipolar
from .oimBox import oimBox
from .oimDisco import oimDisco
from .oimExpRing import oimExpRing

from .oimRadialRings import (
    oimRadialExpRing,
    oimRadialPowRing,
    oimRadialPowRing2)

from .oimFastRotator import (
    oimFastRotator,
    oimFastRotatorLLDD,
    oimFastRotatorMasse,
    oimFastRotatorNLLDD,
    oimFastRotatorQuadLDD,
)
from .oimGaussLorentz import oimGaussLorentz
from .oimInnerRim import oimInnerRim
from .oimKinematicDisk import oimKinematicDisk
from .oimSpiral import oimSpiral
from .oimStarHaloDisc import oimStarHaloGaussLorentz, oimStarHaloIRing
from .oimTempGrad import oimTempGrad
