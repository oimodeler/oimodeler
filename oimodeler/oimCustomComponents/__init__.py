# -*- coding: utf-8 -*-
"""Custom model components added by the community

If you add python files with custom components in the oimCustomComponents
directory add the import here so that all component will be available through
the customComponent module
"""
# TODO: Find a way to automatically add all class that derives from oimComponent
# from all files in the oimCustomComponents directory

from .oimAsymRing import oimAEIRing, oimAERing
from .oimBox import oimBox
from .oimExpRing import oimExpRing
from .oimFastRotator import oimFastRotator
from .oimGaussLorentz import oimGaussLorentz
from .oimKinematicDisk import oimKinematicDisk
from .oimRadialRing import oimRadialRing
from .oimSpiral import oimSpiral
from .oimStarHaloDisk import oimStarHaloGaussLorentz, oimStarHaloIRing, oimStarHaloRing
from .oimTempGradient import oimTempGradient
from .oimKinematicDisk import oimKinematicDisk
