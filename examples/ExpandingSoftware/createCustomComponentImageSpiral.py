# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as units

import oimodeler as oim

path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


class oimSpiral(oim.oimComponentImage):
    name = "Spiral component"
    shorname = "Sp"

    # NOTE: Setting elliptic to True to use the elong parameter/keword for a change of variable
    # to "elongate" objects
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # NOTE: Component parameters. Note that as it inherits from the oimComponentImage class,
        # it already has x, y, f and dim as parameters
        self.params["fwhm"] = oim.oimParam(**oim._standardParameters["fwhm"])
        self.params["P"] = oim.oimParam(
            name="P", value=1, description="Period in mas", unit=units.mas
        )
        self.params["width"] = oim.oimParam(
            name="width",
            value=0.01,
            description="Width as filling factor",
            unit=units.one,
        )

        self._pixSize = 0.05 * units.mas.to(units.rad)

        self._t = np.array([0])  # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model

        # NOTE: Finally, evalutating paramters in the same way as for all other components
        self._eval(**kwargs)

    def _imageFunction(self, xx, yy, wl, t):
        # NOTE: As xx and yy are transformed coordinates, r and phi takes into account
        # the ellipticity and orientation using the pa and elong keywords
        r = np.sqrt(xx**2 + yy**2)
        phi = np.arctan2(yy, xx)

        p = self.params["P"](wl, t)
        sig = self.params["fwhm"](wl, t) / 2.35
        w = self.params["width"](wl, t)

        im = 1 + np.cos(-phi - 2 * np.pi * np.log(r / p + 1))
        im = (im < 2 * w) * np.exp(-(r**2) / (2 * sig**2))
        return im


# NOTE: Create a model including the spiral component
ud = oim.oimUD(d=1, f=0.2)
c = oimSpiral(dim=256, fwhm=5, P=0.1, width=0.2, pa=30, elong=2, x=10, f=0.8)
m = oim.oimModel(c, ud)

# NOTE: Plot images of the model
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
m.showModel(
    256, 0.1, swapAxes=True, fromFT=False, normPow=1, axe=ax[0], colorbar=False
)
m.showModel(
    256, 0.1, swapAxes=True, fromFT=True, normPow=1, axe=ax[1], colorbar=False
)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")
fig.savefig(save_dir / "customCompImageSpiral.png")

# NOTE: Compute and plot the visibilities for various baselines and wavelengths
nB = 5000
nwl = 1
wl = 0.5e-6

B = np.linspace(0, 100, num=nB // 2)
Bx = np.append(B, B * 0)
By = np.append(B * 0, B)

spfx = Bx / wl
spfy = By / wl

vc = m.getComplexCoherentFlux(spfx, spfy)
v = np.abs(vc / vc[0])

fig, ax = plt.subplots(1, 1)
label = [
    "East-West Baselines",
]

ax.plot(
    B / wl / units.rad.to(units.mas),
    v[: nB // 2],
    color="r",
    label="East-West Baselines",
)
ax.plot(
    B / wl / units.rad.to(units.mas),
    v[nB // 2 :],
    color="b",
    label="North-South Baselines",
)

ax.set_xlabel(r"B/$\lambda$ (cycles/mas)")
ax.set_ylabel("Visibility")
ax.legend()
fig.savefig(save_dir / "customCompImageSpiralVis.png")
