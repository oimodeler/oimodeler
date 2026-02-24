# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 14:42:21 2026

@author: ame
"""

from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import oimodeler as oim

plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True

path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# NOTE: Set the oimBinaryOrbit component
orb = oim.oimBinaryOrbit(
    e=0.4,  # Eccentricity
    a=10,  # semi-major axis (mas)
    T0=0,  # Time Periastron passage (MJD)
    T=1,  # Period (in days by default or an astropy unit if specified)
    i=45,  # inclination angle (deg)
    O=40,  # Longitude of ascending node (deg)
    o=-20,  # Argument of periastron
)

# NOTE: Get the separation at various time to plot the orbit
nt2 = 7
t = np.linspace(t0, t0 + T, 100)
t2 = np.linspace(t0, t0 + T, nt2)
x, y = orb.getSeparation(t, mas=True)
x2, y2 = orb.getSeparation(t2, mas=True)

# NOTE: Plot the full orbit
fig, ax = plt.subplots(figsize=(7, 7))
ax.grid(which="major", lw=1, alpha=0.5, zorder=-10)
ax.grid(which="minor", lw=0.5, alpha=0.5, zorder=-10)

ax.plot(x, y, color="lightgrey", alpha=1, lw=4)
ax.scatter(x2[:-1], y2[:-1], color="r", marker=".", zorder=10)

for i in range(nt2 - 1):
    ax.annotate(
        f"t={t2[i]:.2f}",
        (x2[i] + 0.1, y2[i] + 0.2),
        color="r",
        zorder=10,
        fontsize=8,
    )
ax.axis("equal")

ax.scatter(0, 0, marker="o", color="k")
ax.set_xlabel("x (mas)")
ax.set_ylabel("y (mas)")
ax.set_xlim(10, -15)
plt.savefig(save_dir / "ExampleBinary_orbit_plot.png")

# NOTE: Baselines at different times
N = 200
wl = 2.1e-6
B = np.linspace(0.0, 100, num=N)
spf = B / wl

orb.primary = oim.oimUD(d=3)
orb.secondary.params["f"].value = 0.3

morb = oim.oimModel(orb)
fig, ax = plt.subplots(
    3, nt2, figsize=(nt2 * 2.5, 8), sharex="row", sharey="row"
)

mascycle = 1.0 / u.rad.to(u.mas)

for iT in range(nt2):
    for iPA in range(2):
        v = np.abs(
            morb.getComplexCoherentFlux(
                spf * (iPA == 1), spf * (iPA == 0), t=np.ones(N) * t2[iT]
            )
        )
        v = v / np.max(v)
        ax[iPA, iT].plot(spf * mascycle, v, alpha=0.5)
        ax[iPA, iT].set_ylim(0, 1)

    ax[0, iT].text(
        0.5,
        1.01,
        f"sep = ({x2[iT]:.1f},{y2[iT]:.1f}) mas",
        transform=ax[0, iT].transAxes,
        va="bottom",
        ha="center",
    )
    ax[0, iT].text(
        0.5,
        1.10,
        f"t={t2[iT]:.2f} day",
        transform=ax[0, iT].transAxes,
        va="bottom",
        ha="center",
    )
    ax[1, iT].set_xlabel("B/$\\lambda$ (cycle/mas)")
    if iT != 0:
        ax[0, iT].axes.get_yaxis().set_visible(False)
        ax[1, iT].axes.get_yaxis().set_visible(False)

ax[0, 0].set_ylabel("East-West Baseline Vis.")
ax[1, 0].set_ylabel("North-South Baseline Vis.")
dim = 64
fov = 30
pix = fov / dim

for iT in range(nt2):
    morb.showModel(
        dim,
        pix,
        t=t2[iT],
        fromFT=True,
        axe=ax[2, iT],
        normPow=0.5,
        colorbar=False,
    )
    ax[2, iT].plot(x, y, color="r", alpha=0.2, lw=4, ls=":")
    if iT != 0:
        ax[2, iT].axes.get_yaxis().set_visible(False)

plt.tight_layout()
plt.savefig(save_dir / "ExampleBinary_visi_through_full_orbit.png")

# NOTE: Multi spatial frequencies and times for one (simultaneous) "getComplexCoherentFlux" call
nt = 500  # number of epochs in one period
t = np.linspace(t0, t0 + T, nt)
nB = 500  # number of baselines
wl = 2.1e-6  # fixed wavelength
B = np.linspace(0.0, 100, num=nB)
spf = B / wl

Bs = np.tile(B, (nt, 1))
ts = np.transpose(np.tile(t, (nB, 1)))
spf = Bs.flatten() / wl
spf0 = spf * 0

vis1 = np.abs(
    morb.getComplexCoherentFlux(spf, spf * 0, t=ts.flatten())
).reshape(len(t), len(B))
vis1 /= np.outer(np.max(vis1, axis=1), np.ones(nB))
vis2 = np.abs(
    morb.getComplexCoherentFlux(spf * 0, spf, t=ts.flatten())
).reshape(len(t), len(B))
vis2 /= np.outer(np.max(vis2, axis=1), np.ones(nB))

# NOTE: Plot the colormap of the temporal evolution of the visibility during one period
fig, ax = plt.subplots(1, 2, figsize=(15, 8), sharex=True, sharey=True)

ax[0].pcolormesh(Bs, ts, vis1, cmap="plasma")
sc = ax[1].pcolormesh(Bs, ts, vis2, cmap="plasma")

ax[0].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[0].set_ylabel("time (days)")
ax[0].set_title("East-West Baselines")
ax[1].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[1].set_ylabel("time (days)")
ax[1].set_title("North-South Baselines")
plt.tight_layout()
fig.colorbar(sc, ax=ax, label="Visiblity")
plt.savefig(save_dir / "ExampleBinary_temporal_visi_map.png")

# NOTE: Chromatic components with different temperatures
dist = 100

# NOTE: Use the "starWl"" interpolator to set a blackbody flux to each component.
orb.primary = oim.oimUD(
    d=1, f=oim.oimInterp("starWl", temp=3000, dist=dist, radius=0.5)
)
orb.secondary = oim.oimPt(
    f=oim.oimInterp("starWl", temp=30000, dist=dist, radius=0.1)
)
morb = oim.oimModel(orb)

nwl = 400
nspf = 2000

spf_res = 1.22 / (
    orb.primary.params["d"].value * u.mas.to(u.rad)
)  # first zero of the primary

spf = np.linspace(0, spf_res * 1.5, nspf)
wl = np.linspace(1e-6, 4e-6, num=nwl)
wl = np.logspace(-6, -5.4, num=nwl)
spf_2D = np.tile(spf, (nwl, 1)).flatten()
wl_2D = np.transpose(np.tile(wl, (nspf, 1))).flatten()

vis = np.abs(
    morb.getComplexCoherentFlux(spf_2D, spf_2D * 0, wl_2D, t=spf_2D * 0 + 0.1)
).reshape(nwl, nspf)
vis /= np.outer(np.max(vis, axis=1), np.ones(nspf))

fig, ax = plt.subplots(1, 1, figsize=(14, 8))
sc = ax.scatter(spf_2D, vis, c=wl_2D * 1e6, s=0.2, cmap="plasma")

ax.set_xlabel("B/$\\lambda$ (cycles/rad)")
ax.set_ylabel("Visiblity")
ax.margins(0, 0)
plt.tight_layout()
fig.colorbar(sc, ax=ax, label="$\\lambda$ ($\\mu$m)")
plt.savefig(save_dir / "ExampleBinary_chromatic_binary.png")

# NOTE: Plot the flux of both components
fp = orb.primary.params["f"](wl)
fs = orb.secondary.params["f"](wl)

fig, ax = plt.subplots()

ax.plot(
    wl * 1e6,
    fp,
    label=f"primary: T={orb.primary.params['f'].temp.value}K R={orb.primary.params['f'].radius.value}Ro",
)
ax.plot(
    wl * 1e6,
    fs,
    label=f"secondary: T={orb.secondary.params['f'].temp.value}K R={orb.secondary.params['f'].radius.value}Ro",
)
ax.legend()
ax.set_yscale("log")
ax.set_xlabel("$\\lambda$ ($\\mu$m)")
ax.set_ylabel("Flux (Jy)")
plt.savefig(save_dir / "ExampleBinary_chromatic_binary_fluxes.png")

# NOTE: Simulate radial velocity
orb.params["Ka"].value = 5
orb.params["Kb"].value = 50
orb.params["V0"].value = -5

t = np.linspace(0, 0 + 2, nt * 2)
rv = orb.getPrimaryRadialVelocity(t)

rv_a, rv_b = orb.getRadialVelocities(t)
plt.figure()
plt.plot(t, rv_a, label="RV_a")
plt.plot(t, rv_b, label="RV_b")
plt.xlabel("Time")
plt.ylabel("radial velocity (km/s)")
plt.legend()
plt.savefig(save_dir / "ExampleBinary_rv.png")
