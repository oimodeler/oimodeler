# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:55:44 2026

@author: ame
"""
from pathlib import Path
import astropy.units as u
import oimodeler as oim
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling.models import BlackBody

path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    
#%% Create a DISCO model 
disc=oim.oimDisco(incl=70,Td0=15000,powTd=0,dim=256,rho0=2e-8,Rd=50)
m = oim.oimModel(disc)
disc.normalizeImage = False # if the absolute Fluxes are important

#%% plot images of the disk at different wavelengths
wl_im=np.array([0.1,0.5,1,2,4,8,12])*1e-6
pixmas = disc.getPixelSize(mas=True)

fig,ax=plt.subplots(1,len(wl_im),figsize=(15,2.5))
_,_,im=m.showModel(disc.params["dim"].value,pixmas,fromFT=False,wl=wl_im,cmap="hot",
            legend=True,normalize=True,normPow=1,axe=ax,colorbar=False)
fig.tight_layout()
plt.savefig(save_dir / "ExampleDISCO_images.png")
#%% Create a Uniform disk of the same angular size than the primary foer comparison
D=2*disc.params["Rstar"].quantity.to(u.m).value
d= disc.params["dist"].quantity.to(u.m).value
theta=D/d*u.rad.to(u.mas)
ud = oim.oimUD(d=theta)
m2=oim.oimModel(ud)

#%% Computing Chromatic visibilities and reference uniform disk
oim.oimOptions.ft.padding=8 # improve FFT precison (default =4)

nwl = 100
nspf = 1000

# first zero of the primary used to define the spatial-frequency range
spf_res = 1.22 / (ud.params["d"].value * u.mas.to(u.rad)) 
spf = np.linspace(0, spf_res * 1.5, nspf)
wl = np.logspace(-6.5, -5, num=nwl)
spf_2D = np.tile(spf, (nwl, 1)).flatten()
wl_2D = np.transpose(np.tile(wl, (nspf, 1))).flatten()

vis = np.abs(m.getComplexCoherentFlux(spf_2D, spf_2D*0, wl_2D)).reshape(nwl, nspf)
vis /= np.outer(np.max(vis, axis=1), np.ones(nspf))


wl_ref = 0.5e-6
B_ref = np.linspace(0.0, 200, num=200)
spf_ref = B_ref/wl_ref

vis0 = np.abs(m2.getComplexCoherentFlux(spf_ref, spf_ref*0))
vis0 /= np.max(vis0)
#%% Plotting chromatic visibilties and the reference uniform disk
fig, ax = plt.subplots(1, 1, figsize=(14, 8))
sc = ax.scatter(spf_2D, vis, c=wl_2D*1e6, s=0.2, cmap="plasma")
fig.colorbar(sc, ax=ax, label="$\\lambda$ ($\\mu$m)")
ax.set_xlabel("B/$\\lambda$ (cycles/rad)")
ax.set_ylabel("Visiblity")
ax.margins(0, 0)

ax.plot(spf_ref,vis0,color="r",ls=":",lw=5,label=f"Reference UD={ud.params["d"].value:.2f}mas")
ax.legend()

plt.savefig(save_dir / "ExampleDISCO_Vis.png")

#%%

# to compute accurate SED with minimum interpolation we need to moidify 
# the internal wavelength grid of the disc component
disc._wl=np.logspace(-7,-4.8,num=120)  
units=[u.W/u.m**3,u.Jy]
fig,ax = plt.subplots(1,2,figsize=(15,5))
for i,unit in enumerate(units):
    disc.outunit=unit
    wl = disc._wl
    im= m.getImage(disc.params["dim"].value,pixmas,fromFT=False,wl=wl)
    flx= np.sum(im,axis=(1,2))
    
    ax[i].loglog(wl*1e6,flx,label="oimDisco")
    ax[i].set_xlabel("$\\lambda$ ($\\mu$m)")
    ax[i].set_ylabel(f"Flux ({disc.outunit})")
    
    bb = BlackBody(disc.params["Teff"].quantity)
    flx_bb = (bb(wl*u.m)*np.pi*u.sr*(disc.params["Rstar"].quantity)**2/(disc.params["dist"].quantity)**2)
    flx_bb = flx_bb.to(disc.outunit,equivalencies=u.spectral_density(wl*u.m))
    
    ax[i].plot(wl*1e6,flx_bb,label="Blackbody")
    ax[i].set_title(f"DISCO Flux in {unit}")
    
ax[0].legend()   
plt.savefig(save_dir / "ExampleDISCO_flux.png")

#%%
T0 = 0
T  = 1
orb = oim.oimBinaryOrbit(e=0.4, # Eccentricity
                         a=7,  # semi-major axis (mas)
                         T0 = T0,     # Time Periastron passage (MJD by default or decimal year)
                         T  = T,      # Period (in days by default or any compatible astropy unit if specified)
                         i=45,   # inclination angle (deg)
                         O=40,   # Longitude of ascending node (deg)
                         o=-20    # Argument of periastron
                         )


orb.primary=disc
orb.secondary=oim.oimPt(f=oim.oimInterp("starWl", temp=12000,  radius=4))
orb.secondary.params["f"].dist = orb.primary.params["dist"]
morb = oim.oimModel(orb)
disc.outunit=u.Jy

#%%Plot the flux of both components
fp = flx
fs = orb.secondary.params["f"](wl)

fig, ax = plt.subplots()

ax.loglog( wl * 1e6, fp,label="primary")
ax.plot(wl * 1e6, fs, label="secondary")
ax.legend()
ax.set_yscale("log")
ax.set_xlabel("$\\lambda$ ($\\mu$m)")
ax.set_ylabel("Flux (Jy)")
plt.savefig(save_dir / "ExampleDISCO_orbit_flux.png")
#%%

fig, ax = plt.subplots()
ax.plot(wl * 1e6, fs/(fs+fp), label="secondary")
ax.set_ylim(0,1)
ax.set_xlabel("$\\lambda$ ($\\mu$m)")
ax.set_ylabel("Secondary Relative flux")


#%%
disc.params["dim"].value=32
pixmas = disc.getPixelSize(mas=True)

fov_model = 20
dim_model=int(fov_model/pixmas)
pix_model = pixmas

wl_im=np.array([0.5,1.5,4,10])*1e-6
t_im=np.array([0,0.25,0.5,0.75,1])
nwl=wl_im.size
nt=t_im.size
fig,ax=plt.subplots(nwl,nt,figsize=(2.2*nt,2.1*nwl))
disc.params["f"].value=1
fig,ax,im =morb.showModel(dim_model,pix_model,wl=wl_im,t=t_im,
                          fromFT=True,axe=ax,colorbar=False,normPow=1,
                          legend=True,normalize=True,cmap="hot",swapAxes=True)
disc.params["f"].value=1
t = np.linspace(T0, T0 + T, 100)
x, y = orb.getSeparation(t, mas=True)
for iwl in range(len(wl_im)):
    for it in range(len(t_im)):
        ax[it, iwl].plot(x, y, color="b", alpha=1, lw=1, ls="--")
        ax[it, iwl].get_xaxis().set_visible(False)
        ax[it, iwl].get_yaxis().set_visible(False)
fig.tight_layout()
plt.savefig(save_dir / "ExampleDISCO_orbit.png")



#%% Multi spatial frequencies and times for one (simultaneous) getComplexCoherentFlux call
nt = 500  # number of epochs in one period
t = np.linspace(T0, T0 + T, nt)

nB = 500  # number of baselines
wl = 2.1e-6  # fixed wavelength
B = np.linspace(0.0, 100, num=nB)
spf = B / wl

Bs = np.tile(B, (nt, 1))
ts = np.transpose(np.tile(t, (nB, 1)))
spf = Bs.flatten() / wl
spf0 = spf * 0

wls=wl+spf0

vis1 = np.abs(morb.getComplexCoherentFlux(spf,spf0, wl=wls,t=ts.flatten())).reshape(len(t), len(B))
vis1 /= np.outer(np.max(vis1, axis=1), np.ones(nB))
vis2 = np.abs(morb.getComplexCoherentFlux(spf0, spf,  wl=wls,t=ts.flatten())).reshape(len(t), len(B))
vis2 /= np.outer(np.max(vis2, axis=1), np.ones(nB))

#%% Plot the colormap of the temporal evolution of the visibility during one period
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

#%%

nwl = 100
nspf = 1000

# first zero of the primary used to define the spatial-frequency range
spf_res = 1.22 / (ud.params["d"].value * u.mas.to(u.rad)) 
spf = np.linspace(0, spf_res * 1.5, nspf)
wl = np.logspace(-6.5, -5, num=nwl)
spf_2D = np.tile(spf, (nwl, 1)).flatten()
wl_2D = np.transpose(np.tile(wl, (nspf, 1))).flatten()

ts=[0,0.25,0.5,0.75]
nt=len(ts)
viss=[]
for i in range(nt):
    vis = np.abs(morb.getComplexCoherentFlux(spf_2D, spf_2D*0, wl_2D,spf_2D*0+ts[i])).reshape(nwl, nspf)
    vis /= np.outer(np.max(vis, axis=1), np.ones(nspf))
    viss.append(vis)

fig, ax = plt.subplots(1, nt, figsize=(3*nt, 2.5))
for i in range(nt):
    sc = ax[i].scatter(spf_2D, viss[i], c=wl_2D*1e6, s=0.2, cmap="plasma")
    ax[i].set_xlabel("B/$\\lambda$ (cycles/rad)")
    ax[i].margins(0, 0)
    ax[i].set_title(f"Time = {ts[i]} T")
ax[0].set_ylabel("Visiblity")
plt.tight_layout()
fig.colorbar(sc, ax=ax, label="$\\lambda$ ($\\mu$m)")
plt.savefig(save_dir / "exampleDISCO_binary_visi.png")
