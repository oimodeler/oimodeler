# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:04:48 2023

@author: Ame
"""

import numpy as np
import oimodeler as oim
import matplotlib.pyplot as plt
import astropy.units as u

import matplotlib.colors as colors
import matplotlib.cm as cm

# setting the FFT parameters and backend. As the kinematicsDisk model as no 
# outer sharp edge we decide not to zero-pad the images during the FFT process
# no oder to save time and without noticeable effect on the simulated visibilities
oim.oimOptions['FTpaddingFactor']=2
oim.oimOptions['FTBackend']=oim.FFTWBackend

# Global simulation parameters
dim   = 256	            # size of the simulation in pixel
fov   = 20              # field of view in D*
nwl   = 101             # number of wavelengths
wl0   = 2.1656e-6       # line central wavelength in m
dwl   = 0.5e-10         # delta lam in m 
res   = 1.8e-10         # spectral resolution in m 

# Central Star
incl  = 40              # inclination angle 
pa    = 0               # disk minor-axis postion angle in degrees 
Rstar = 6               # stellat radius in solar radii 
dist  = 150             # distance in parsecs

# Disk geometry in the continuum 
fwhmCont     = 3            # major-axis FWHM or diameter in D* 
fluxDiskCont = 0.5          # disk flux ratio in the continuum 

# Disk geometry in the line
fwhmLine      = 8.          # circumstellar disk FWHM in the line in D*
EW            = 20.         # line equivalent width in Angstroms


# disk Kinematics
vrot  =  300               # stellar rotational velocity (in km/s  )
beta  = -0.5               # exponent of the rotational velocity law 
v0    =  0                 # radial velocity at the stellar surface (in km/s)
vinf  =  0                 # terminal radial velocity (in km/s)
gamma =  0.86              # exponent of the CAK type of radial velocity law

disk = oim.oimKinematicDisk(dim=dim, fov=fov, nwl=nwl, wl0=wl0, dwl=dwl, res=res,
                       incl=incl, pa=pa, Rstar=Rstar, dist=dist,
                       fwhmCont=fwhmCont, fluxDiskCont=fluxDiskCont, 
                       EW=EW, fwhmLine=fwhmLine, vrot=vrot, beta=beta, 
                       v0=v0, vinf=vinf, gamma=gamma)

m = oim.oimModel(disk)
#%% get the internal image-cube, wavelength table and pixSize for this model

im0 = disk._internalImage()
wl  = disk._wl                          # wl in meters
pixSize = disk._pixSize*u.rad.to(u.mas) # pixel size in rad converted to mas

print(f"(ntime, nwl, dim, dim) = {im0.shape}")
print(f"pixSize = {pixSize:.3f} mas")
print(f"wl = [ {wl[0]*1e9} ... {wl[-1]*1e9}] nm")

#%% plot the line profile
spectra = np.sum(im0,axis=(0, 2, 3))

fig, ax = plt.subplots(1,1)
ax.plot(wl*1e9, spectra)
ax.set_xlabel("$\\lambda ($nm)")
ax.set_ylabel("Normalized flux")
ax.set_title( "Line profile")

#%% plot some images through the line
c = 3e5 # light speed in km/s

dopplerShift = (wl - wl0)/wl0 * c

# indexzes of the the walvengths for which we will plot images
index  = [6, 21, 35, 50, 65, 79, 94]
nindex = len(index)

fig, ax = plt.subplots(1,nindex,figsize=(15,15/nindex))
for i in range(nindex):
    ax[i].imshow(im0[0,index[i],:,:]**0.5)   
    ax[i].text(dim/2,0,f"{wl[index[i]]*1e9:.2f} nm",ha="center",va="top",color="w")
    ax[i].text(dim/2,dim,f"{dopplerShift[index[i]]:.2f} km/s",ha="center",va="bottom",color="w")
    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_visible(False)



#%% same plots but using a model class that allow interpolation between wavelengths and arbitary pixel size
n = 7
dopplerShift2 = np.linspace( -300, 300, n)
wl2 = wl0 + dopplerShift2/c*wl0 


ims = m.getImage(dim, pixSize, wl2)

fig, ax = plt.subplots(1,nindex,figsize=(15,15/nindex))
for i in range(n):
    ax[i].imshow(ims[i,:,:]**0.5)   
    ax[i].text(dim/2,0,f"{wl2[i]*1e9:.2f} nm",ha="center",va="top",color="w")
    ax[i].text(dim/2,dim,f"{dopplerShift2[i]:.2f} km/s",ha="center",va="bottom",color="w")
    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_visible(False)

#%% Creating spatial frequencies for two sets of perpendicular baselines 
# ranging from 0 to 100m and with the previouslt defined wavelength-range

nB = 20
wls = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
B = np.linspace(0, 200, num=nB//2)

Bx = np.append(B, B*0)  # 1st half of B array are baseline in the East-West orientation
By = np.append(B*0, B)  # 2nd half are baseline in the North-South orientation

Bx_arr = np.tile(Bx[None, :], (nwl, 1)).flatten()
By_arr = np.tile(By[None, :], (nwl,  1)).flatten()
wl_arr = np.tile(wl[:, None], (1, nB)).flatten()

spfx = Bx_arr/wl_arr
spfy = By_arr/wl_arr

#%%
vc = m.getComplexCoherentFlux(spfx, spfy, wl_arr)

v = np.abs(vc.reshape(nwl, nB))
v = v/np.tile(v[:, 0][:, None], (1, nB))

phi = np.angle(vc.reshape(nwl, nB))
phi = np.rad2deg(phi - np.tile(phi[0, :][None, :], (nwl,1)))
#%%

fig, ax = plt.subplots(1, 2, figsize=(12, 4),sharey=True)
titles = ["East-West Baselines", "North-South Baselines"]

for iB in range(nB):
    cB = (iB % (nB//2))/(nB//2-1)
    ax[2*iB//nB].plot(wl*1e9, v[:, iB],color=plt.cm.plasma(cB))

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel(r"$\lambda$ (nm)")
ax[0].set_ylabel("Visibility")
ax[1].get_yaxis().set_visible(False)

norm = colors.Normalize(vmin=np.min(B), vmax=np.max(B))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="B (m)")


#%%

fig, ax = plt.subplots(1, 2, figsize=(12, 4),sharey=True)
titles = ["East-West Baselines", "North-South Baselines"]

for iB in range(nB):
    cB = (iB % (nB//2))/(nB//2-1)
    ax[2*iB//nB].plot(wl*1e9, phi[:, iB],color=plt.cm.plasma(cB))

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel(r"$\lambda$ (nm)")
ax[0].set_ylabel("Differential Phase (deg)")

ax[1].get_yaxis().set_visible(False)

norm = colors.Normalize(vmin=np.min(B), vmax=np.max(B))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="B (m)")



#%%

fig, ax = plt.subplots()

maxphi=np.max(phi[:,:nB//2],axis=0)
maxi=np.max(maxphi)

plt.plot(B[:nB//2],maxphi,label="simulated data")
#plt.plot(pa,maxi*np.abs(np.sin(np.deg2rad(pa))),label="abs of sine function",ls=":")
ax.set_ylabel("Maximum of phase variation (deg)")
ax.set_xlabel("Baseline length (m)")
ax.legend()





#%% Creating spatial frequencies for a set of 50m with various orientation 
B0=50

nB = 200
wls = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
pa = np.linspace(-180, 180, num=nB)

Bx = B0*np.sin(np.deg2rad(pa)) 
By = B0*np.cos(np.deg2rad(pa)) 

Bx_arr = np.tile(Bx[None, :], (nwl, 1)).flatten()
By_arr = np.tile(By[None, :], (nwl,  1)).flatten()
wl_arr = np.tile(wl[:, None], (1, nB)).flatten()

spfx = Bx_arr/wl_arr
spfy = By_arr/wl_arr

#%%
vc = m.getComplexCoherentFlux(spfx, spfy, wl_arr)

v = np.abs(vc.reshape(nwl, nB))
v = v/np.tile(spectra[:, None], (1, nB))

phi = np.angle(vc.reshape(nwl, nB))
phi = np.rad2deg(phi - np.tile(phi[0, :][None,:], (nwl, 1)))

#%%

fig, ax = plt.subplots(1, 2, figsize=(12, 4),sharex=True)

for iB in range(nB):
    cB = (pa[iB])/90.
    ax[0].plot(wl*1e9, v[:, iB],color=plt.cm.plasma(cB))
    
for iB in range(nB):
    cB = (pa[iB])/90.
    ax[1].plot(wl*1e9, phi[:, iB],color=plt.cm.plasma(cB))
    
for i in range(2):
    ax[i].set_xlabel(r"$\lambda$ (nm)")
ax[0].set_ylabel("Visibility")
ax[1].set_ylabel("Differential Phase (deg)")
fig.suptitle(f"Visiblity and phase and phase \nfor a {B0}m baseline depending on its orientation")
norm = colors.Normalize(vmin=np.min(pa), vmax=np.max(pa))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="P.A. (deg)")

#%%

fig, ax = plt.subplots()

maxphi=np.max(phi,axis=0)
maxi=np.max(maxphi)

plt.plot(pa,maxphi,label="simulated data")
plt.plot(pa,maxi*np.abs(np.sin(np.deg2rad(pa))),label="abs of sine function",ls=":")
ax.set_ylabel("Maximum of phase variation (deg)")
ax.set_xlabel("disk major-axis P.A. (deg)")
ax.legend()

