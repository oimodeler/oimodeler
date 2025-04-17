from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
import time
from astropy.io import fits

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


#%% getting the list of all Fourier-based components currently available.
print(oim.listComponents(componentType="fourier"))

#%% Getting information on a component

help(oim.oimUD)

#%% getting the list of all image-based components currently available.

print(oim.listComponents(componentType="image"))

#%% creating and showing images from a few components

spiral = oim.oimSpiral(dim=128, fwhm=20, P=0.01, width=0.1, pa=30, elong=2)
mspiral = oim.oimModel(spiral)

frot = oim.oimFastRotator(dpole=5, dim=128, incl=-50,rot=0.99, Tpole=20000, beta=0.25,pa=20)
mfrot = oim.oimModel(frot)


radmc3D_fname = product_dir / "radmc3D_model.fits"
radmc3D = oim.oimComponentFitsImage(radmc3D_fname,pa=180)
mradmc3D = oim.oimModel(radmc3D)

figimages,aximages = plt.subplots(1,3,figsize=(15,4))

mspiral.showModel(256,0.3,axe=aximages[0],colorbar=False,normPow=1)
mfrot.showModel(256,0.05,wl=1e-6,axe=aximages[1],colorbar=False,normPow=1)
mradmc3D.showModel(256,0.5,wl=3e-6,axe=aximages[2],colorbar=False,normPow=0.22)

plt.tight_layout()
figimages.savefig(save_dir / "componentImages_images.png")


#%% check active and available FT backends
print(oim.oimOptions.ft.backend.available)
print(oim.oimOptions.ft.backend.active)

#%% changing FT backend
oim.oimOptions.ft.backend.active = oim.FFTWBackend # using oimOptions
oim.setFTBackend("fftw") # or with the setFTBackend specifying  the alias 

#%% Zero padding
oim.oimOptions.ft.padding = 8

#%% Padding effect on accuracxy of FFT

#creating the spiral model
spiral = oim.oimSpiral(dim=128, fwhm=20, P=1, width=0.1, pa=0, elong=1)
mspiral = oim.oimModel(spiral)

#creating a set of baselines from 0 to 100m in K band
wl = 2.1e-6
B = np.linspace(0, 100, num=1000)
spf = B/wl

#computing the reference model with padding of 32
oim.oimOptions.ft.padding = 32
ccf01 = mspiral.getComplexCoherentFlux(spf, spf*0)
v01 = np.abs(ccf01/ccf01[0])

start = time.time()
ccf02 = mspiral.getComplexCoherentFlux(spf*0, spf)
v02 = np.abs(ccf02/ccf02[0])
end = time.time()
dt0 = end -start

#%%  computing FFT with different padding

padding=[16,8,4,2,1]
figpad,axpad = plt.subplots(2,2, figsize=(10,5),sharey="row",sharex=True)

axpad[0,0].plot(spf, v01, color="k",lw=4)
axpad[0,1].plot(spf, v02, color="k",lw=4,label=f"padding=32x ({dt0*1000:.0f}ms)")

for pi in padding :
    oim.oimOptions.ft.padding = pi
    ccf1 = mspiral.getComplexCoherentFlux(spf, spf*0)
    v1 = np.abs(ccf1/ccf1[0])
    start = time.time()
    ccf2 = mspiral.getComplexCoherentFlux(spf*0, spf)
    v2 = np.abs(ccf2/ccf2[0])
    end = time.time()
    dt = end -start
    axpad[0,0].plot(spf, v1)
    axpad[0,1].plot(spf, v2,label=f"padding={pi}x ({dt*1000:.0f}ms)")
    axpad[1,0].plot(spf, (v1-v01)/v01*100)
    axpad[1,1].plot(spf, (v2-v02)/v02*100)
    

for i in range(2):
    axpad[1,i].set_xlabel("spatial frequency (cycles/rad)")
    axpad[1,i].set_yscale("symlog")
axpad[0,0].set_title("East-West baselines")
axpad[0,1].set_title("North-South baselines")
axpad[0,0].set_ylabel("Visbility")
axpad[0,1].legend()
axpad[1,0].set_ylabel("Residual (%)")



figpad.savefig(save_dir / "componentImages_padding.png")
oim.oimOptions.ft.padding = 4
