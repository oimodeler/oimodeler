from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


#%% getting the list of all Fourier-based component currently available.
print(oim.listComponents(componentType="fourier"))

#%%

print(oim.listComponents(componentType="image"))

#%%

spiral = oim.oimSpiral(dim=256, fwhm=20, P=0.1, width=0.2, pa=30, elong=2)
mspiral = oim.oimModel(spiral)

frot = oim.oimFastRotator(dpole=5, dim=128, incl=-50,rot=0.99, Tpole=20000, beta=0.25,pa=20)
mfrot = oim.oimModel(frot)


radmc3D_fname = product_dir / "radmc3D_model.fits"
radmc3D = oim.oimComponentFitsImage(radmc3D_fname,pa=180)
mradmc3D = oim.oimModel(radmc3D)

figimages,aximages = plt.subplots(1,3,figsize=(15,4))

mspiral.showModel(256,0.2,normPow=1,axe=aximages[0],colorbar=False)
mfrot.showModel(256,0.05,normPow=1,axe=aximages[1],colorbar=False)
mradmc3D.showModel(256,0.5,wl=3e-6,normPow=0.25,axe=aximages[2],colorbar=False)

plt.tight_layout()
figimages.savefig(save_dir / "componentImages_images.png")


#%%

