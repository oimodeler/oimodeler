from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim

oimodeler_dir = Path(oim.__file__).parent

# Creating a model including the spiral component
pl = oim.oimAsymmetricRadialPowerLaw(dim=128, din=2, dout=5,
                                     pa=30, elong=2, a=0.5, phi=45, f=0.8)
m = oim.oimModel(pl)

# Plotting images of the model
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
m.showModel(128, 0.1, swapAxes=True, fromFT=False,
            normPow=1, axe=ax[0], colorbar=False)
m.showModel(128, 0.1, swapAxes=True, fromFT=True,
            normPow=1, axe=ax[1], colorbar=False)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")
plt.show()
