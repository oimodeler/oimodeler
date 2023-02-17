import os

import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy import units as units


if __name__ == "__main__":
    #%% creating a model including the spiral component
    pl = oim.oimAsymmetricRadialPowerLaw(dim=128, din=2, dout=5,
                                         pa=30, elong=2, a=0.5, phi=45, f=0.8)
    m = oim.oimModel(pl)

    #%% plotting images of the model
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    m.showModel(128, 0.1, swapAxes=True, fromFT=False, normPow=1, axe=ax[0], colorbar=False)
    m.showModel(128, 0.1, swapAxes=True, fromFT=True, normPow=1, axe=ax[1], colorbar=False)
    ax[1].get_yaxis().set_visible(False)
    ax[0].set_title("Direct Image")
    ax[1].set_title("From FFT")
    plt.show()
    # fig.savefig(os.path.join(path, os.pardir, "images", "customCompImageSpiral.png"))

    #%% Computing and plotting visibilities for various baselines and walvelengths
    # nB, nwl, wl=5000
    # nwl=1
    # wl=0.5e-6

    # B=np.linspace(0,100,num=nB//2)
    # Bx=np.append(B,B*0)
    # By=np.append(B*0,B)

    # spfx=Bx/wl
    # spfy=By/wl

    # vc=m.getComplexCoherentFlux(spfx,spfy)
    # v=np.abs(vc/vc[0])

    # fig,ax=plt.subplots(1,1)
    # label=["East-West Baselines",]

    # ax.plot(B/wl/units.rad.to(units.mas),v[:nB//2],color="r",label="East-West Baselines")
    # ax.plot(B/wl/units.rad.to(units.mas),v[nB//2:],color="b",label="North-South Baselines")  

    # ax.set_xlabel("B/$\lambda$ (cycles/mas)")
    # ax.set_ylabel("Visibility")    
    # ax.legend()

    # fig.savefig(os.path.join(path,os.pardir,"images","customCompImageSpiralVis.png"))
