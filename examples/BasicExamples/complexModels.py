import oimodeler as oim
import numpy as np
import matplotlib.pyplot as plt
import os

path = os.path.dirname(oim.__file__)

fromFT=False
nB = 500  # number of baselines
nwl = 100  # number of walvengths

# Create some spatial frequencies
wl = np.linspace(3e-6, 4e-6, num=nwl)
B = np.linspace(1, 400, num=nB)
Bs = np.tile(B, (nwl, 1)).flatten()
wls = np.transpose(np.tile(wl, (nB, 1))).flatten()
spf = Bs/wls
spf0 = spf*0


# %%
g = oim.oimGauss(fwhm=oim.oimInterpWl([3e-6, 4e-6], [2, 8]))
print(g.params['fwhm']([3e-6, 3.5e-6, 4e-6, 4.5e-6]))

# %%
mg = oim.oimModel([g])

figGim, axGim,im = mg.showModel(256, 0.1, wl=[3e-6, 3.5e-6, 4e-6, 4.5e-6],swapAxes=True, figsize=(2.5, 2.5),fromFT=fromFT,
   savefig=os.path.join(path, os.pardir, "images", "complexModel_chromaticGaussian.png"))


# %%
vis = np.abs(mg.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))

figGv, axGv = plt.subplots(1, 1, figsize=(14, 8))
sc = axGv.scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
figGv.colorbar(sc, ax=axGv, label="$\\lambda$ ($\\mu$m)")
axGv.set_xlabel("B/$\\lambda$ (cycles/rad)")
axGv.set_ylabel("Visiblity")
axGv.margins(0, 0)

plt.savefig(os.path.join(path, os.pardir, "images",
            "complexModel_chromaticGaussianVis.png"))


# %%

ud = oim.oimUD(d=0.5, f=oim.oimInterpWl([3e-6, 4e-6], [2, 0.2]))

m2 = oim.oimModel([ud, g])
fig2im, ax2im,im = m2.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],swapAxes=True, normPow=0.2, figsize=(2.5, 2.5),fromFT=fromFT,
   savefig=os.path.join(path, os.pardir, "images", "complexModel_UDAndGauss.png"))

vis = np.abs(m2.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))

fig2v, ax2v = plt.subplots(1, 1, figsize=(14, 8))
sc = ax2v.scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
fig2v.colorbar(sc, ax=ax2v, label="$\\lambda$ ($\\mu$m)")
ax2v.set_xlabel("B/$\\lambda$ (cycles/rad)")
ax2v.set_ylabel("Visiblity")
ax2v.margins(0, 0)
ax2v.set_ylim(0, 1)
plt.savefig(os.path.join(path, os.pardir, "images",
            "complexModel_UDAndGaussVis.png"))

# %%
eg = oim.oimEGauss(fwhm=oim.oimInterpWl([3e-6, 4e-6], [2, 8]), elong=2, pa=90)
el = oim.oimEllipse(d=0.5, f=oim.oimInterpWl([3e-6, 4e-6], [2, 0.2]), elong=2, pa=90)

m3 = oim.oimModel([el, eg])
fig3im, ax3im,im = m3.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6], figsize=(2.5, 2.5), normPow=0.2,fromFT=fromFT,
   savefig=os.path.join(path, os.pardir, "images", "complexModel_Elong.png"),swapAxes=True)


# %%
fig3v, ax3v = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)

# East-West
vis = np.abs(m3.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
ax3v[0].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
ax3v[0].set_title("East-West Baselines")
ax3v[0].margins(0, 0)
ax3v[0].set_ylim(0, 1)
ax3v[0].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax3v[0].set_ylabel("Visiblity")

# North-South
vis = np.abs(m3.getComplexCoherentFlux(
    spf*0, spf, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
sc = ax3v[1].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
ax3v[1].set_title("North-South Baselines")
ax3v[1].set_xlabel("B/$\\lambda$ (cycles/rad)")
fig3v.colorbar(sc, ax=ax3v.ravel().tolist(), label="$\\lambda$ ($\\mu$m)")



    
plt.savefig(os.path.join(path, os.pardir,
            "images", "complexModel_ElongVis.png"))


#%%

print(m3.getFreeParameters())

#%%

eg.params['elong']=el.params['elong']
eg.params['pa']=el.params['pa']

print(m3.getFreeParameters())


#%%
    
er = oim.oimERing()

er.params['elong']=eg.params['elong']
er.params['pa']=oim.oimParamLinker(eg.params["pa"],"add",90)
er.params['din']=oim.oimParamLinker(el.params["d"],"mult",2)
er.params['dout']=oim.oimParamLinker(el.params["d"],"mult",4)

m4= oim.oimModel([el, eg,er])

fig4im, ax4im,im = m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6], figsize=(2.5, 2.5), normPow=0.2,swapAxes=True,fromFT=fromFT,
   savefig=os.path.join(path, os.pardir, "images", "complexModel_link.png"))    
    

print(m4.getFreeParameters())

#%%
el.params['d'].value=4
el.params['pa'].value=45
    
fig5im, ax5im,im = m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],  figsize=(2.5, 2.5),swapAxes=True, normPow=0.2,fromFT=fromFT,
       savefig=os.path.join(path, os.pardir, "images", "complexModel_linkRotScale.png"))    
     


#%%

  
 
  
