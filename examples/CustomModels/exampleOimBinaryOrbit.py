# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 14:42:21 2026

@author: ame
"""
import oimodeler as oim
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

t0=0 # should be MJD if data model comparison (arbitrary if not)
e = 0.5
a = 10
T = 1. # days by default but can be changed modifying  orb.param["T"].unit 
i = 45.
o = 0.
O = 0.

nt2=13
t = np.linspace(t0,t0+T,100)
t2 = np.linspace(t0,t0+T,nt2)

orb = oim.oimBinaryOrbit(e=e, a=a, T=T, i=i, o=o,O=O)

x,y = orb._getSeparation(t,mas=True)
x2,y2 = orb._getSeparation(t2,mas=True)        
#%% plot the full orbit
fig,ax = plt.subplots(figsize=(7,7))
ax.plot(x,y,color='k',alpha=0.5)
ax.scatter(x2[:-1],y2[:-1],color='r',marker=".",zorder=10)

for i in range(nt2-1):
    ax.annotate(f"  t={t2[i]:.2f}",(x2[i],y2[i]),color='r',zorder=10,fontsize=8)
ax.axis('equal')


ax.scatter(0,0,marker='o',color='k')
ax.set_xlabel('x (mas)')
ax.set_ylabel('y (mas)')
ax.set_xlim([10,-10])

#%% Baselines at different time
N=200
wl = 2.1e-6
B = np.linspace(0.0, 100, num=N)
spf = B/wl

orb = oim.oimBinaryOrbit(e=e, a=a, T=T, i=i, o=o,O=O)
orb.primary = oim.oimUD(d=3)
orb.secondary.params["f"].value=0.3

morb = oim.oimModel(orb)
fig,ax = plt.subplots(3,nt2,figsize=(nt2*2.5,8))


mascycle=1./u.rad.to(u.mas)

for iT in range(nt2):
    for iPA in range(2):
        v = np.abs(morb.getComplexCoherentFlux(spf*(iPA==1),spf*(iPA==0),t=np.ones(N)*t2[iT]))
        v = v/np.max(v)
        ax[iPA,iT].plot(spf*mascycle,v,alpha=0.5)
        ax[iPA,iT].set_ylim(0,1)

    ax[0,iT].text(0.5,1.01,f"x={x2[iT]:.0f}mas  y={y2[iT]:.0f}mas", transform=ax[0,iT].transAxes,va="bottom",ha="center")
    ax[0,iT].text(0.5,1.10,f"PA={iPA*90}°  t={t2[iT]:.2f}day", transform=ax[0,iT].transAxes,va="bottom",ha="center")
    ax[1,iT].set_xlabel("B/$\\lambda$ (cycle/mas)")
    if iT!=0:
        ax[0,iT].axes.get_yaxis().set_visible(False)
        ax[1,iT].axes.get_yaxis().set_visible(False)        

ax[0,0].set_ylabel("East-West Vis.")
ax[1,0].set_ylabel("North-South Vis.")
dim=64
fov=30
pix=fov/dim

for iT in range(nt2):
    morb.showModel(dim,pix,t=t2[iT],fromFT=True,axe=ax[2,iT],normPow=0.5,colorbar=False)
    if iT!=0:
        ax[2,iT].axes.get_yaxis().set_visible(False)

plt.tight_layout()

#%% multi spafreq and time simulatenously in one getComplexCoherentFlux call

nt=500
t = np.linspace(t0,t0+T,nt)
nB = 500  # number of baselines
wl = 2.1e-6
B = np.linspace(0.0, 100, num=nB)
spf = B/wl

Bs = np.tile(B, (nt, 1))
ts = np.transpose(np.tile(t, (nB, 1)))
spf = Bs.flatten()/wl
spf0 = spf*0

vis1 = np.abs(morb.getComplexCoherentFlux(spf, spf*0, t=ts.flatten())).reshape(len(t), len(B))
vis1 /= np.outer(np.max(vis1, axis=1), np.ones(nB))
vis2 = np.abs(morb.getComplexCoherentFlux(spf*0, spf, t=ts.flatten())).reshape(len(t), len(B))
vis2 /= np.outer(np.max(vis2, axis=1), np.ones(nB))
#%%
fig, ax = plt.subplots(1,2,figsize=(15,8),sharex=True,sharey=True)

ax[0].pcolormesh(Bs,ts,vis1, cmap="plasma")
sc =ax[1].pcolormesh(Bs,ts,vis2, cmap="plasma")

fig.colorbar(sc, ax=ax, label="Visiblity")

ax[0].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[0].set_ylabel("time (days)")
ax[0].set_title("East-West Baselines")
ax[1].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[1].set_ylabel("time (days)")
ax[1].set_title("North-South Baselines")