# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import oimodeler as oim
from scipy.interpolate import interp1d

import os
np.random.seed(1)

path = os.path.dirname(oim.__file__)

#%%
class oimParamLinearRangeWl(oim.oimParamInterpolator):

    def _init(self, param, wl0=2e-6, dwl=1e-9,values=[], kind="linear",**kwargs):
       
        self.kind=kind

        n = len(values)
        self.wl0 = (oim.oimParam(**oim._standardParameters["wl"]))
        self.wl0.name="wl0"
        self.wl0.description="Initial wl of the range"
        self.wl0.value=wl0
        self.wl0.free=False
        
        self.dwl = (oim.oimParam(**oim._standardParameters["wl"]))
        self.dwl.name="dwl"
        self.dwl.description="wl step in range"
        self.dwl.value=dwl
        self.dwl.free=False

        self.values = []
        
        for i in range(n):
            self.values.append(oim.oimParam(name=param.name, value=values[i],
                                        mini=param.min,maxi=param.max,
                                        description=param.description,
                                        unit=param.unit, free=param.free,
                                        error=param.error))
            
    def _interpFunction(self, wl, t):

        vals=np.array([vi.value for vi in self.values])
        nwl=vals.size   
        wl0=np.linspace(self.wl0.value,self.wl0.value+self.dwl.value*nwl,num=nwl)
        
        return interp1d (wl0,vals,kind=self.kind,fill_value="extrapolate")(wl)

    
    def _getParams(self):
        params=[]
        params.extend(self.values)
        params.append(self.wl0)
        params.append(self.dwl)
        return params

oim._interpolator["rangeWl"]=oimParamLinearRangeWl

#%%
nref=10

c = oim.oimUD(d=oim.oimInterp('rangeWl',wl0=2e-6,kind="cubic",
                               dwl=5e-8,values=np.random.rand(nref)*3+4))
m   = oim.oimModel(c)



print(m.getParameters())
print(m.getFreeParameters())


#%%
nB=200
B=np.linspace(0,60,num=nB)
nwl=1000
wl=np.linspace(2.0e-6,2.5e-6,num=nwl)
Bx_arr=np.tile(B[None,:], (nwl, 1)).flatten()
wl_arr=np.tile(wl[:,None], (1, nB)).flatten()
spfx_arr=Bx_arr/wl_arr
spfy_arr=spfx_arr*0

#%%


wl0=wl0=np.linspace(c.params['d'].wl0.value,
                    c.params['d'].wl0.value+c.params['d'].dwl.value*nref,num=nref)

vals=np.array([vi.value for vi in c.params['d'].values])


v=np.abs(m.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr).reshape(nwl,nB))

fig,ax=plt.subplots(2,1)
ax[0].plot(wl*1e6,c.params['d'](wl,0),color="r",label="interpolated param")
ax[0].scatter(wl0*1e6,vals,marker=".",color="k",label="reference values")

ax[0].set_ylabel("UD (mas)")
ax[0].get_xaxis().set_visible(False)
ax[0].legend()

for iB in range(1,nB):
    ax[1].plot(wl*1e6,v[:,iB]/v[:,0],color=plt.cm.plasma(iB/(nB-1)))
    
   
ax[1].set_xlabel("$\lambda$ ($\mu$m)")   
ax[1].set_ylabel("Visibility")


norm = colors.Normalize(vmin=np.min(B[1:]),vmax=np.max(B))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="Baseline Length (m)")


plt.savefig(os.path.join(path, os.pardir, "images","createInterp1.png"))

