# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:00:52 2024

@author: ame
"""

import oimodeler as oim
from pprint import pprint as print
ud = oim.oimUD(d=10, f=0.95)
pt = oim.oimPt(x=15, y=7, f=0.05)

mbin = oim.oimModel(pt,ud)

#%%
print(ud.params)
print(pt.params)
#%%
print(mbin.getParameters())
#%%

im = mbin.showModel(32,1.2,fromFT=False,normPow=0.5,normalize=True,cmap="hot")
