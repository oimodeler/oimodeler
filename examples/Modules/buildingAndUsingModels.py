# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:00:52 2024

@author: ame
"""

import oimodeler as oim


ud1 = oim.oimUD(x=15, y=7, d=5, f=0.5)
ud2 = oim.oimUD(d=10, f=0.5)

model = oim.oimModel(ud1,ud2)
#%%
model.components
print(model)
print(model.shortname)
#%%

im = model.showModel(512,0.1,fromFT=True,normPow=1,normalize=True)
