# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:37:28 2024

@author: ame
"""

import oimodeler as oim
from pathlib import Path

path = Path(__file__).parent
path_save= path.parent / "source"

#%%
res=oim.listDataFilters(details=True,save2csv=path_save /"table_dataFilter.csv")
res=oim.listComponents(details=True,save2csv=path_save /"table_components_fourier.csv",componentType="fourier")
res=oim.listComponents(details=True,save2csv=path_save /"table_components_image.csv",componentType="image")
res=oim.listComponents(details=True,save2csv=path_save /"table_components_radial.csv",componentType="radial")
res=oim.listFitters(details=True,save2csv=path_save /"table_fitters.csv")
res=oim.listParamInterpolators(details=True,save2csv=path_save /"table_interpolators.csv")

