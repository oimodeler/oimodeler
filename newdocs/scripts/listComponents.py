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
oim.listDataFilters(details=True,save2csv=path_save /"table_dataFilter.csv")
oim.listComponents(details=True,save2csv=path_save /"table_components_fourier.csv",componentType="fourier")
oim.listComponents(details=True,save2csv=path_save /"table_components_image.csv",componentType="image")
oim.listComponents(details=True,save2csv=path_save /"table_components_radial.csv",componentType="radial")
oim.listFitters(details=True,save2csv=path_save /"table_fitters.csv")
oim.listParamInterpolators(details=True,save2csv=path_save /"table_interpolators.csv")

