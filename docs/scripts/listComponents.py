# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:37:28 2024

@author: ame
"""

from pathlib import Path

import oimodeler as oim

path = Path(__file__).parent.parent / "source"

# %%
res = oim.listDataFilters(details=True, save2csv=path / "table_dataFilter.csv")
res = oim.listComponents(
    details=True,
    save2csv=path / "table_components_fourier.csv",
    componentType="fourier",
)
res = oim.listComponents(
    details=True,
    save2csv=path / "table_components_image.csv",
    componentType="image",
)
res = oim.listComponents(
    details=True,
    save2csv=path / "table_components_radial.csv",
    componentType="radial",
)
res = oim.listFitters(details=True, save2csv=path / "table_fitters.csv")
res = oim.listParamInterpolators(
    details=True, save2csv=path / "table_interpolators.csv"
)
