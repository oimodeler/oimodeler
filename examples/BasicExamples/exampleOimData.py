# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim

path = Path(oim.__file__).parent.parent
pathData = path / Path().parent / "examples" / "testData" / "FSCMa_MATISSE"

# TODO: After pathlib change of all `oimodeler` modules, remove str here
files = list(map(str, pathData.glob("*.fits")))

# TODO: After pathlib change of all `oimodeler` modules, remove str here
data = oim.oimData(files)
pprint(data.data)

data.prepareData()
pprint(data.vect_u)
pprint(data.vect_v)
pprint(data.vect_wl)
pprint(data.vect_u.shape)
