# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim


path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "FSCMa_MATISSE"

files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)
pprint(data.data)

data.prepareData()
pprint(data.vect_u)
pprint(data.vect_v)
pprint(data.vect_wl)
pprint(data.vect_u.shape)
