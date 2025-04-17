# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 10:22:47 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
import time
from astropy.io import fits

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


#%% getting the list of all radial-profile-based components currently available.
print(oim.listComponents(componentType="radial"))