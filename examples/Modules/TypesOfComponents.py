from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


#%% getting the list of all Fourier-based component currently available.
print(oim.listComponents(componentType="fourier"))