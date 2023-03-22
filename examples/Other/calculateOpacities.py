from pathlib import Path

import oimodeler as oim

# Get wavelength solution from file
path = Path(oim.__file__).parent.parent
# files = list(map(str, (path / "examples" / "testData" / "FSCMa_MATISSE").glob("*.fits")))
file = path / "examples" / "testData" / "Optool" / "dustkappa.dat"
# data = oim.oimData(files[0])
# data.prepareData()

test = oim.oimOptoolBackend(grains="pyr")
breakpoint()
