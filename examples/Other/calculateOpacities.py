import oimodeler as oim
test = oim.oimOptoolBackend("pyr", computational_method="mmf",
                            monomer_radius=0.01, gs_min=0.05, gs_max=3000)
print(test.kabs)
