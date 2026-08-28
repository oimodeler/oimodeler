"""
Tests for the oimodeler.oimCustomComponents sub-package.
"""

import copy
import json
import warnings
from pathlib import Path

import pytest

from oimodeler.oimData import oimData
from oimodeler.oimDataFilter import (
    oimDataFilter,
    oimKeepDataTypeFilter,
    oimWavelengthBinningFilter,
    oimWavelengthRangeFilter,
)
from oimodeler.oimModel import oimModel
from oimodeler.oimSimulator import oimSimulator


class TestOimTempGrad:

    @pytest.fixture
    def data(self, global_data_dir: Path) -> oimData:
        """Data suited for an `oimCustomComponents.oimTempGrad` fit."""
        data = oimData(
            sorted((global_data_dir / "AS209_MATISSE").glob("*.fits"))
        )
        f1 = oimWavelengthRangeFilter(targets=[0, 2], wlRange=[3.2e-6, 3.8e-6])
        filt_bin_L = oimWavelengthBinningFilter(
            targets=[0, 2], bin=5, normalizeError=False
        )
        filt_bin_N = oimWavelengthBinningFilter(
            targets=1, bin=7, normalizeError=False
        )
        f2 = oimKeepDataTypeFilter(dataType=["FLUXDATA", "VISAMP"])
        data.setFilter(oimDataFilter([f1, f2, filt_bin_L, filt_bin_N]))
        return data

    @pytest.fixture
    def model(self, global_data_dir: Path) -> oimModel:
        """A model with a central star plus a temperature gradient disc."""
        with (
            warnings.catch_warnings(record=True),
            open(global_data_dir / "TempGrad_model.json", "r") as f,
        ):
            return oimModel.deserialize(json.load(f))

    @pytest.mark.parametrize("radius_cut", (1.5, 3, 5, 10))
    def test_single_vs_multi_component(
        self, data: oimData, model: oimModel, radius_cut: float
    ) -> None:
        """Tests the output of one `oimCustomComponents.oimTempGrad` to be similar to
        that of two that are adjoined to each other in radii."""
        complex_vis = model.getComplexCoherentFlux(
            data.vect_u, data.vect_v, data.vect_wl
        )
        sim = oimSimulator(data=data, model=model)
        sim.compute(computeChi2=True, computeSimulatedData=True)

        two_model = copy.deepcopy(model)
        tg2 = copy.deepcopy(two_model.components[-1])
        two_model.components[-1].rout.value = radius_cut
        tg2.rin.value = radius_cut
        tg2.elong.free = tg2.pa.free = False
        tg2.q.free = tg2.T0.free = False
        tg2.Mdust.free = tg2.p.free = False
        two_model.components.append(tg2)
        two_sim = oimSimulator(data=data, model=two_model)
        two_sim.compute(computeChi2=True, computeSimulatedData=True)
        breakpoint()
