"""
Tests for the oimodeler.oimCustomComponents sub-package.
"""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import NDArray

from oimodeler.oimBasicFourierComponents import oimPt
from oimodeler.oimCustomComponents import oimTempGrad
from oimodeler.oimData import oimData
from oimodeler.oimDataFilter import (
    oimDataFilter,
    oimKeepDataTypeFilter,
    oimWavelengthBinningFilter,
    oimWavelengthRangeFilter,
)
from oimodeler.oimModel import oimModel
from oimodeler.oimParam import oimInterp
from oimodeler.oimSimulator import oimSimulator


class TestOimTempGrad:
    """Tests `oimCustomComponents.oimTempGrad`."""

    # TODO: Should this data be saved in the `test_data_dir`
    @pytest.fixture(scope="module")
    def data(self, global_data_dir: Path) -> oimData:
        """Data suited for temperature gradient."""
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

    @pytest.fixture(scope="module")
    def vis_base(self, test_data_dir: Path) -> NDArray[np.float64]:
        """Baseline correlated fluxes for temperature gradient
        to detect if model executes correctly. Should be equal to test."""
        return np.load(test_data_dir / "TempGradVis.npy")

    @pytest.fixture(scope="module")
    def star(self) -> oimPt:
        """A central star."""
        return oimPt(f=oimInterp("starWl", T=4600, L=2.23, dist=121))

    @pytest.fixture(scope="module")
    def temp_grad_kwargs(
        self, global_data_dir: Path
    ) -> dict[str, float | oimInterp]:
        """Parameters for the base temperature gradient."""
        opac_file = (
            global_data_dir
            / "FSCMa_MATISSE"
            / "dustkappa_olivine_graphite_1_20.inp"
        )
        op_wl, op = np.loadtxt(opac_file, usecols=[0, 1], unpack=True)
        return {
            "dim": 128,
            "dist": 121,
            "T0": 350,
            "rin": 0.1,
            "rout": 30,
            "p": -0.5841548005057318,
            "q": -0.4493209018837913,
            "pa": 72.52114779136605,
            "elong": 2.5123019579326464,
            "Mdust": 3.412703968363078e-08,
            "kappa_abs": oimInterp("wl", wl=op_wl * 1e-6, values=op),
        }

    @pytest.fixture(scope="module")
    def temp_grad_model(
        self, star: oimPt, temp_grad_kwargs: dict[str, float | oimInterp]
    ) -> oimModel:
        """Base temperature gradient model."""
        comp_kwargs = copy.deepcopy(temp_grad_kwargs)
        tg = oimTempGrad(**comp_kwargs)
        tg.rin.free = tg.rout.free = tg.T0.free = False
        return oimModel([star, tg])

    def add_disc(self, model: oimModel, rcut: float) -> oimModel:
        """Helper function to add multiple components to continuous disc
        but retain the same overall structure."""
        tg = copy.deepcopy(model.components[-1])
        for param in tg.params.values():
            param.free = False

        model.components[-1].rout.value = tg.rin.value = rcut
        model.components.append(tg)
        return model

    # TODO: Save the complex vis as well to compare if something changes in the computation over tim
    @pytest.mark.parametrize(
        "params",
        (
            {"cosi": 0.3980413249460237, "flat": True},
            {
                "log_sigma0": np.log10(0.0005537958068277064),
                "compute_sigma0": False,
            },
            {
                "cosi": 0.3980413249460237,
                "flat": True,
                "log_sigma0": np.log10(0.0005537958068277064),
                "compute_sigma0": False,
            },
        ),
    )
    def test_init(
        self,
        data: oimData,
        vis_base: NDArray[np.float64],
        star: oimPt,
        temp_grad_model: oimModel,
        temp_grad_kwargs: dict[str, float | oimInterp],
        global_data_dir: Path,
        params: dict[str, float],
    ) -> None:
        """Tests if models initialised with `compute_sigma0` and `flat` are
        identical to base."""
        comp_kwargs = copy.deepcopy(temp_grad_kwargs)
        sim_temp_grad = oimSimulator(data=data, model=temp_grad_model)
        sim_temp_grad.compute(computeChi2=True, computeSimulatedData=True)
        vis_temp_grad = np.abs(
            temp_grad_model.getComplexCoherentFlux(
                data.vect_u, data.vect_v, data.vect_wl
            )
        )
        assert np.allclose(vis_temp_grad, vis_base, atol=1e-1)

        if "cosi" in params:
            del comp_kwargs["elong"]

        if "log_sigma0" in params:
            del comp_kwargs["Mdust"]

        tg = oimTempGrad(**comp_kwargs, **params)
        tg.rin.free = tg.rout.free = tg.T0.free = False
        model = oimModel([star, tg])
        sim = oimSimulator(data=data, model=model)
        sim.compute(computeChi2=True, computeSimulatedData=True)
        vis = np.abs(
            model.getComplexCoherentFlux(
                data.vect_u, data.vect_v, data.vect_wl
            )
        )

        assert sim_temp_grad.chi2r == pytest.approx(sim.chi2r)
        assert np.allclose(vis_temp_grad, vis)

    # NOTE: For more than one disc, dim of >= 2**10 is required
    @pytest.mark.slow
    @pytest.mark.parametrize("dim", (1024, 2048))
    @pytest.mark.parametrize(
        "cut_radii", ([1.5], [3], [5], [10], [20], [1.5, 5], [1.5, 5, 15])
    )
    def test_single_vs_multi_component(
        self,
        data: oimData,
        star: oimPt,
        temp_grad_kwargs: dict[str, float | oimInterp],
        cut_radii: list[float],
        dim: int,
    ) -> None:
        """Tests the output of a continuous disc to be similar/identical
        to one that has multiple components that simulate a continuous disc."""
        comp_kwargs = copy.deepcopy(temp_grad_kwargs)
        comp_kwargs["dim"] = dim
        del comp_kwargs["Mdust"]
        del comp_kwargs["elong"]

        tg = oimTempGrad(
            cosi=0.3980413249460237,
            log_sigma0=np.log10(0.0005537958068277064),
            compute_sigma0=False,
            flat=True,
            **comp_kwargs,
        )
        tg.rin.free = tg.rout.free = tg.T0.free = False
        model = oimModel([star, tg])

        sim = oimSimulator(data=data, model=model)
        sim.compute(computeChi2=True, computeSimulatedData=True)

        vis = np.abs(
            model.getComplexCoherentFlux(
                data.vect_u, data.vect_v, data.vect_wl
            )
        )

        model2 = copy.deepcopy(model)
        for rcut in cut_radii:
            model2 = self.add_disc(model2, rcut)

        sim2 = oimSimulator(data=data, model=model2)
        sim2.compute(computeChi2=True, computeSimulatedData=True)

        vis2 = np.abs(
            model2.getComplexCoherentFlux(
                data.vect_u, data.vect_v, data.vect_wl
            )
        )
        assert np.isclose(sim.chi2r, sim2.chi2r, atol=1e-2)
        assert np.allclose(vis, vis2, atol=1e-2)

        # NOTE: Tests congurency of the logarithmic grids (i.e. same ratio + endpoints)
        radii = model.components[1]._r / 1e3 * model.components[1].dist.value
        radii2 = [
            model2.components[1]._r / 1e3 * model2.components[1].dist.value
        ]
        for i in range(2, len(cut_radii) + 2):
            tmp_radii = (
                model2.components[i]._r / 1e3 * model2.components[i].dist.value
            )
            assert radii2[-1][-1] == pytest.approx(tmp_radii[0])
            radii2.append(tmp_radii)

        assert radii[0] == pytest.approx(radii2[0][0])
        assert radii[-1] == pytest.approx(radii2[-1][-1])
        assert np.allclose(
            np.diff(np.log(radii)), np.diff(np.log(radii2)).sum(0)
        )
