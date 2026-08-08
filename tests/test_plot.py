"""
Runs basic plots (from the oimodeler.oimModel, oimodeler.oimPlot, and
oimodeler.oimSimulator modules) to test if anything is broken.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from oimodeler.oimBasicFourierComponents import (
    oimGauss,
    oimIRing,
    oimPt,
    oimUD,
)
from oimodeler.oimComponent import oimComponent
from oimodeler.oimData import oimData
from oimodeler.oimModel import oimModel
from oimodeler.oimPlots import oimWlTemplatePlots
from oimodeler.oimSimulator import oimSimulator


@pytest.fixture
def data(global_data_dir: Path) -> oimData:
    """Read-in dataset."""
    return oimData(list((global_data_dir / "ASPRO_MATISSE2").glob("*.fits")))


@pytest.fixture
def wl() -> float:
    """Wavelength."""
    return 3.5e-6


@pytest.fixture
def spf(wl: float) -> np.ndarray:
    """Spatial frequencies."""
    B = np.linspace(0.0, 300, num=200)
    return B / wl


class TestOimModel:
    """Tests the plotting methods of the oimModel class."""

    @pytest.fixture
    def components(self) -> list[oimComponent]:
        """List of basic components."""
        return [
            oimPt(f=0.1),
            oimUD(d=10, f=0.5),
            oimGauss(fwhm=5, f=1),
            oimIRing(d=5, f=0.5),
        ]

    @pytest.fixture
    def models(self, components: list[oimComponent]) -> list[oimModel]:
        """List of models containing basic components."""
        return [oimModel(c) for c in components]

    def test_showModel(
        self, wl: float, spf: np.ndarray, models: list[oimModel]
    ) -> None:
        """Tests the showModel plotting method."""
        _, ax = plt.subplots(
            2,
            len(models),
            figsize=(3 * len(models), 6),
            sharex="row",
            sharey="row",
        )

        for i, model in enumerate(models):
            model.showModel(
                512, 0.1, normPow=0.2, axe=ax[0, i], colorbar=False
            )
            v = np.abs(model.getComplexCoherentFlux(spf, spf * 0))
            v = v / v.max()
            ax[1, i].plot(spf, v)
            ax[0, i].set_title(model.components[-1].name)
            ax[1, i].set_xlabel("Sp. freq. (cycles/rad)")

        plt.close()

    def test_showFourier(
        self, wl: float, spf: np.ndarray, models: list[oimModel]
    ) -> None:
        """Tests the showFourier plotting method."""
        _, ax = plt.subplots(
            2,
            len(models),
            figsize=(3 * len(models), 6),
            sharex="row",
            sharey="row",
        )
        for i, model in enumerate(models):
            model.showFourier(
                512,
                0.1,
                wl=wl,
                axe=ax[0, i],
                colorbar=False,
                display="amp",
            )
            model.showFourier(
                512,
                0.1,
                wl=wl,
                axe=ax[1, i],
                colorbar=False,
                display="phase",
            )
            ax[0, i].set_title(model.components[-1].name)

        plt.close()


class TestOimPlot:
    """Tests the plotting functionalities of the oimPlot module."""

    def test_uvplot(self, data: oimData) -> None:
        """Tests oimPlot.uvplot."""
        _, ax = plt.subplots(
            3, 2, subplot_kw={"projection": "oimAxes"}, figsize=(8, 10)
        )
        ax[0, 0].uvplot(
            data.data,
            color="byBaseline",
            marker=".",
            legendkwargs={"fontsize": 6},
        )
        ax[0, 1].uvplot(
            data.data,
            color="byFile",
            facecolor="w",
            legendkwargs={"fontsize": 3.8},
        )
        ax[1, 0].uvplot(
            data.data, color="byConfiguration", colorTab=["r", "g", "b"]
        )
        ax[1, 1].uvplot(
            data.data, color="byArrname", colorTab=["r", "g", "b"], marker="+"
        )
        ax[2, 0].uvplot(data.data, label="custom label", unit="km")
        ax[2, 1].uvplot(
            data.data,
            unit="cycle/rad",
            cunit="micron",
            lw=3,
            color="byConfiguration",
        )
        plt.close()

    @pytest.mark.parametrize("column", ["VISAMP", "VIS2DATA", "T3PHI"])
    def test_oiplot(self, data: oimData, column: str) -> None:
        """Tests oimPlot.oiplot."""
        _, (ax, bx) = plt.subplots(1, 2, subplot_kw={"projection": "oimAxes"})
        ax.oiplot(
            data,
            "EFF_WAVE",
            column,
            xunit="micron",
            color="byConfiguration",
            errorbar=True,
            kwargs_error={"alpha": 0.3},
        )
        ax.legend()

        bx.oiplot(
            data,
            "SPAFREQ",
            column,
            xunit="cycle/rad",
            errorbar=True,
            lw=2,
            ls=":",
            color="byFile",
        )
        bx.legend()
        plt.close()

    def test_template(self, data: oimData) -> None:
        """Tests template plots."""
        fig = plt.figure(FigureClass=oimWlTemplatePlots, figsize=(12, 7))
        fig.autoShape(
            data.data, shape=[["VIS2DATA", None], ["VISPHI", "T3PHI"]]
        )
        fig.set_xunit("micron")
        fig.plot(
            data.data,
            plotFunction=plt.Axes.errorbar,
            plotFunctionkwarg={"color": "gray", "alpha": 0.3},
        )
        fig.plot(data.data, plotFunctionkwarg={"color": "tab:blue", "lw": 0.5})
        fig.set_ylim(["VISPHI", "T3PHI"], -25, 25)
        fig.set_ylim(["VIS2DATA"], 0, 1.2)
        fig.set_xlim(3, 3.6)
        fig.set_legends(
            0.5,
            0.1,
            "$BASELINE$ $LENGTH$m $PA$$^o$",
            ["VIS2DATA", "VISPHI"],
            fontsize=12,
            ha="center",
        )
        fig.set_legends(
            0.5, 0.1, "$BASELINE$", ["T3PHI"], fontsize=12, ha="center"
        )
        fig.tight_layout()
        plt.close()


class TestOimSimulator:
    """Tests the plotting methods of the oimSimulator class."""

    @pytest.fixture
    def model(self) -> oimModel:
        return oimModel([oimPt(f=0.5), oimUD(d=3, f=1, x=10, y=20)])

    @pytest.fixture
    def sim(self, data: oimData, model: oimModel) -> oimSimulator:
        return oimSimulator(data=data, model=model)

    def test_plot(self, sim: oimSimulator) -> None:
        """Tests oimSimulator.plot."""
        sim.plot(["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"])
        plt.close()

    def test_plotResiduals(self, sim: oimSimulator) -> None:
        """Tests oimSimulator.plotResiduals."""
        sim.plotResiduals(["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"])
        plt.close()

    def test_plotWithResiduals(self, sim: oimSimulator) -> None:
        """Tests oimSimulator.plotResiduals."""
        sim.plotWithResiduals(
            ["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"]
        )
        plt.close()
