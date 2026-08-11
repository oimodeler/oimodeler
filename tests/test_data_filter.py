"""
Tests for the oimodeler.oimDataFilter module.
"""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import oifits
import pytest
from astropy.io import fits
from numpy.typing import ArrayLike, NDArray

import oimodeler.oimUtils as utils
from oimodeler.oimData import oimData
from oimodeler.oimDataFilter import (
    oimDataTypeFilter,
    oimKeepDataTypeFilter,
    oimRemoveArrayFilter,
    oimRemoveInsnameFilter,
    oimResetFlagsFilter,
    oimSetMinErrFilter,
    oimWavelengthRangeFilter,
)

# TODO: Potentially move the tested function from oimUtils to oimDataFilter module?
# Concerning those that find no usage in other scripts.

COLUMN_TO_ERROR: dict[str, str] = {
    "VISAMP": "VISAMPERR",
    "VISPHI": "VISPHIERR",
    "VIS2DATA": "VIS2ERR",
    "T3PHI": "T3PHIERR",
    "T3AMP": "T3AMPERR",
}

TABLE_TO_COLUMN: dict[str, list[str]] = {
    "OI_VIS": ["VISAMP", "VISPHI"],
    "OI_VIS2": ["VIS2DATA"],
    "OI_T3": ["T3PHI", "T3AMP"],
}
COLUMN_TO_TABLE: dict[str, str] = {
    "VISAMP": "OI_VIS",
    "VISPHI": "OI_VIS",
    "VIS2DATA": "OI_VIS2",
    "T3PHI": "OI_T3",
    "T3AMP": "OI_T3",
}


@pytest.fixture(scope="module")
def oifits_files(global_data_dir: Path) -> list[Path]:
    """List of OIFITS files."""
    return sorted((global_data_dir / "ASPRO_MATISSE2").glob("*.fits"))


@pytest.fixture(scope="module")
def data(oifits_files: list[Path]) -> oimData:
    """Read-in dataset."""
    return oimData(oifits_files)


@pytest.fixture(scope="module")
def oifits_inputs(
    oifits_files: list[Path], data: oimData
) -> list[Path | fits.HDUList]:
    """A list containing either paths to OIFITS files or the
    read-in files themselves."""
    return [oifits_files[0], *data.data]


@pytest.fixture(scope="module")
def wl() -> NDArray[np.float32]:
    """Wavelengths."""
    return np.array(
        [
            2.98966e-06,
            3.01416e-06,
            3.04099e-06,
            3.25596e-06,
            3.29906e-06,
            3.33992e-06,
            3.41283e-06,
            3.45752e-06,
            3.49175e-06,
        ],
        dtype=np.float32,
    )


@pytest.fixture(scope="module")
def hdu(wl: NDArray[np.float32]) -> fits.BinTableHDU:
    """Create a sample fits.BinTableHDU of the "OI_VIS" extension."""
    shape = (6, len(wl))
    visamp = np.random.rand(*shape) * 3
    visamperr = np.random.rand(*shape)
    visphi = np.random.rand(*shape) * 180
    visphierr = np.random.rand(*shape) * 20

    columns = [
        fits.Column(
            name="VISAMP", format=f"{shape[-1]}E", unit="Jy", array=visamp
        ),
        fits.Column(
            name="VISAMPERR",
            format=f"{shape[-1]}E",
            unit="Jy",
            array=visamperr,
        ),
        fits.Column(
            name="VISPHI", format=f"{shape[-1]}E", unit="deg", array=visphi
        ),
        fits.Column(
            name="VISPHIERR",
            format=f"{shape[-1]}E",
            unit="deg",
            array=visphierr,
        ),
        fits.Column(
            name="FLAG",
            format=f"{shape[-1]}L",
            unit="",
            array=~visamp.astype(bool),
        ),
    ]

    hdu = fits.BinTableHDU.from_columns(fits.ColDefs(columns))
    hdu.header["EXTNAME"] = "OI_VIS"
    return hdu


class TestOimDataFilterComponent:
    """Tests the oimDataFilter.oimDataFilterComponent class."""


class TestOimDataFilter:
    """Tests the oimDataFilter.TestOimDataFilter class."""


@pytest.mark.parametrize(
    "arr", (["OI_VIS"], ["OI_VIS2"], ["OI_T3"], ["OI_VIS2", "OI_T3"])
)
def test_oimRemoveArrayFilter(data: oimData, arr: list[str]) -> None:
    """Tests the oimDataFilter.oimRemoveArrayFilter class."""
    data = copy.deepcopy(data)
    data.setFilter(oimRemoveArrayFilter(targets="all", arr=arr))

    for hdul in data.data:
        oif = oifits.open(copy.deepcopy(hdul), quiet=True)
        assert oif.wavelength

        for table in ["OI_VIS", "OI_VIS2", "OI_T3"]:
            oif_table = getattr(oif, table.split("_")[-1].lower())
            if table in arr:
                assert oif_table.size == 0
            else:
                assert oif_table.size != 0


def test_oimRemoveInsnameFilter(data: oimData) -> None:
    """Tests the oimDataFilter.oimRemoveInsnameFilter class."""
    data = copy.deepcopy(data)
    data.setFilter(
        oimRemoveInsnameFilter(
            targets="all", insname="MATISSE_LM_2.86542-4.18239-62ch"
        )
    )

    for hdul in data.data:
        oif = oifits.open(copy.deepcopy(hdul), quiet=True)
        assert not oif.wavelength

        for table in ["vis", "vis2", "t3"]:
            assert getattr(oif, table).size == 0


@pytest.mark.parametrize(
    "dataType",
    (
        ["VISAMP"],
        ["VISPHI"],
        ["VIS2DATA"],
        ["T3PHI"],
        ["VISAMP", "T3PHI"],
        ["T3AMP"],
        ["T3PHI", "T3AMP"],
    ),
)
def test_oimDataTypeFilter(data: oimData, dataType: list[str]) -> None:
    """Tests the oimDataFilter.oimDataTypeFilter class."""
    data = copy.deepcopy(data)
    data.setFilter(oimDataTypeFilter(targets="all", dataType=dataType))

    for hdul in data.data:
        oif = oifits.open(copy.deepcopy(hdul), quiet=True)
        assert oif.wavelength

        for table in ["vis", "vis2", "t3"]:
            for elem in getattr(oif, table):
                for dtype in dataType:
                    if f"_{dtype.lower()}" in dir(elem):
                        assert not getattr(elem, dtype.lower()).any()
                        assert getattr(
                            elem, COLUMN_TO_ERROR[dtype].lower()
                        ).any()


@pytest.mark.parametrize(
    "dataType",
    (
        ["VISAMP"],
        ["VISPHI"],
        ["VIS2DATA"],
        ["T3PHI"],
        ["VISAMP", "T3PHI"],
        ["T3AMP"],
        ["T3PHI", "T3AMP"],
    ),
)
def test_oimKeepDataTypeFilter(data: oimData, dataType: list[str]) -> None:
    """Tests the oimDataFilter.oimKeepDataTypeFilter class."""
    data = copy.deepcopy(data)
    data.setFilter(oimKeepDataTypeFilter(targets="all", dataType=dataType))

    for hdul in data.data:
        oif = oifits.open(copy.deepcopy(hdul), quiet=True)
        assert oif.wavelength

        dtypes = [COLUMN_TO_TABLE[dtype].lower() for dtype in dataType]
        for table in ["vis", "vis2", "t3"]:
            oif_table = getattr(oif, table)
            if f"oi_{table}" not in dtypes:
                assert oif_table.size == 0

            for elem in oif_table:
                for column in TABLE_TO_COLUMN[f"oi_{table}".upper()]:
                    if column in dataType:
                        assert getattr(elem, column.lower()).any()
                        assert getattr(
                            elem, COLUMN_TO_ERROR[column].lower()
                        ).any()
                    else:
                        assert not getattr(elem, column.lower()).any()
                        assert getattr(
                            elem, COLUMN_TO_ERROR[column].lower()
                        ).any()


class TestFlagWithExpressionFilter:
    """Tests all functionality related to the expression-based filtering."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oifitsFlagWithExpression(self) -> None:
        """Tests the oimUtils.oifitsFlagWithExpression function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oimFlagWithExpressionFilter(self) -> None:
        """Tests the oimDataFilter.oimFlagWithExpressionFilter class."""


class TestWavelengthFilters:
    """Tests all functionality related to wavelength filtering."""

    # TODO: Further test kwargs: addCut
    @pytest.mark.parametrize(
        "wlRange", ([3e-6, 4e-6], [3.5e-6, 4e-6], [3e-6, 3.5e-6])
    )
    def test_cutWavelengthRange(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        wlRange: list[float],
    ) -> None:
        """Tests the oimUtils.cutWavelengthRange."""
        for data in oifits_inputs:
            oif = oifits.open(copy.deepcopy(data), quiet=True)
            oif_filt = oifits.open(
                utils.cutWavelengthRange(copy.deepcopy(data), wlRange=wlRange),
                quiet=True,
            )

            for table in ["vis", "vis2", "flux", "t3"]:
                oif_table = getattr(oif, table)
                oif_filt_table = getattr(oif_filt, table)

                for column, column_filt in zip(oif_table, oif_filt_table):
                    assert (
                        column.wavelength.eff_wave.size
                        > column_filt.wavelength.eff_wave.size
                    )
                    assert np.all(
                        (column_filt.wavelength.eff_wave >= wlRange[0])
                        & (column_filt.wavelength.eff_wave <= wlRange[-1])
                    )

        # TODO: Further test kwargs: addCut, method
        @pytest.mark.parametrize(
            "wlRange", ([3e-6, 4e-6], [3.5e-6, 4e-6], [3e-6, 3.5e-6])
        )
        def test_oimWavelengthRangeFilter(
            self, data: oimData, wlRange: list[float]
        ) -> None:
            """Tests the oimDataFilter.oimWavelengthRangeFilter class."""
            oif = oifits.open(copy.deepcopy(data), quiet=True)

            data = copy.deepcopy(data)
            data.setFilter(
                oimWavelengthRangeFilter(targets="all", wlRange=wlRange)
            )

            oif_filt = oifits.open(copy.deepcopy(data), quiet=True)
            for table in ["vis", "vis2", "flux", "t3"]:
                oif_table = getattr(oif, table)
                oif_filt_table = getattr(oif_filt, table)

                for column, column_filt in zip(oif_table, oif_filt_table):
                    assert (
                        column.wavelength.eff_wave.size
                        > column_filt.wavelength.eff_wave.size
                    )
                    assert np.all(
                        (column_filt.wavelength.eff_wave >= wlRange[0])
                        & (column_filt.wavelength.eff_wave <= wlRange[-1])
                    )

        @pytest.mark.skip(reason="Doesn't currently work?")
        def test_shiftWavelength(self) -> None:
            """Tests the oimUtils.shiftWavelength function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_oimWavelengthShiftFilter(self) -> None:
            """Tests the oimDataFilter.oimWavelengthShiftFilter class."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_spectralSmoothing(self) -> None:
            """Tests the oimUtils.spectralSmoothing function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_oimWavelengthSmoothingFilter(self) -> None:
            """Tests the oimDataFilter.oimWavelengthSmoothingFilter class."""

    class TestBinningFilter:
        """Tests all functionality related to the
        oimDataFilter.oimWavelengthBinningFilter class."""

        # NOTE: Currently kind="mean" and kind="circular" are identical.
        # Change test parameters when the underlying function changes.
        @pytest.mark.parametrize(
            "array,binSize,kind,expected",
            [
                (np.array([1, 2, 3, 4]), 1, "mean", [1.0, 2.0, 3.0, 4.0]),
                (np.array([1, 2, 3, 4]), 1, "median", [1.0, 2.0, 3.0, 4.0]),
                (np.array([1, 2, 3, 4]), 1, "circular", [1.0, 2.0, 3.0, 4.0]),
                (np.array([1, 2, 100, 4, 5, 6]), 2, "mean", [1.5, 52.0, 5.5]),
                (
                    np.array([1, 2, 100, 4, 5, 6]),
                    2,
                    "median",
                    [1.5, 52.0, 5.5],
                ),
                (
                    np.array([1, 2, 100, 4, 5, 6]),
                    2,
                    "circular",
                    [1.5, 52.0, 5.5],
                ),
                (np.array([1, 2, 6, 4, 5, 9]), 3, "mean", [3.0, 6.0]),
                (np.array([1, 2, 6, 4, 5, 9]), 3, "median", [2.0, 5.0]),
                (np.array([1, 2, 6, 4, 5, 9]), 3, "circular", [3.0, 6.0]),
            ],
        )
        def test__rebin(
            self,
            array: np.ndarray,
            binSize: int,
            kind: str,
            expected: list[float],
        ) -> None:
            """Tests the oimUtils._rebin."""
            assert np.array_equal(
                utils._rebin(array, binSize=binSize, kind=kind), expected
            )

        # NOTE: Currently kind="mean" and kind="circular" are identical.
        # Change test parameters when the underlying function changes.
        @pytest.mark.skip(reason="Test not yet finished.")
        def test__rebinHDU(
            self, hdu: fits.BinTableHDU, binSize: int, exception: list[str]
        ) -> None:
            """Tests the oimUtils._rebinHDU function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_binWavelength(self) -> None:
            """Tests the oimUtils.binWavelength function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_oimWavelengthBinninFilter(self) -> None:
            """Tests the oimDataFilter.oimWavelengthBinningFilter class."""

    class TestIntpBinFilter:
        """Tests all functionality related to the
        oimDataFilter.oimWavelengthIntpBinFilter class."""

        # NOTE: Currently kind="mean" and kind="circular" are identical.
        # Change test parameters when the underlying function changes.
        # TODO: Add kwarg tests: values, and nSpecChannels
        @pytest.mark.parametrize(
            "array,binMasks,binEdgeValues,kind,expected",
            (
                (
                    np.array([1, 2, 5, 4]),
                    [[True, True, True, False], [False, False, True, True]],
                    [[0, 0], [0, 0]],
                    "mean",
                    np.array([1.6, 2.25]),
                ),
                (
                    np.array([1, 2, 5, 4]),
                    [[True, True, True, False], [False, False, True, True]],
                    [[0.5, 2.5], [1.5, 3.5]],
                    "mean",
                    np.array([2.2, 3.5]),
                ),
                (
                    np.array([1, 2, 5, 4]),
                    [[True, True, True, False], [False, True, True, True]],
                    [[0, 3], [1, 6]],
                    "median",
                    np.array([2.0, 4.0]),
                ),
                (
                    np.array([1, 2, 5, 4]),
                    [[True, True, True, False], [False, False, True, True]],
                    [[0.5, 2.5], [1.5, 3.5]],
                    "circular",
                    np.array([2.2, 3.5]),
                ),
            ),
        )
        def test__intpBinning(
            self,
            array: np.ndarray,
            binMasks: ArrayLike,
            binEdgeValues: ArrayLike,
            kind: str,
            expected: list[float],
        ) -> None:
            """Tests the oimUtils._intpBinning function."""
            assert np.allclose(
                utils._intpBinning(
                    array,
                    binMasks,
                    binEdgeValues,
                    kind=kind,
                ),
                expected,
            )

        @pytest.mark.parametrize(
            "exception",
            (
                [],
                ["FLAG"],
                ["VISAMP", "VISAMPERR", "FLAG"],
                ["VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR", "FLAG"],
            ),
        )
        def test__interpolateBinHDU(
            self,
            wl: NDArray[np.float32],
            hdu: fits.BinTableHDU,
            exception: list[str],
        ) -> None:
            """Tests the oimUtils._interpolateBinHDU function."""
            hdu = copy.deepcopy(hdu)
            binGrid = np.array([3e-6, 3.3e-6, 3.45e-6])
            binEdgeValues = [
                (grid - 0.05e-6, grid + 0.05e-6) for grid in binGrid
            ]
            binMasks = np.array(
                [
                    (wl >= lower) & (wl <= upper)
                    for (lower, upper) in binEdgeValues
                ]
            )
            result_hdu = utils._interpolateBinHDU(
                hdu, binGrid, binMasks, binEdgeValues, wl, exception
            )
            for column in hdu.columns.names:
                orig_col = hdu.data[column]
                res_col = result_hdu.data[column]
                if column in exception:
                    assert np.array_equal(res_col, orig_col)
                    assert res_col.shape == orig_col.shape
                else:
                    assert res_col.shape != orig_col.shape
                    assert res_col.shape == (orig_col.shape[0], binGrid.size)

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_intpBinWavelength(self, data: oimData) -> None:
            """Tests the oimUtils.intpBinWavelength function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_oimWavelengthIntpBinFilter(self, data: oimData) -> None:
            """Tests the oimUtils.oimWavelengthIntpBinFilter class."""


class TestBaselinesFilters:
    """Tests all functionality related to the baseline filtering."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oifitsKeepBaselines(self) -> None:
        """Tests the oimUtils.oifitsKeepBaselines function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oimKeepBaselinesFilter(self) -> None:
        """Tests the oimDataFilter.test_oimKeepBaselinesFilter function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oifitsRemoveBaselines(self) -> None:
        """Tests the oimUtils.oifitsKeepBaselines function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oimRemoveBaselinesFilter(self) -> None:
        """Tests the oimDataFilter.test_oimRemoveBaselinesFilter function."""


class TestTelescopesFilters:
    """Tests all functionality related to the telescope filtering."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oifitsKeepTelescopes(self) -> None:
        """Tests the oimUtils.oifitsKeepTelescopes function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oimKeepTelescopesFilter(self) -> None:
        """Tests the oimDataFilter.oimKeepTelescopesFilter class."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oifitsRemoveTelescopes(self) -> None:
        """Tests the oimUtils.oifitsRemoveTelescopes function."""

    @pytest.mark.skip(reason="Test not yet finished.")
    def test_oimRemoveTelescopesFilter(self) -> None:
        """Tests the oimDataFilter.oimRemoveTelescopesFilter class."""


@pytest.mark.parametrize(
    "arr",
    (
        ["OI_VIS"],
        ["OI_VIS2"],
        ["OI_T3"],
        ["OI_FLUX"],
        ["OI_VIS2", "OI_T3"],
        ["OI_VIS", "OI_FLUX", "OI_T3"],
    ),
)
def test_oimResetFlagsFilter(data: oimData, arr: list[str]) -> None:
    """Tests the oimDataFilter.oimResetFlagsFilter class."""
    data = copy.deepcopy(data)
    data.setFilter(oimResetFlagsFilter(targets="all", arr=arr))

    for hdul in data.data:
        oif = oifits.open(copy.deepcopy(hdul), quiet=True)
        assert oif.wavelength

        for table in arr:
            for elem in getattr(oif, table.split("_")[-1].lower()):
                assert not elem.flag.any()


class TestErrorFilters:
    """Tests all functionality related to filtering errors."""

    class TestDiffErrFilter:
        """Tests all functionality related to the
        oimDataFilter.oimDiffErrFilter class."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_computeDifferentialError(self) -> None:
            """Tests the oimUtils.computeDifferentialError function."""

        @pytest.mark.skip(reason="Test not yet finished.")
        def test_oimDiffErrFilter(self) -> None:
            """Tests the oimUtils.oimDiffErrFilter class."""

    class TestSetMinErr:
        """Tests all functionality related to the
        oimDataFilter.oimSetMinErrFilter class."""

        @pytest.mark.parametrize(
            "dataTypes,values",
            (
                [["VISAMP"], [5]],
                [["VISAMP"], [10]],
                [["VISPHI"], [5]],
                [["VISPHI"], [10]],
                [["VIS2DATA"], [5]],
                [["VIS2DATA"], [10]],
                [["T3PHI"], [5]],
                [["T3PHI"], [10]],
                [["VISAMP", "T3PHI"], [5, 5]],
                [["VISAMP", "T3PHI"], [10, 10]],
                [["T3AMP"], [5]],
                [["T3AMP"], [10]],
                [["T3PHI", "T3AMP"], [5, 5]],
                [["T3PHI", "T3AMP"], [10, 10]],
            ),
        )
        def test_setMinimumError(
            self,
            data: oimData,
            dataTypes: list[str],
            values: list[float],
        ) -> None:
            """Tests the oimUtils.oimSetMiniumError function."""
            for hdul in copy.deepcopy(data).data:
                utils.setMinimumError(hdul, dataTypes, values)
                oif = oifits.open(copy.deepcopy(hdul), quiet=True)
                for column, value in zip(dataTypes, values):
                    table = COLUMN_TO_TABLE[column].split("_")[-1].lower()
                    err_name = COLUMN_TO_ERROR[column].lower()

                    for elem in getattr(oif, table):
                        err = getattr(elem, err_name).data
                        fraction = np.abs(
                            err / getattr(elem, column.lower()).data
                        )
                        if "PHI" in column:
                            assert np.all(
                                (err >= value) | np.isclose(err, value)
                            )
                        else:
                            assert np.all(
                                (fraction >= (value / 100))
                                | np.isclose(fraction, value / 100)
                            )

        @pytest.mark.parametrize(
            "dataTypes,values",
            (
                [["VISAMP"], [5]],
                [["VISAMP"], [10]],
                [["VISPHI"], [5]],
                [["VISPHI"], [10]],
                [["VIS2DATA"], [5]],
                [["VIS2DATA"], [10]],
                [["T3PHI"], [5]],
                [["T3PHI"], [10]],
                [["VISAMP", "T3PHI"], [5, 5]],
                [["VISAMP", "T3PHI"], [10, 10]],
                [["T3AMP"], [5]],
                [["T3AMP"], [10]],
                [["T3PHI", "T3AMP"], [5, 5]],
                [["T3PHI", "T3AMP"], [10, 10]],
            ),
        )
        def test_oimSetMinErr(
            self, data: oimData, dataTypes: list[str], values: list[float]
        ) -> None:
            """Tests the oimDataFilter.setMinErrFilter class."""
            data = copy.deepcopy(data)
            data.setFilter(
                oimSetMinErrFilter(
                    targets="all", dataType=dataTypes, values=values
                )
            )

            for hdul in data.data:
                oif = oifits.open(copy.deepcopy(hdul), quiet=True)
                for column, value in zip(dataTypes, values):
                    table = COLUMN_TO_TABLE[column].split("_")[-1].lower()
                    err_name = COLUMN_TO_ERROR[column].lower()

                    for elem in getattr(oif, table):
                        err = getattr(elem, err_name).data
                        fraction = np.abs(
                            err / getattr(elem, column.lower()).data
                        )
                        if "PHI" in column:
                            assert np.all(
                                (err >= value) | np.isclose(err, value)
                            )
                        else:
                            assert np.all(
                                (fraction >= (value / 100))
                                | np.isclose(fraction, value / 100)
                            )
