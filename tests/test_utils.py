"""
Tests for the oimodeler.oimUtils module.
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Iterable
from itertools import permutations
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import oifits
import pytest
from astropy.io import fits
from astropy.modeling.physical_models import BlackBody

import oimodeler.oimUtils as utils
from oimodeler.oimData import oimData
from oimodeler.oimOptions import constants as const


def generic_method(self=None) -> str:
    """A test method."""
    return f"I am now attached to {self.__class__.__name__}"


def unique_permutations(
    iterable: Iterable, length: int
) -> list[tuple[Any, ...]]:
    """Returns unique permutations of a certain length."""
    return list(set(permutations(set(iterable), length)))


@pytest.fixture(scope="module")
def oifits_files(global_data_dir: Path) -> list[Path]:
    """List of OIFITS files."""
    return sorted((global_data_dir / "ASPRO_MATISSE2").glob("*.fits"))


@pytest.fixture(scope="module")
def data(oifits_files: list[Path]) -> oimData:
    """Read-in dataset."""
    return oimData(oifits_files)


@pytest.fixture(scope="module")
def toml_file(tmp_path_factory: Path) -> Path:
    """Create a temp TOML file that gets cleaned up automatically."""
    tmp_dir = tmp_path_factory.mktemp("data")
    toml_path = tmp_dir / "parameters.toml"
    toml_path.write_text(
        """
        [default]
        name = ""
        value = 0
        mini = "-inf"
        maxi = "inf"
        description = ""
        unit = ""
        free = true
        error = 0
        """,
        encoding="utf-8",
    )
    return toml_path


@pytest.mark.skip(reason="Test not yet implemented.")
def test__pickle() -> None:
    """Tests the oimUtils._pickel function."""


@pytest.mark.skip(reason="Test not yet implemented.")
def test__unpickle() -> None:
    """Tests the oimUtils._unpickel function."""


@pytest.mark.skip(reason="Test not yet implemented.")
def test__serialize_function() -> None:
    """Tests the oimUtils._serialize_function function."""


@pytest.mark.skip(reason="Test not yet implemented.")
def test__deserialize_function() -> None:
    """Tests the oimUtils._deserialize_function function."""


def test_load_toml(toml_file: Path) -> None:
    """Tests the oimUtils.load_toml function."""
    param_dict = utils.load_toml(toml_file)
    assert "default" in param_dict

    param = param_dict["default"]
    assert not param["name"]
    assert param["value"] == 0
    assert param["mini"] == -np.inf
    assert param["maxi"] == np.inf
    assert not param["description"]
    assert param["unit"] == u.one
    assert param["free"]
    assert param["error"] == 0


@pytest.mark.parametrize(
    "method",
    (generic_method, [generic_method], {"generic_method": generic_method}),
)
def test_attach_methods(
    method: Callable | list | dict[str, Callable],
) -> None:
    """Tests the oimUtils.attach_methods function."""

    @utils.attach_methods(method)
    class GenericClass:
        pass

    assert "generic_method" in dir(GenericClass)

    test = GenericClass()
    assert "I am now attached to GenericClass" == test.generic_method()


class TestComputations:
    """Tests computation functions contained in the oimUtils module."""

    @pytest.mark.parametrize(
        "wl", (2.1e-6, 3.5e-6, 10.5e-6, np.array((2.1e-6, 3.5e-6, 10.5e-6)))
    )
    @pytest.mark.parametrize(
        "T", (1500, 3000, 5000, np.array((1500, 3000, 5000)))
    )
    def test_blackbody(
        self, wl: float | np.ndarray, T: float | np.ndarray
    ) -> None:
        """Tests the oimUtils.blackbody function."""
        expected = BlackBody(T * u.K)(wl * u.m).value
        assert expected == pytest.approx(utils.blackbody(T, const.c / wl))

    # FIXME: Currently only tests if there are errors in function
    # execution. Needs to actually test something more.
    def test_spectral_index(self, data: oimData) -> None:
        """Tests the oimUtils.compute_photometric_slope function."""
        utils.spectral_index(data, 6500)


class TestImageOperations:
    """Tests image altering functions of the oimUtils module."""

    @pytest.fixture(scope="module")
    def image(self) -> np.ndarray:
        """An image with only "hot" pixels."""
        return np.ones((1, 1, 8, 8))

    # TODO: Test more edge cases here? For instance multiple nt or nwl?
    # Or non-integer factors for the correct error raised? Or non-square image?
    @pytest.mark.parametrize("factor", (2, 4, 6, 8))
    def test_pad(self, image: np.ndarray, factor: int) -> None:
        """Tests the oimUtils.pad_image function."""
        im = utils.pad_image(image, factor)
        dim_x, dim_y = image.shape[-2], image.shape[-1]

        new_dim_x = 2 ** (dim_x * factor - 1).bit_length()
        new_dim_y = 2 ** (dim_y * factor - 1).bit_length()
        expected_shape = (1, 1, new_dim_x, new_dim_y)
        assert im.shape == expected_shape

        centre_x = (new_dim_x - dim_x) // 2
        centre_y = (new_dim_y - dim_y) // 2
        padded_region = im[
            ..., centre_x : centre_x + dim_x, centre_y : centre_y + dim_y
        ]
        assert np.allclose(padded_region, image)

        top_padding = im[..., 0:centre_x, :]
        assert np.all(top_padding == 0)

        bottom_padding = im[..., centre_x + dim_x :, :]
        assert np.all(bottom_padding == 0)

        left_padding = im[..., :, 0:centre_y]
        assert np.all(left_padding == 0)

        right_padding = im[..., :, centre_y + dim_y :]
        assert np.all(right_padding == 0)


# TODO: Extend coverage UT configurations and maybe multiple
# files simultaneously?
class TestOIFITSOperations:
    """Tests readout functions contained in the oimUtils module."""

    @pytest.fixture(scope="module")
    def oifits_inputs(
        self, oifits_files: list[Path], data: oimData
    ) -> list[Path | fits.HDUList]:
        """A list containing either paths to OIFITS files or the
        read-in files themselves."""
        return [oifits_files[0], *data.data]

    @pytest.mark.parametrize(
        "dataType",
        (
            "VIS2DATA",
            "VISAMP",
            "VISPHI",
            "T3AMP",
            "T3PHI",
            "FLUXDATA",
            "",
            None,
        ),
    )
    def test_getDataTypeIsAnalysisComplex(self, dataType: str | None) -> None:
        """Tests the oimUtils.getDataTypeIsAnalysisComplex function."""
        if not dataType or dataType is None:
            with pytest.raises(
                TypeError, match=f"{dataType} not a valid OIFITS2 datatype"
            ):
                utils.getDataTypeIsAnalysisComplex(dataType)
        else:
            iscomplex = utils.getDataTypeIsAnalysisComplex(dataType)
            if "PHI" in dataType:
                assert iscomplex
            else:
                assert not iscomplex

    @pytest.mark.parametrize(
        "dataType",
        (
            "VIS2DATA",
            "VISAMP",
            "VISPHI",
            "T3AMP",
            "T3PHI",
            "FLUXDATA",
            "",
            None,
        ),
    )
    def test_getDataArrname(self, dataType: str | None) -> None:
        """Tests the oimUtils.getDataArrname function."""
        if not dataType or dataType is None:
            with pytest.raises(
                TypeError, match=f"{dataType} not a valid OIFITS2 datatype"
            ):
                utils.getDataArrname(dataType)
        else:
            arr_name = utils.getDataArrname(dataType)

            expected = ""
            if "VIS2" in dataType:
                expected = "OI_VIS2"
            elif "VIS" in dataType:
                expected = "OI_VIS"
            elif "T3" in dataType:
                expected = "OI_T3"
            elif "FLUX" in dataType:
                expected = "OI_FLUX"

            assert arr_name == expected

    @pytest.mark.parametrize(
        "dataArrname", ("OI_VIS2", "OI_VIS", "OI_T3", "OI_FLUX", "", None)
    )
    def test_getDataType(self, dataArrname: str | None) -> None:
        """Tests the oimUtils.getDataType function."""
        dataType = set(utils.getDataType(dataArrname))
        if not dataArrname or dataArrname is None:
            assert not dataType
        else:
            expected = ""
            if "VIS2" in dataArrname:
                expected = {"VIS2DATA"}
            elif "VIS" in dataArrname:
                expected = {"VISAMP", "VISPHI"}
            elif "T3" in dataArrname:
                expected = {"T3PHI", "T3AMP"}
            elif "FLUX" in dataArrname:
                expected = {"FLUXDATA"}

            assert dataType == expected

    @pytest.mark.parametrize(
        "dataArrname", ("OI_VIS2", "OI_VIS", "OI_T3", "OI_FLUX", "", None)
    )
    def test_getDataTypeError(self, dataArrname: str | None) -> None:
        """Tests the oimUtils.getDataTypeError function."""
        dataTypeError = set(utils.getDataTypeError(dataArrname))
        if not dataArrname or dataArrname is None:
            assert not dataTypeError
        else:
            expected = ""
            if "VIS2" in dataArrname:
                expected = {"VIS2ERR"}
            elif "VIS" in dataArrname:
                expected = {"VISAMPERR", "VISPHIERR"}
            elif "T3" in dataArrname:
                expected = {"T3PHIERR", "T3AMPERR"}
            elif "FLUX" in dataArrname:
                expected = {"FLUXERR"}

            assert dataTypeError == expected

    # TODO: Add kwarg tests: extver, and squeeze.
    @pytest.mark.parametrize("arr", ("OI_VIS2", "OI_VIS", "OI_T3"))
    @pytest.mark.parametrize("length", (False, True))
    @pytest.mark.parametrize("angle", (False, True))
    def test_getBaselineName(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        arr: str,
        length: bool,
        angle: bool,
    ) -> None:
        """Tests the oimUtils.test_getBaselineName function."""
        for data in oifits_inputs:
            oif_arr = getattr(
                oifits.open(copy.deepcopy(data), quiet=True),
                arr.split("_")[-1].lower(),
            )

            expected = []
            for elem in oif_arr:
                baseline = "-".join([e.sta_name for e in elem.station])

                if "VIS" in arr:
                    if length:
                        baseline += (
                            f" {np.hypot(elem.ucoord, elem.vcoord):.0f}m"
                        )
                    if angle:
                        baseline += f" {np.rad2deg(np.arctan2(elem.ucoord, elem.vcoord)):.0f}$^o$"

                expected.append(baseline)

            assert (
                utils.getBaselineName(data, arr, length=length, angle=angle)
                == expected
            )

    # TODO: Add kwarg tests: extver.
    @pytest.mark.parametrize("arr", ("OI_VIS2", "OI_VIS", "OI_T3"))
    @pytest.mark.parametrize("squeeze", (False, True))
    def test_getConfigName(
        self, oifits_inputs: list[Path | fits.HDUList], arr: str, squeeze: bool
    ) -> None:
        """Tests the oimUtils.test_getConfigName function."""
        for data in oifits_inputs:
            oif = oifits.open(copy.deepcopy(data), quiet=True)
            oif_arr = getattr(oif, arr.split("_")[-1].lower())
            sta_names = {e.sta_name for elem in oif_arr for e in elem.station}

            expected = "-".join(
                [
                    elem.sta_name
                    for elem in oif.array["VLTI"].station
                    if elem.sta_name in sta_names
                ]
            )
            if not squeeze:
                expected = [expected]

            assert utils.getConfigName(data, arr, squeeze=squeeze) == expected

    # TODO: Add kwarg tests: extver, squeeze, and showFlagged.
    # FIXME: Check if this test accounts for multiple extensions of the same name
    @pytest.mark.parametrize("arr", ("OI_VIS2", "OI_VIS", "OI_T3"))
    @pytest.mark.parametrize("T3Max", (False, True))
    @pytest.mark.parametrize("returnUV", (False, True))
    def test_getBaselineLengthAndPA(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        arr: str,
        T3Max: bool,
        returnUV: bool,
    ) -> None:
        """Tests the oimUtils.test_getBaselineLengthAndPA function."""
        for data in oifits_inputs:
            oif_arr = getattr(
                oifits.open(copy.deepcopy(data), quiet=True),
                arr.split("_")[-1].lower(),
            )

            if "VIS" in arr:
                expected_uvcoord = np.transpose(
                    [(elem.ucoord, elem.vcoord) for elem in oif_arr]
                )
            else:
                expected_ucoord = np.transpose(
                    [
                        (
                            elem.u1coord,
                            elem.u2coord,
                            -(elem.u1coord + elem.u2coord),
                        )
                        for elem in oif_arr
                    ]
                )
                expected_vcoord = np.transpose(
                    [
                        (
                            elem.v1coord,
                            elem.v2coord,
                            -(elem.v1coord + elem.v2coord),
                        )
                        for elem in oif_arr
                    ]
                )
                expected_uvcoord = np.array([expected_ucoord, expected_vcoord])

            expected_baselines = np.hypot(*expected_uvcoord)
            expected_angles = np.rad2deg(np.arctan2(*expected_uvcoord))
            if "T3" in arr:
                expected_angles = np.max(expected_angles, axis=0)
                if T3Max:
                    expected_baselines = np.max(expected_baselines, axis=0)

                    # TODO: Change this when PAs of triangle are
                    # determined in tested function.
                    expected_angles *= 0

            result = utils.getBaselineLengthAndPA(
                data, arr, T3Max=T3Max, returnUV=returnUV
            )

            baselines, angles = result[:2]
            assert np.array_equal(baselines, expected_baselines)
            assert np.allclose(angles, expected_angles)

            if returnUV:
                ucoord, vcoord = result[2:]
                expected_ucoord, expected_vcoord = expected_uvcoord
                assert np.array_equal(ucoord, expected_ucoord)
                assert np.array_equal(vcoord, expected_vcoord)

    # TODO: Further test kwargs: extver.
    @pytest.mark.parametrize("arr", ("OI_VIS2", "OI_VIS", "OI_T3"))
    @pytest.mark.parametrize("squeeze", (False, True))
    @pytest.mark.parametrize(
        "unit,factor",
        (
            (None, 1),
            ("cycles/mas", u.mas.to(u.rad)),
            ("cycles/arcsec", u.arcsec.to(u.rad)),
            ("Mlam", 1e-6),
        ),
    )
    def test_getSpaFreq(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        arr: str,
        squeeze: bool,
        unit: str | None,
        factor: float,
    ) -> None:
        """Tests the oimUtils.test_getSpaFreq function."""
        for data in oifits_inputs:
            oif_arr = getattr(
                oifits.open(copy.deepcopy(data), quiet=True),
                arr.split("_")[-1].lower(),
            )
            baselines, _ = utils.getBaselineLengthAndPA(
                data, arr, squeeze=False, T3Max="T3" in arr
            )
            expected = [
                baseline / elem.wavelength.eff_wave * factor
                for baseline, elem in zip(baselines[0], oif_arr)
            ]

            assert np.allclose(
                utils.getSpaFreq(data, arr, squeeze=squeeze, unit=unit),
                expected,
            )

    # TODO: Further test kwargs: extver.
    @pytest.mark.parametrize("arr", ["OI_VIS2", "OI_VIS", "OI_T3"])
    @pytest.mark.parametrize("squeeze", (False, True))
    @pytest.mark.parametrize(
        "unit,factor",
        (
            (None, 1),
            ("cycles/mas", u.mas.to(u.rad)),
            ("cycles/arcsec", u.arcsec.to(u.rad)),
            ("Mlam", 1e-6),
        ),
    )
    def test_get2DSpaFreq(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        arr: str,
        squeeze: bool,
        unit: str,
        factor: float,
    ) -> None:
        """Tests the oimUtils.test_get2DSpaFreq function."""
        for data in oifits_inputs:
            if "T3" in arr:
                with pytest.raises(
                    TypeError,
                    match="get2DSpaFreq does not accept OI_T3 extension",
                ):
                    utils.get2DSpaFreq(data, arr, unit=unit, squeeze=squeeze)
            else:

                oif_arr = getattr(
                    oifits.open(copy.deepcopy(data), quiet=True),
                    arr.split("_")[-1].lower(),
                )
                *_, ucoord, vcoord = utils.getBaselineLengthAndPA(
                    data, arr, squeeze=False, returnUV=True
                )

                # HACK: The reshape here is a bit arbitrary right now.
                # FIXME: Understand the proper reshape in the function.
                expected_spaFreqU = np.array(
                    [
                        baseline / elem.wavelength.eff_wave * factor
                        for baseline, elem in zip(ucoord[0], oif_arr)
                    ],
                )[np.newaxis, ...]
                expected_spaFreqV = np.array(
                    [
                        baseline / elem.wavelength.eff_wave * factor
                        for baseline, elem in zip(vcoord[0], oif_arr)
                    ],
                )[np.newaxis, ...]
                expected = (
                    expected_spaFreqU[0] if squeeze else expected_spaFreqU,
                    expected_spaFreqV[0] if squeeze else expected_spaFreqV,
                )

                assert np.allclose(
                    utils.get2DSpaFreq(data, arr, unit=unit, squeeze=squeeze),
                    expected,
                )

    # TODO: Further test kwargs: extver
    @pytest.mark.parametrize("returnBand", (False, True))
    @pytest.mark.parametrize("arr", ("OI_VIS2", "OI_VIS", "OI_T3"))
    def test_getWlFromOifits(
        self,
        oifits_inputs: list[Path | fits.HDUList],
        arr: str,
        returnBand: bool,
    ) -> None:
        """Tests the oimUtils.test_getWlFromOifits function."""
        for data in oifits_inputs:
            oif_wl = next(
                iter(
                    oifits.open(
                        copy.deepcopy(data), quiet=True
                    ).wavelength.values()
                )
            )
            expected_wavelength = oif_wl.eff_wave

            result = utils.getWlFromOifits(data, arr, returnBand=returnBand)
            if returnBand:
                wavelength, bandwidth = result
                expected_bandwidth = oif_wl.eff_band
                assert np.array_equal(bandwidth, expected_bandwidth)
            else:
                wavelength = result

            assert np.array_equal(wavelength, expected_wavelength)

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_getWlFromFitsImageCube(self) -> None: ...

    # TODO: Is this even needed or would copy.deepcopy suffice?
    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_hdulistDeepCopy(self) -> None: ...

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test__createOiTab(self) -> None: ...

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_createOiTargetFromSimbad(self) -> None: ...
