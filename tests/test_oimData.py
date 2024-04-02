from pathlib import Path

from astropy.io import fits

import oimodeler as oim


def test_oimDataGetWl(real_data_dir: Path) -> None:
    """Tests the oimDataGetWl function."""
    files = list((real_data_dir / "MATISSE" / "binary75Vir").glob("*"))
    # data = [fits.open(fits_file) for fits_file in files][0]
    # t = oim.oimDataGetWl(data, data[0], False)


def test_oimDataType() -> None:
    ...


def test_oimGetDataValErrAndTypeFlag() -> None:
    ...


def test_oimDataCheckData() -> None:
    ...


def test_oimDataGetVectCoord() -> None:
    ...


def test_oimData_init() -> None:
    ...


def test_oimData_str_() -> None:
    ...


def test_oimData_addData() -> None:
    ...


def test_oimData_removeData() -> None:
    ...


def test_oimData_setFilter() -> None:
    ...


def test_oimData_applyFilter() -> None:
    ...


def test_oimData_useFilter() -> None:
    ...


def test_oimData_analyzeOIFitFile() -> None:
    ...


def test_oimData_prepareData() -> None:
    ...


def test_oimData_writeto() -> None:
    ...
