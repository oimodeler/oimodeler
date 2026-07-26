import copy

import numpy as np
import pytest

from oimodeler.oimComponent import oimComponent
from oimodeler.oimModel import oimModel
from oimodeler.oimParam import oimParam, oimParamInterpolator


# TODO: Should __eq__ for oimParam, oimComponent, etc. be implemented?
def assert_param_equal(original: oimParam, restored: oimParam) -> bool:
    """Compares an original to a restored oimParam."""
    results = []
    for key, value_org in original.__dict__.items():
        value_res = getattr(restored, key)

        if key in ["value", "error"]:
            value_res = pytest.approx(value_res)
        elif key == "unit":
            value_org, value_res = (
                value_org.to_string(),
                value_res.to_string(),
            )

        results.append(value_org == value_res)

    result = all(results)
    assert result
    return result


def assert_intp_equal(
    original: oimParamInterpolator, restored: oimParamInterpolator
):
    """Compares an original to a restored oimParamInterpolator."""
    original_dict = copy.deepcopy(original.__dict__)
    for key, value in original_dict.items():
        value_res = getattr(restored, key)
        if isinstance(value, (tuple, list, np.ndarray)):
            assert isinstance(value_res, (tuple, list, np.ndarray))
            for i, v in enumerate(value):
                if isinstance(v, oimParam):
                    assert_param_equal(v, value_res[i])
                else:
                    assert v == value_res[i]

        elif isinstance(value, oimParam):
            assert_param_equal(value, value_res)
        else:
            assert value == value_res


def assert_component_equal(
    original: oimComponent, restored: oimComponent
) -> bool:
    """Compares an original to a restored oimComponent."""
    original_dict = copy.deepcopy(original.__dict__)
    del original_dict["params"]

    results = []
    for key, value in original_dict.items():
        # HACK: Combating python adding the class from which it inherited
        # in front of private attributes
        key = key.replace("_oimComponent_", "")

        value_res = getattr(restored, key)
        if isinstance(value, oimParam):
            res = assert_param_equal(value, value_res)
        else:
            res = value == value_res

        results.append(res)

    for name, param in original.params.items():
        results.append(
            assert_param_equal(param, restored.params.get(name, oimParam()))
        )

    result = all(results)
    assert result
    return result


def assert_model_equal(original: oimModel, restored: oimModel) -> bool:
    """Compares an original to a restored oimModel."""
    results = []
    for comp_org, comp_res in zip(original.components, restored.components):
        results.append(assert_component_equal(comp_org, comp_res))

    for param_org, param_res in zip(original.extParams, restored.extParams):
        results.append(param_org == param_res)

    return all(results)
