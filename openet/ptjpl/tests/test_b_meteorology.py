import ee
import pytest

import openet.ptjpl.meteorology as meteorology
import openet.ptjpl.utils as utils


@pytest.mark.parametrize(
    'Ta_C, expected',
    [
        [0.0, 0.611],
        [10.0, 1.226980],
        [20.0, 2.334178],
        [30.0, 4.232179],
    ]
)
def test_SVP_from_Ta(Ta_C, expected, tol=0.000001):
    output = utils.getinfo(meteorology.SVP_from_Ta(ee.Number(Ta_C)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'Ta_C, expected',
    [
        [0.0, 0.044450],
        [10.0, 0.082190],
        [20.0, 0.144439],
        [30.0, 0.242659],
    ]
)
def test_delta_from_Ta(Ta_C, expected, tol=0.000001):
    output = utils.getinfo(meteorology.delta_from_Ta(ee.Number(Ta_C)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'T_K, expected',
    [
        [300, 26.85],
        [273.15, 0.0],
    ]
)
def test_kelvin_to_celsius(T_K, expected, tol=0.000001):
    output = utils.getinfo(meteorology.kelvin_to_celsius(ee.Number(T_K)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'P_Pa, expected',
    [
        [10000.0, 10.0],
    ]
)
def test_pascal_to_kilopascal(P_Pa, expected, tol=0.000001):
    output = utils.getinfo(meteorology.pascal_to_kilopascal(ee.Number(P_Pa)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'Ta_C, Ea_kPa, expected',
    [
        [0, 1.22, 1.0],
        [6.5, 1.22, 1.0],  # Test .max(1)
        [10, 1.22, 1.226980],
        [20, 1.22, 2.334178],
    ]
)
def test_meteorology_svp(Ta_C, Ea_kPa, expected, tol=0.000001):
    # CGM - This test uses images due to .updateMask() call in meteorology
    SVP, VPD, RH = meteorology.meteorology(ee.Image.constant(Ta_C), ee.Image.constant(Ea_kPa))
    output = utils.constant_image_value(SVP.rename(['svp']))
    assert abs(output['svp'] - expected) <= tol


@pytest.mark.parametrize(
    'Ta_C, Ea_kPa, expected',
    [
        [20, 1.22, 2.334178 - 1.22],
        [20, 3.0, None],  # Trigger negative VPD masking
    ]
)
def test_meteorology_vpd(Ta_C, Ea_kPa, expected, tol=0.000001):
    # CGM - This test uses images due to .updateMask() call in meteorology
    SVP, VPD, RH = meteorology.meteorology(ee.Image.constant(Ta_C), ee.Image.constant(Ea_kPa))
    output = utils.constant_image_value(VPD.rename(['vpd']))
    if expected is None and output['vpd'] is None:
        assert True
    else:
        assert abs(output['vpd'] - expected) <= tol


@pytest.mark.parametrize(
    'Ta_C, Ea_kPa, expected',
    [
        [20, 1.22, 0.522668],
        [20, 2.34, 1.0],  # Test .min(1)
    ]
)
def test_meteorology_rh(Ta_C, Ea_kPa, expected, tol=0.000001):
    # CGM - This test uses images due to .updateMask() call in meteorology
    SVP, VPD, RH = meteorology.meteorology(ee.Image.constant(Ta_C), ee.Image.constant(Ea_kPa))
    output = utils.constant_image_value(RH.rename(['rh']))
    assert abs(output['rh'] - expected) <= tol


@pytest.mark.parametrize(
    'RH, expected',
    [
        [0.5, 0.5 ** 4],
    ]
)
def test_fwet_from_RH(RH, expected, tol=0.000001):
    output = utils.getinfo(meteorology.fwet_from_RH(ee.Number(RH)))
    assert abs(output - expected) <= tol
