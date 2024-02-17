import ee
import pytest

import openet.ptjpl.vegetation as vegetation
import openet.ptjpl.utils as utils


def test_savi_from_ndvi(ndvi=0.5, expected=0.5 * 0.45 + 0.132, tol=0.000001):
    output = utils.constant_image_value(vegetation.savi_from_ndvi(ee.Image.constant(ndvi)))
    assert abs(output['constant'] - expected) <= tol


@pytest.mark.parametrize(
    'savi, expected',
    [
        [-0.1, 0],
        [0.5, 0.5 * 1.3632 - 0.048],
        [1.0, 1],
    ]
)
def test_fAPAR_from_savi(savi, expected, tol=0.000001):
    output = utils.constant_image_value(vegetation.fAPAR_from_savi(ee.Image.constant(savi)))
    assert abs(output['constant'] - expected) <= tol


@pytest.mark.parametrize(
    'ndvi, expected', [[0.5, (0.45 * 0.5 + 0.132) * 1.3632 - 0.048]])
def test_fAPAR_from_ndvi(ndvi, expected, tol=0.000001):
    output = utils.constant_image_value(vegetation.fAPAR_from_ndvi(ee.Image.constant(ndvi)))
    assert abs(output['constant'] - expected) <= tol


@pytest.mark.parametrize(
    'ndvi, expected',
    [
        [-0.1, 0],
        [0.5, 0.45],
        [1.1, 1],
    ]
)
def test_fIPAR_from_ndvi(ndvi, expected, tol=0.000001):
    output = utils.constant_image_value(vegetation.fIPAR_from_ndvi(ee.Image.constant(ndvi)))
    assert abs(output['constant'] - expected) <= tol
