import datetime

import ee
import pytest

import openet.ptjpl.landsat as landsat
import openet.ptjpl.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
# SCENE_ID = 'LC08_042035_20150713'
SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_TIME = utils.millis(SCENE_DT)


def sr_image(blue=0.2, green=0.2, red=0.2, nir=0.7, swir1=0.2, swir2=0.2, bt=300):
    """Construct a fake Landsat 8 image with renamed bands"""
    return (
        ee.Image.constant([blue, green, red, nir, swir1, swir2, bt])
        .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir'])
        .set({
            'system:time_start': ee.Date(SCENE_DATE).millis(),
            'k1_constant': ee.Number(607.76),
            'k2_constant': ee.Number(1260.56),
        })
    )


# Test the static methods of the class first
# Do these need to be inside the TestClass?
def test_Image_ndvi_band_name():
    output = landsat.ndvi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'ndvi'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 9.0 / 55, -0.1],
        [0.2, 0.2, 0.0],
        [0.1, 11.0 / 90,  0.1],
        [0.2, 0.3, 0.2],
        [0.1, 13.0 / 70, 0.3],
        [0.3, 0.7, 0.4],
        [0.2, 0.6, 0.5],
        [0.2, 0.8, 0.6],
        [0.1, 17.0 / 30, 0.7],
        [0.2, 0.7, 0.55555555],
    ]
)
def test_Image_ndvi_calculation(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.ndvi(sr_image(red=red, nir=nir)))
    assert abs(output['ndvi'] - expected) <= tol


def test_Image_ndwi_band_name():
    output = landsat.ndwi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'ndwi'


def test_Image_ndwi_calculation(green=0.2, nir=0.2, expected=0, tol=0.000001):
    output = utils.constant_image_value(landsat.ndwi(sr_image(green=green, nir=nir)))
    assert abs(output['ndwi'] - expected) <= tol


def test_Image_mndwi_band_name():
    output = landsat.mndwi(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'mndwi'


@pytest.mark.parametrize(
    'green, swir2, expected',
    [
        [0.2, 0.2, 0.0],
        [0.3, 0.2, 0.2],
    ]
)
def test_Image_mndwi_calculation(green, swir2, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.mndwi(sr_image(green=green, swir2=swir2)))
    assert abs(output['mndwi'] - expected) <= tol


def test_Image_wri_band_name():
    output = landsat.wri(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'wri'


@pytest.mark.parametrize(
    'green, red, nir, swir2, expected',
    [
        [0.2, 0.2, 0.7, 0.2, 0.444444],
        [0.2, 0.2, 0.2, 0.2, 1.0],
    ]
)
def test_Image_wri_calculation(green, red, nir, swir2, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.wri(
        sr_image(green=green, red=red, nir=nir, swir2=swir2)
    ))
    assert abs(output['wri'] - expected) <= tol


def test_Image_emissivity_ptjpl_band_name():
    output = landsat.emissivity_ptjpl(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'emissivity'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 9.0 / 55, 0.985],      # -0.1
        [0.2, 0.2,  0.977],          # 0.0
        [0.1, 11.0 / 90,  0.977],    # 0.1
        [0.2, 0.2999, 0.977],        # 0.3- (0.3 NIR isn't exactly an NDVI of 0.2)
        [0.2, 0.3001, 0.986335],     # 0.3+
        [0.1, 13.0 / 70, 0.986742],  # 0.3
        [0.3, 0.7, 0.987964],        # 0.4
        [0.2, 0.6, 0.99],            # 0.5
        [0.2, 0.8, 0.99],            # 0.6
        [0.1, 17.0 / 30, 0.99],      # 0.7
        [0.1, 0.9, 0.99],            # 0.8
    ]
)
def test_Image_emissivity_ptjpl_calculation(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.emissivity_ptjpl(sr_image(red=red, nir=nir)))
    assert abs(output['emissivity'] - expected) <= tol


def test_Image_emissivity_metric_band_name():
    output = landsat.emissivity_metric(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'emissivity'


@pytest.mark.parametrize(
    'red, nir, expected',
    [
        [0.2, 9.0 / 55, 0.99],         # -0.1
        [0.2, 0.2,  0.99],             # 0.0
        [0.1, 11.0 / 90,  0.9700233],  # 0.1
        [0.2, 0.3, 0.970187],          # 0.2 (ish)
        [0.1, 13.0 / 70, 0.970630],    # 0.3
        [0.3, 0.7, 0.971493],          # 0.4
        [0.2, 0.6, 0.972917],          # 0.5
        [0.2, 0.8, 0.975040],          # 0.6
        [0.1, 17.0 / 30, 0.978003],    # 0.7
        [0.1, 0.9, 0.98],              # 0.8
    ]
)
def test_Image_emissivity_metric_calculation(red, nir, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.emissivity_metric(sr_image(red=red, nir=nir)))
    assert abs(output['emissivity'] - expected) <= tol


def test_Image_lst_band_name():
    output = landsat.lst(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'lst'


@pytest.mark.parametrize(
    'red, nir, bt, expected',
    [
        [0.2, 0.7, 300, 304.487567],  # METRIC Emissivity
        # [0.2, 0.7, 300, 303.471031],  # SSEBop Emissivity
    ]
)
def test_Image_lst_calculation(red, nir, bt, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.lst(sr_image(red=red, nir=nir, bt=bt)))
    assert abs(output['lst'] - expected) <= tol


def test_Image_albedo_band_name():
    output = landsat.albedo_metric(sr_image()).getInfo()['bands'][0]['id']
    assert output == 'albedo'


@pytest.mark.parametrize(
    'blue, green, red, nir, swir1, swir2, expected',
    [
        [0.2, 0.2, 0.2, 0.7, 0.2, 0.2, 0.3555],
    ]
)
def test_Image_albedo_calculation(blue, green, red, nir, swir1, swir2, expected, tol=0.000001):
    output = utils.constant_image_value(landsat.albedo_metric(
        sr_image(blue=blue, green=green, red=red, nir=nir, swir1=swir1, swir2=swir2)
    ))
    assert abs(output['albedo'] - expected) <= tol
