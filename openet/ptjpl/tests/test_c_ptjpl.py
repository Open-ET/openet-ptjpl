import datetime
# import logging

import ee
import pytest

import openet.ptjpl.ptjpl as ptjpl
import openet.ptjpl.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
SCENE_DT = datetime.datetime.utcfromtimestamp(SCENE_TIME / 1000.0)
# SCENE_DT = datetime.datetime.strptime(SCENE_ID[-8:], '%Y%m%d')
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
# SCENE_TIME = utils.millis(SCENE_DT)


@pytest.mark.parametrize(
    'Rn, LAI, water_mask, expected',
    [
        [10, 0, 0, 10],
        [10, 1, 0, 5.488],
        # [10, -1, 0, 0],  # Do we need a test for negative LAI values?
        [10, 3, 1, None],
    ]
)
def test_Model_Rns_calculation(Rn, LAI, water_mask, expected, tol=0.001):
    output = utils.constant_image_value(ptjpl.Rns(
        Rn=ee.Image.constant(Rn), LAI=ee.Image.constant(LAI),
        water_mask=ee.Image.constant(water_mask)
    ))
    if expected is None:
        assert output['constant'] is None
    else:
        assert abs(output['constant'] - expected) <= tol


# @pytest.mark.parametrize(
#     'Rn, fIPAR, Rns, W, water_mask, expected',
#     [
#         [-1, 0.25, 2, 0, 0, 0],  # G is clamped > 0
#         [10, 0.25, 2, 0, 0, 0.35 * 2],  # G is clamped to Rns * GMAX
#         [10, 0.25, 20, 0, 0, 2.4875],  # Set Rns large to test G calculation
#         [10, 0.25, 2, 5, 1, 5],  # Water pixels use W directly
#     ]
# )
# def test_Model_G_calculation(Rn, fIPAR, Rns, W, water_mask, expected, tol=0.0001):
#     output = utils.constant_image_value(ptjpl.G(
#         Rn=ee.Image.constant(Rn), fIPAR=ee.Image.constant(fIPAR),
#         Rns=ee.Image.constant(Rns), W=ee.Image.constant(W),
#         water_mask=ee.Image.constant(water_mask)))
#     assert abs(output['constant'] - expected) <= tol


# @pytest.mark.parametrize(
#     'LST, emissivity, NDVI, albedo, Ta_K, Ea_Pa, SWin, '
#     'Topt, fAPARmax, cloud_mask, datetime, latitude, expected',
#     [
#         [300, 0.99, 0.8, 0.2, 300, 1200, 900,
#          16, 0.6, 0, 1500230731090, 39.05, 11.5328]
#     ]
# )
# def test_Model_et_calculation(LST, emissivity, NDVI, albedo, Ta_K, Ea_Pa, SWin,
#                               Topt, fAPARmax, cloud_mask, datetime, latitude,
#                               expected, tol=0.0001):
#     output_image = ptjpl.et(
#         LST=ee.Image.constant(LST),
#         emissivity=ee.Image.constant(emissivity),
#         NDVI=ee.Image.constant(NDVI),
#         albedo=ee.Image.constant(albedo),
#         Ta_K=ee.Image.constant(Ta_K),
#         Ea_Pa=ee.Image.constant(Ea_Pa),
#         SWin=ee.Image.constant(SWin),
#         Topt=ee.Image.constant(Topt),
#         fAPARmax=ee.Image.constant(fAPARmax),
#         cloud_mask=ee.Image.constant(cloud_mask),
#         datetime=ee.Date(datetime),
#         latitude=ee.Image.constant(latitude),
#     )
#     output = utils.constant_image_value(output_image)
#     logging.debug(f'\n  Target values: {expected}')
#     logging.debug(f'  Output values: {output['constant']}')
#     assert abs(output['constant'] - expected) <= tol
