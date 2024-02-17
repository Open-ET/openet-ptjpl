import ee
import pytest

import openet.ptjpl.daylight_hours as dh
import openet.ptjpl.utils as utils


@pytest.mark.parametrize(
    'doy, expected',
    [
        [197, 3.373984],
    ]
)
def test_day_angle_rad_from_doy(doy, expected, tol=0.000001):
    output = utils.getinfo(dh.day_angle_rad_from_doy(ee.Number(doy)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'day_angle_rad, expected',
    [
        [3.373984, 21.507808],
    ]
)
def test_solar_dec_deg_from_day_angle_rad(day_angle_rad, expected, tol=0.000001):
    output = utils.getinfo(dh.solar_dec_deg_from_day_angle_rad(ee.Number(day_angle_rad)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'doy, latitude, expected',
    [
        [197, 39.0, 108.6091],
    ]
)
def test_sha_deg_from_doy_lat(doy, latitude, expected, tol=0.000001):
    output = utils.constant_image_value(dh.sha_deg_from_doy_lat(
        ee.Number(doy), ee.Image.constant(latitude)
    ))
    # output = utils.getinfo(dh.sha_deg_from_doy_lat(ee.Number(doy), ee.Number(latitude)))
    assert abs(output['constant'] - expected) <= tol


@pytest.mark.parametrize(
    'sha_deg, expected',
    [
        [108.6091, 4.759393],
    ]
)
def test_sunrise_from_sha(sha_deg, expected, tol=0.000001):
    output = utils.getinfo(dh.sunrise_from_sha(ee.Number(sha_deg)))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'sha_deg, expected',
    [
        [108.6091, 14.481213],
    ]
)
def test_daylight_from_sha(sha_deg, expected, tol=0.000001):
    output = utils.getinfo(dh.daylight_from_sha(ee.Number(sha_deg)))
    assert abs(output - expected) <= tol
