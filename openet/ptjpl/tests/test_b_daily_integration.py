import ee
import pytest

import openet.ptjpl.daily_integration as daily_integration
import openet.ptjpl.utils as utils


@pytest.mark.parametrize(
    'rn, hour_of_day, sunrise_hour, daylight_hours, expected',
    [
        [1000.0, 10.0, 4.759393, 14.481213, 561.3],
        [1000.0, 11.0, 4.759393, 14.481213, 521.5],
        [1000.0, 12.0, 4.759393, 14.481213, 509.3],
    ]
)
def test_day_angle_rad_from_doy_number(rn, hour_of_day, sunrise_hour,
                                       daylight_hours, expected, tol=0.1):
    # Test the math for all number inputs
    output = utils.getinfo(daily_integration.daily_integration(
        instantaneous_net_radiation=ee.Number(rn),
        hour_of_day=ee.Number(hour_of_day),
        sunrise_hour=ee.Number(sunrise_hour),
        daylight_hours=ee.Number(daylight_hours),
    ))
    assert abs(output - expected) <= tol


@pytest.mark.parametrize(
    'rn, hour_of_day, sunrise_hour, daylight_hours, expected',
    [
        [1000.0, 11.0, 4.759393, 14.481213, 521.5],
    ]
)
def test_day_angle_rad_from_doy_image(rn, hour_of_day, sunrise_hour,
                                      daylight_hours, expected, tol=0.1):
    # Test that hour can be a number when all other inputs are images
    output = utils.constant_image_value(daily_integration.daily_integration(
        instantaneous_net_radiation=ee.Image.constant(rn),
        hour_of_day=ee.Number(hour_of_day),
        sunrise_hour=ee.Image.constant(sunrise_hour),
        daylight_hours=ee.Image.constant(daylight_hours),
    ))
    assert abs(output['constant'] - expected) <= tol


@pytest.mark.parametrize(
    'LE_daily, daylight_hours, expected',
    [
        [500.0, 12.0, 8.8],
        [-10.0, 12.0, 0.0],  # Test the .max(0)
    ]
)
def test_solar_dec_deg_from_day_angle_rad(LE_daily, daylight_hours, expected, tol=0.1):
    output = utils.getinfo(daily_integration.calculate_vapor(
        ee.Number(LE_daily), ee.Number(daylight_hours)
    ))
    assert abs(output - expected) <= tol
