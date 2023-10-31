import ee
import pytest

import openet.ptjpl.solar_azimuth as solar_azimuth
import openet.ptjpl.utils as utils


@pytest.mark.parametrize(
    'solar_dec_deg, sza_deg, hour, expected',
    [
        [0, 0, 0, 0],
    ]
)
def test_calculate_solar_azimuth(solar_dec_deg, sza_deg, hour, expected, tol=0.000001):
    output = utils.getinfo(solar_azimuth.calculate_solar_azimuth(
        ee.Number(solar_dec_deg), ee.Number(sza_deg), ee.Number(hour)
    ))
    assert abs(output - expected) <= tol
