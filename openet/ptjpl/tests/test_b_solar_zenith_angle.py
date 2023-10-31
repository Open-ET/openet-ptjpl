import ee
import pytest

import openet.ptjpl.solar_zenith_angle as sza
import openet.ptjpl.utils as utils


@pytest.mark.parametrize(
    'latitude, solar_dec_deg, hour, expected',
    [
        [0, 0, 0, 180],
    ]
)
def test_sza(latitude, solar_dec_deg, hour, expected, tol=0.000001):
    output = utils.getinfo(sza.sza_deg_from_lat_dec_hour(
        ee.Number(latitude), ee.Number(solar_dec_deg), ee.Number(hour)
    ))
    assert abs(output - expected) <= tol
