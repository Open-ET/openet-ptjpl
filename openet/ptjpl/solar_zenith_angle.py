"""
This module calculates solar zenith angle.
Gregory Halverson, Jet Propulsion Laboratory
"""
import math

__author__ = 'Gregory Halverson'


def sza_deg_from_lat_dec_hour(latitude, solar_dec_deg, hour):
    """
    This function calculates solar zenith angle from longitude, solar declination, and solar time.
    SZA calculated by this function matching SZA provided by MOD07 to within 0.4 degrees.

    :param latitude: latitude in degrees
    :param solar_dec_deg: solar declination in degrees
    :param hour: solar time in hours
    :return: solar zenith angle in degrees
    """
    # convert angles to radians
    latitude_rad = latitude.multiply(math.pi / 180)
    solar_dec_rad = solar_dec_deg.multiply(math.pi / 180)
    hour_angle_deg = hour.multiply(15.0).subtract(180.0)
    # hour_angle_deg = hour * 15.0 - 180.0
    hour_angle_rad = hour_angle_deg.multiply(math.pi / 180)
    sza_rad = (
        latitude_rad.sin().multiply(solar_dec_rad.sin())
        .add(latitude_rad.cos().multiply(solar_dec_rad.cos()).multiply(hour_angle_rad.cos()))
        .acos()
    )
    # sza_rad = arccos(
    #     sin(latitude_rad) * sin(solar_dec_rad) +
    #     cos(latitude_rad) * cos(solar_dec_rad) * cos(hour_angle_rad)
    # )
    sza_deg = sza_rad.multiply(180 / math.pi)

    return sza_deg
