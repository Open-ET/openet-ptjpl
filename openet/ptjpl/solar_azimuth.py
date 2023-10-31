"""
This module calculates solar azimuth.
Gregory Halverson, Jet Propulsion Laboratory
"""
import math

__author__ = 'Gregory Halverson'


def calculate_solar_azimuth(solar_dec_deg, sza_deg, hour):
    """
    This function calculates solar azimuth.
    :param latitude: latitude in degrees
    :param solar_dec_deg: solar declination in degrees
    :param sza_deg: solar zenith angle in degrees
    :return: solar azimuth in degrees
    """
    # convert angles to radians
    solar_dec_rad = solar_dec_deg.multiply(math.pi / 180)
    sza_rad = sza_deg.multiply(math.pi / 180)
    hour_angle_deg = hour.multiply(15.0).subtract(180.0)
    # hour_angle_deg = hour * 15.0 - 180.0
    hour_angle_rad = hour_angle_deg.multiply(math.pi / 180)

    # https://en.wikipedia.org/wiki/Solar_azimuth_angle
    # solar_azimuth_rad = arccos(
    #     (sin(solar_dec_rad) - cos(sza_rad) * sin(latitude_rad)) /
    #     (sin(sza_rad) * cos(latitude_rad)))
    solar_azimuth_rad = (
        hour_angle_rad.sin().multiply(solar_dec_rad.cos())
        .multiply(-1).divide(sza_rad.sin()).asin()
    )
    solar_azimuth_deg = solar_azimuth_rad.multiply(180 / math.pi)

    return solar_azimuth_deg
