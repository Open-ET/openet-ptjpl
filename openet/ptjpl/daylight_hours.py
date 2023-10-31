"""
This module calculates sunrise hour and daylight hours.
Gregory Halverson, Jet Propulsion Laboratory
"""
import math

__author__ = "Gregory Halverson"


def day_angle_rad_from_doy(doy):
    """
    Calculate day angle in radians from day of year between 1 and 365.
    :param doy: day of year
    :type doy: ee.Number

    """
    return doy.subtract(1).multiply(2 * math.pi / 365)
    # return (2 * pi * (doy - 1)) / 365


def solar_dec_deg_from_day_angle_rad(day_angle_rad):
    """
    This function calculates solar declination in degrees from day angle in radians.
    """
    # CGM - Moved the subtract operation (minus signs) into the multiplication
    return (
        day_angle_rad.cos().multiply(-0.399912)
        .add(day_angle_rad.sin().multiply(0.070257))
        .add(day_angle_rad.multiply(2).cos().multiply(-0.006758))
        .add(day_angle_rad.multiply(2).sin().multiply(0.000907))
        .add(day_angle_rad.multiply(3).cos().multiply(-0.002697))
        .add(day_angle_rad.multiply(3).sin().multiply(0.00148))
        .add(0.006918).multiply(180 / math.pi)
    )
    # return (0.006918
    #         - 0.399912 * cos(day_angle_rad)
    #         + 0.070257 * sin(day_angle_rad)
    #         - 0.006758 * cos(2 * day_angle_rad)
    #         + 0.000907 * sin(2 * day_angle_rad)
    #         - 0.002697 * cos(3 * day_angle_rad)
    #         + 0.00148 * sin(3 * day_angle_rad)) * (180 / math.pi)


def sha_deg_from_doy_lat(doy, latitude):
    """
    This function calculates sunrise hour angle in degrees from latitude in
    degrees and day of year between 1 and 365.
    """
    # calculate day angle in radians
    day_angle_rad = day_angle_rad_from_doy(doy)

    # calculate solar declination in degrees
    solar_dec_deg = solar_dec_deg_from_day_angle_rad(day_angle_rad)

    # convert latitude to radians
    latitude_rad = latitude.multiply(math.pi / 180)
    # latitude_rad = np.radians(latitude)

    # convert solar declination to radians
    solar_dec_rad = solar_dec_deg.multiply(math.pi / 180)
    # solar_dec_rad = np.radians(solar_dec_deg)

    # calculate cosine of sunrise angle at latitude and solar declination
    # need to keep the cosine for polar correction
    sunrise_cos = latitude_rad.tan().multiply(-1).multiply(solar_dec_rad.tan())
    # sunrise_cos = -np.tan(latitude_rad) * np.tan(solar_dec_rad)

    # calculate sunrise angle in radians from cosine
    sunrise_rad = sunrise_cos.acos()

    # convert to degrees
    sunrise_deg = sunrise_rad.multiply(180 / math.pi)

    # apply polar correction
    # TODO: If possible convert to pure math operations (i.e. don't use .where())
    sunrise_deg = sunrise_deg.where(sunrise_cos.gte(1), 0)
    sunrise_deg = sunrise_deg.where(sunrise_cos.lte(-1), 180)

    return sunrise_deg


def sunrise_from_sha(sha_deg):
    """
    This function calculates sunrise hour from sunrise hour angle in degrees.
    """
    return sha_deg.divide(15.0).multiply(-1).add(12)
    # return 12.0 - (sha_deg / 15.0)


def daylight_from_sha(sha_deg):
    """
    This function calculates daylight hours from sunrise hour angle in degrees.
    """
    return sha_deg.multiply(2.0 / 15.0)
    # return (2.0 / 15.0) * sha_deg
