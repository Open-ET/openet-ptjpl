"""
This module calculates instantaneous to daily time integration for PT-JPL.
Gregory Halverson, Jet Propulsion Laboratory
"""
import math

__author__ = "Gregory Halverson"


def daily_integration(instantaneous_net_radiation, hour_of_day, sunrise_hour, daylight_hours):
    # calculate daily net radiation using solar parameters
    # this is the average rate of energy transfer from sunrise to sunset
    # in watts per square meter
    # watts are joules per second
    # to get the total amount of energy transferred, factor seconds out of joules
    # the number of seconds for which this average is representative is (daylight_hours * 3600)
    # documented in verma et al, bisht et al, and lagouarde et al

    # CGM - Swapped order of subtraction so that hour_of_day can be a number
    daily_net_radiation = (
        instantaneous_net_radiation
        .multiply(1.6 / math.pi)
        .divide(sunrise_hour.multiply(-1).add(hour_of_day).multiply(math.pi)
                .divide(daylight_hours).sin())
    )
    # daily_net_radiation = (
    #     1.6 * instantaneous_net_radiation /
    #     (math.pi * np.sin(math.pi * (hour_of_day - sunrise_hour) / (daylight_hours))
    # )

    return daily_net_radiation


def calculate_vapor(LE_daily, daylight_hours):
    # convert length of day in hours to seconds
    daylight_seconds = daylight_hours.multiply(3600.0)

    # constant latent heat of vaporization for water: the number of joules of
    # energy it takes to evaporate one kilogram
    LATENT_VAPORIZATION_JOULES_PER_KILOGRAM = 2450000.0

    # factor seconds out of watts to get joules and divide by latent heat of
    # vaporization to get kilograms
    ET_daily_kg = (
        LE_daily.multiply(daylight_seconds)
        .divide(LATENT_VAPORIZATION_JOULES_PER_KILOGRAM)
        .max(0)
    )
    # ET_daily_kg = np.clip(
    #     LE_daily * daylight_seconds / LATENT_VAPORIZATION_JOULES_PER_KILOGRAM,
    #     0, None
    # )

    return ET_daily_kg
