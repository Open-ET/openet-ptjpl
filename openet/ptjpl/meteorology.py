"""
This module calculates near-surface meteorology for PT-JPL.
Gregory Halverson, Jet Propulsion Laboratory
"""

__author__ = "Gregory Halverson"


def SVP_from_Ta(Ta_C):
    """
    Calculates saturation vapor pressure for air temperature.
    :param Ta_C: air temperature in celsius
    :return: saturation vapor pressure in kilopascals
    """
    return Ta_C.multiply(17.27).divide(Ta_C.add(237.7)).exp().multiply(0.611)
    # return Ta_C.expression('0.611 * exp((Ta_C * 17.27) / (Ta_C + 237.7))', {'Ta_C': Ta_C})
    # return 0.611 * np.exp((Ta_C * 17.27) / (Ta_C + 237.7))


def delta_from_Ta(Ta_C):
    """
    Calculates slope of saturation vapor pressure to air temperature.
    :param Ta_C: air temperature in celsius
    :return: slope of the vapor pressure curve in kilopascals/kelvin
    """
    return (
        Ta_C.multiply(17.27).divide(Ta_C.add(237.7)).exp()
        .multiply(0.6108 * 4098).divide(Ta_C.add(237.3).pow(2.0))
    )
    # return Ta_C.expression(
    #     '4098 * 0.6108 * exp(17.27 * Ta_C / (237.7 + Ta_C)) / (Ta_C + 237.3) ** 2',
    #     {'Ta_C': Ta_C}
    # )
    # return 4098 * (0.6108 * np.exp(17.27 * Ta_C / (237.7 + Ta_C))) / (Ta_C + 237.3) ** 2


def kelvin_to_celsius(temperature_K):
    """
    Reduces temperature in kelvin to celsius.
    :param temperature_K: temperature in kelvin
    :return: temperature in celsius
    """
    return temperature_K.subtract(273.15)


def pascal_to_kilopascal(pressure_Pa):
    """
    Scales order of magnitude of pressure in pascals to kilopascals.
    :param pressure_Pa: pressure in pascals
    :return: pressure in kilopascals
    """
    return pressure_Pa.divide(1000.0)


# CGM - It might be cleaner to split this into three separate functions
def meteorology(Ta_C, Ea_kPa):
    """
    Calculates saturation vapor pressure, vapor pressure deficit, and relative
    humidity from air temperature and water vapor pressure.
    :param Ta_C: air temperature in celsius
    :type Ta_C: ee.Image
    :param Ea_kPa: water vapor pressure in kilopascals
    :return:
        3-tuple of:
            saturation vapor pressure in kilopascals
            vapor pressure deficit in kilopascals
            relative humidity as a proportion between 0 and 1
    """
    # calculate saturation vapor pressure in kPa from air temperature in celsius
    SVP = SVP_from_Ta(Ta_C)

    # floor saturation vapor pressure at 1
    SVP = SVP.max(1)

    # calculate vapor pressure deficit from water vapor pressure
    VPD = SVP.subtract(Ea_kPa)

    # CGM - I think the .gte(0) is correct since you are telling updateMask()
    #   which pixels are valid
    # TODO: Check if the updateMask() can can be skipped (and then change tests from images to numbers)
    # lower bound of vapor pressure deficit is zero, negative values replaced with nodata
    VPD = VPD.updateMask(VPD.gte(0))
    # VPD = np.where(VPD < 0, np.nan, VPD)

    # calculate relative humidity from water vapor pressure and saturation vapor pressure
    RH = Ea_kPa.divide(SVP)

    # upper bound of relative humidity is one, results higher than one are capped at one
    RH = RH.min(1)

    return SVP, VPD, RH


def fwet_from_RH(RH):
    """
    Calculates relative surface wetness from relative humidity.
    :param RH: relative humidity as a proportion between 0 and 1
    :return: relative surface wetness as a proportion between 0 and 1
    """
    return RH.pow(4)
