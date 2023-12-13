import ee


# TODO: Change this to a Landsat class similar to DisALEXI
def albedo_disalexi(landsat_image):
    """Total shortwave broadband albedo following [Liang2001]

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.

    Returns
    -------
    albedo : ee.Image

    Notes
    -----
    The Python DisALEXI code had the following line and comment:
        "bands = [1, 3, 4, 5, 7]  # dont use blue"
    IDL code and [Liang2001] indicate that the green band is not used.
    Coefficients were derived for Landsat 7 ETM+, but were found to be
        "suitable" to Landsat 4/5 TM also.

    References
    ----------
    .. [Liang2001] Shunlin Liang (2001).  Narrowband to broadband conversions
        of land surface albedo - I Algorithms,
        Remote Sensing of Environment, Volume 76, Issue2, Pages 213-238,
        http://doi.org/10.1016/S0034-4257(00)00205-4

    """
    albedo_img = (
        landsat_image
        .select(['blue', 'red', 'nir', 'swir1', 'swir2'])
        .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
    )
    return (
        albedo_img.select([0])
        .add(albedo_img.select([1])).add(albedo_img.select([2]))
        .add(albedo_img.select([3])).add(albedo_img.select([4]))
        .subtract(0.0018)
        .rename(['albedo'])
    )

    # # CGM - Using a sum reducer was returning an unbounded image
    # return (
    #     ee.Image(self.input_image)
    #     .select(['blue', 'red', 'nir', 'swir1', 'swir2'])
    #     .multiply([0.356, 0.130, 0.373, 0.085, 0.072])
    #     .reduce(ee.Reducer.sum())
    #     .subtract(0.0018)
    #     .rename(['albedo'])
    # )


def albedo_metric(landsat_image):
    """Total shortwave broadband albedo following [Tasumi2008]

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.

    Returns
    -------
    albedo : ee.Image

    References
    ----------
    .. [Tasumi2008] Tasumi, M., Allen, R., and Trezza, R. (2008).
        At-surface reflectance and albedo from satellite for operational
        calculation of land surface energy balance.
        Journal of Hydrologic Engineering, 13(2), 51-63.
        https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

    """
    albedo_img = (
        landsat_image
        .select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
        .multiply([0.254, 0.149, 0.147, 0.311, 0.103, 0.036])
    )
    return (
        albedo_img.select([0]).add(albedo_img.select([1]))
        .add(albedo_img.select([2])).add(albedo_img.select([3]))
        .add(albedo_img.select([4])).add(albedo_img.select([5]))
        .rename(['albedo'])
    )

    # # CGM - Using a sum reducer was returning an unbounded image
    # return (
    #     ee.Image(landsat_image)
    #     .select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
    #     .multiply([0.254, 0.149, 0.147, 0.311, 0.103, 0.036])
    #     .reduce(ee.Reducer.sum())
    #     .rename(['albedo'])
    # )


def emissivity_metric(landsat_image):
    """Emissivity as a function of NDVI

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    # TODO: Add METRIC emissivity reference(s)
    # TODO: Add METRIC LAI reference

    """
    ndvi_img = ndvi(landsat_image)
    lai_img = ndvi_img.pow(3).multiply(7.0).clamp(0, 6).rename(['lai'])

    # Initial values are for NDVI > 0 and LAI <= 3
    return (
        lai_img.divide(300).add(0.97)
        .where(ndvi_img.lte(0), 0.99)
        .where(ndvi_img.gt(0).And(lai_img.gt(3)), 0.98)
        .rename(['emissivity'])
    )


def emissivity_ptjpl(landsat_image):
    """Emissivity as a function of NDVI

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.

    Returns
    -------
    ee.Image

    References
    ----------
    .. [Sobrino2004] Sobrino, J., J. Jiménez-Muñoz, & L. Paolini (2004).
        Land surface temperature retrieval from LANDSAT TM 5.
        Remote Sensing of Environment, 90(4), 434-440.
        https://doi.org/10.1016/j.rse.2004.02.003

    """
    ndvi_img = ndvi(landsat_image)

    # Vegetation proportion
    Pv = ndvi_img.expression('((ndvi - 0.2) / 0.3) ** 2', {'ndvi': ndvi_img})
    # ndviRangevalue = ndvi_image.where(
    #     ndvi_image.gte(0.2).And(ndvi_image.lte(0.5)), ndvi_image
    # )
    # Pv = ndviRangevalue.expression(
    #     '(((ndviRangevalue - 0.2) / 0.3) ** 2', {'ndviRangevalue': ndviRangevalue}
    # )

    # Assuming typical Soil Emissivity of 0.97 and Veg Emissivity of 0.99
    #   and shape Factor mean value of 0.553
    dE = Pv.expression('(1 - 0.97) * (1 - Pv) * (0.55 * 0.99)', {'Pv': Pv})
    RangeEmiss = dE.expression('(0.99 * Pv) + (0.97 * (1 - Pv)) + dE', {'Pv': Pv, 'dE': dE})

    # RangeEmiss = 0.989
    return (
        ndvi_img
        .where(ndvi_img.lt(0), 0.985)
        .where(ndvi_img.gte(0).And(ndvi_img.lt(0.2)), 0.977)
        .where(ndvi_img.gt(0.5), 0.99)
        .where(ndvi_img.gte(0.2).And(ndvi_img.lte(0.5)), RangeEmiss)
        .clamp(0.977, 0.99)
        .rename(['emissivity'])
    )


def lst(landsat_image):
    """Emissivity corrected land surface temperature (LST) from brightness
    temperature

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.
        Image must have 'k1_constant' and 'k2_constant' properties.

    Returns
    -------
    ee.Image

    Notes
    -----
    The corrected radiation coefficients were derived from a small number
    of scenes in southern Idaho [Allen2007] and may not be appropriate for
    other areas.

    References
    ----------
    .. [Allen2007] R. Allen, M. Tasumi, R. Trezza (2007),
        Satellite-Based Energy Balance for Mapping Evapotranspiration with
        Internalized Calibration (METRIC) Model,
        Journal of Irrigation and Drainage Engineering, Vol 133(4),
        http://dx.doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)

    """
    # Get properties from image
    k1 = ee.Number(ee.Image(landsat_image).get('k1_constant'))
    k2 = ee.Number(ee.Image(landsat_image).get('k2_constant'))

    ts_brightness = ee.Image(landsat_image).select(['tir'])
    emissivity_img = emissivity_metric(landsat_image)

    # First back out radiance from brightness temperature
    # Then recalculate emissivity corrected Ts
    thermal_rad_toa = ts_brightness.expression(
        'k1 / (exp(k2 / ts_brightness) - 1)',
        {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2}
    )

    # tnb = 0.866   # narrow band transmissivity of air
    # rp = 0.91     # path radiance
    # rsky = 1.32   # narrow band clear sky downward thermal radiation
    rc = thermal_rad_toa.expression(
        '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
        {
            'thermal_rad_toa': thermal_rad_toa,
            'emiss': emissivity_img,
            'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32
        }
    )
    lst = rc.expression(
        'k2 / log(emiss * k1 / rc + 1)',
        {'emiss': emissivity_img, 'rc': rc, 'k1': k1, 'k2': k2}
    )

    return lst.rename(['lst'])


def ndvi(landsat_image):
    """Normalized Difference Vegetation Index (NDVI)"""
    return ee.Image(landsat_image).normalizedDifference(['nir', 'red']).rename(['ndvi'])


def ndwi(landsat_image):
    """Normalized Difference Water Index (NDWI)"""
    return ee.Image(landsat_image).normalizedDifference(['green', 'nir']).rename(['ndwi'])


def mndwi(landsat_image):
    """Modified Normalized Difference Water Index (MNDWI)"""
    return ee.Image(landsat_image).normalizedDifference(['green', 'swir2']).rename(['mndwi'])


def wri(landsat_image):
    """Water Ratio Index (WRI)"""
    image = ee.Image(landsat_image)
    WRI = image.expression(
        '(green + red) / (NIR + SWIR)',
        {
            'green': image.select(['green']),
            'red': image.select(['red']),
            'NIR': image.select(['nir']),
            'SWIR': image.select(['swir2'])
        }
    )

    return WRI.rename(['wri'])
