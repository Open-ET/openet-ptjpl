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


def normalized_difference(image, b1_name, b2_name):
    """Generic normalized difference index function to handle negative reflectance values"""

    # Force the input values to be at greater than or equal to zero
    #   since C02 surface reflectance values can be negative
    #   but the normalizedDifference function will return nodata
    ndi = (
        image.select([b1_name, b2_name]).max(0)
        .normalizedDifference([b1_name, b2_name])
    )

    # Assume that low reflectance values are unreliable for computing the index
    # If both reflectance values are below the minimum, set the output to 0
    # If either of the reflectance values was below 0, set the output to 0
    b1 = image.select([b1_name])
    b2 = image.select([b2_name])
    ndi = ndi.where(b1.lt(0.01).And(b2.lt(0.01)), 0)
    # ndi = ndi.where(b1.lt(0).Or(b2.lt(0)), 0)

    return ndi.rename('ndi')


def ndvi(landsat_image):
    """Normalized Difference Vegetation Index (NDVI)"""
    return normalized_difference(landsat_image, 'nir', 'red').rename(['ndvi'])
    # return landsat_image.normalizedDifference(['nir', 'red']).rename(['ndvi'])


def ndwi(landsat_image):
    """Normalized Difference Water Index (NDWI)"""
    return normalized_difference(landsat_image, 'green', 'nir').rename(['ndwi'])
    # return landsat_image.normalizedDifference(['green', 'nir']).rename(['ndwi'])


def mndwi(landsat_image):
    """Modified Normalized Difference Water Index (MNDWI)"""
    return normalized_difference(landsat_image, 'green', 'swir2').rename(['mndwi'])
    # return landsat_image.normalizedDifference(['green', 'swir2']).rename(['mndwi'])


def wri(landsat_image):
    """Water Ratio Index (WRI)"""
    # Force minimum value to be > 0 to avoid divide by zero
    output = landsat_image.max(0.001).expression('(b("green") + b("red")) / (b("nir") + b("swir2"))')
    # # TODO: Should the WRI value be overwritten to indicate water (>1)
    #    if all reflectance values are low?
    # output = output.where(
    #     landsat_image.select(['nir']).lte(0).And(landsat_image.select(['swir2']).lte(0)), 2
    # )
    return output.rename(['wri'])


def water_mask(landsat_image, gsw_extent_flag=True):
    """Water pixel identification

    Parameters
    ----------
    landsat_image : ee.Image
        "Prepped" Landsat image with standardized band names.
    gsw_extent_flag : boolean
        If True, apply the global surface water extent mask to the QA_PIXEL water mask
        The default is True.

    Returns
    -------
    ee.Image

    """
    # Start with a combination of indices to identify water
    # Adding an albedo threshold limit to help catch saturated pixels that have
    #   water index values (but it may be better to use the QA_RADSAT band instead)
    # The conditionals should probably all be changed to gte or lte
    #   since normalized_difference function is setting them to zero for low reflectance,
    #   but leaving them as is for now to match original implementation,
    #   and the QA_PIXEL water mask below should catch most water pixels
    ptjl_water_mask = (
        ndwi(landsat_image).gt(0)
        .And(mndwi(landsat_image).gt(0))
        .And(wri(landsat_image).gt(1))
        .And(ndvi(landsat_image).lt(0))
        .And(albedo_metric(landsat_image).lt(0.3))
    )

    # Using the water mask in the QA_PIXEL band to include additional water pixels
    #   since the index values are not always reliable for identifying water when the
    #   surface reflectance values are low
    qa_water_mask = landsat_image.select(['QA_PIXEL']).rightShift(7).bitwiseAnd(1)

    # Including the dynamic surface water max extent layer helps to exclude shadows
    #   that are misclassified as water
    if gsw_extent_flag:
        gsw_mask = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select(['max_extent']).gte(1)
        qa_water_mask = qa_water_mask.And(gsw_mask)

    return ptjl_water_mask.Or(qa_water_mask).rename(['water_mask'])

    # TODO: Test out assigning water when 2 of the 3 masks are True
    #   Not sure if this is necessary but it would be nice to be able to catch
    #   very obvious water that is in a new location not in the GSW layer
    # return ptjl_water_mask.add(qa_water_mask).add(gsw_mask).gte(2).rename(['water_mask'])
