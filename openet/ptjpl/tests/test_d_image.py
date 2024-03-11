import datetime
# import logging
# import pprint

import ee
import pytest

import openet.ptjpl as ptjpl
import openet.ptjpl.utils as utils
# TODO: import utils from openet.core
# import openet.core.utils as utils


# TODO: Try moving to conftest and/or make a fixture
COLL_ID = 'LANDSAT/LC08/C02/T1_L2/'
SCENE_ID = 'LC08_044033_20170716'
SCENE_TIME = 1500230731090
SCENE_DT = datetime.datetime.utcfromtimestamp(SCENE_TIME / 1000.0)
SCENE_DATE = SCENE_DT.strftime('%Y-%m-%d')
SCENE_DOY = int(SCENE_DT.strftime('%j'))
SCENE_0UTC_DT = datetime.datetime.strptime(SCENE_DATE, '%Y-%m-%d')
# TEST_POINT = (-121.5265, 38.7399)
TEST_POINT = [-120.113, 36.336]


# Should these be test fixtures instead?
# I'm not sure how to make them fixtures and allow input parameters
# CGM - This function is not currently used
# def sr_image(blue=0.2, green=0.2, red=0.2, nir=0.7, swir1=0.2, swir2=0.2, bt=300):
#     """Construct a fake Landsat 8 image with renamed bands"""
#     return ee.Image.constant([blue, green, red, nir, swir1, swir2, bt]) \
#         .rename(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst']) \
#         .set({
#             'system:time_start': SCENE_TIME,
#             'k1_constant': ee.Number(607.76),
#             'k2_constant': ee.Number(1260.56)})


def default_image(albedo=0.2, emissivity=0.99, lst=300, ndvi=0.8, ndwi=0.0, mndwi=0.0, wri=0.45):
    # First construct a fake 'prepped' input image
    return (
        ee.Image.constant([albedo, emissivity, lst, ndvi, ndwi, mndwi, wri])
        .rename(['albedo', 'emissivity', 'lst', 'ndvi', 'ndwi', 'mndwi', 'wri'])
        .set({
            'system:index': SCENE_ID,
            'system:time_start': SCENE_TIME,
            'system:id': COLL_ID + SCENE_ID,
        })
    )


# Setting etr_source and etr_band on the default image to simplify testing
#   but these do not have defaults in the Image class init
def default_image_args(
        albedo=0.2,
        emissivity=0.99,
        lst=300,
        ndvi=0.8,
        ndwi=0.0,
        mndwi=0.0,
        wri=0.45,
        ea_source=1000,
        LWin_source=440,
        rs_source=900,
        ta_source=315,
        windspeed_source=0,
        topt_source=42,
        faparmax_source=0.6,
        latitude=36,
        longitude=-120,
        floor_Topt=True,
        et_reference_source=15,
        et_reference_band='etr',
        et_reference_factor=0.85,
        et_reference_resample='nearest',
        crop_pm_adjust_flag=False,
        crop_pm_adjust_source=1,
        crop_pm_adjust_band=None,
        crop_type_source='USDA/NASS/CDL',
        crop_type_remap='CDL',
        ):
    return {
        'image': default_image(albedo=albedo, emissivity=emissivity, lst=lst,
                               ndvi=ndvi, ndwi=ndwi, mndwi=mndwi, wri=wri),
        'ea_source': ea_source,
        'LWin_source': LWin_source,
        'rs_source': rs_source,
        'ta_source': ta_source,
        'windspeed_source': windspeed_source,
        'topt_source': topt_source,
        'faparmax_source': faparmax_source,
        'latitude': latitude,
        'longitude': longitude,
        'floor_Topt': floor_Topt,
        'et_reference_source': et_reference_source,
        'et_reference_band': et_reference_band,
        'et_reference_factor': et_reference_factor,
        'et_reference_resample': et_reference_resample,
        'crop_pm_adjust_flag': crop_pm_adjust_flag,
        'crop_pm_adjust_source': crop_pm_adjust_source,
        'crop_pm_adjust_band': crop_pm_adjust_band,
        'crop_type_source': crop_type_source,
        'crop_type_remap': crop_type_remap,
    }


def default_image_obj(
        albedo=0.2,
        emissivity=0.99,
        lst=300,
        ndvi=0.8,
        ndwi=0.0,
        mndwi=0.0,
        wri=0.45,
        ea_source=1000,
        LWin_source=440,
        rs_source=900,
        ta_source=315,
        windspeed_source=0,
        topt_source=42,
        faparmax_source=0.6,
        latitude=36,
        longitude=-120,
        floor_Topt=True,
        et_reference_source=15,
        et_reference_band='etr',
        et_reference_factor=0.85,
        et_reference_resample='nearest',
        crop_pm_adjust_flag=False,
        crop_pm_adjust_source=1,
        crop_pm_adjust_band=None,
        crop_type_source='USDA/NASS/CDL',
        crop_type_remap='CDL',
        ):
    return ptjpl.Image(**default_image_args(
        albedo=albedo,
        emissivity=emissivity,
        lst=lst,
        ndvi=ndvi,
        ndwi=ndwi,
        mndwi=mndwi,
        wri=wri,
        ea_source=ea_source,
        LWin_source=LWin_source,
        rs_source=rs_source,
        ta_source=ta_source,
        windspeed_source=windspeed_source,
        topt_source=topt_source,
        faparmax_source=faparmax_source,
        latitude=latitude,
        longitude=longitude,
        floor_Topt=floor_Topt,
        et_reference_source=et_reference_source,
        et_reference_band=et_reference_band,
        et_reference_factor=et_reference_factor,
        et_reference_resample=et_reference_resample,
        crop_pm_adjust_flag=crop_pm_adjust_flag,
        crop_pm_adjust_source=crop_pm_adjust_source,
        crop_pm_adjust_band=crop_pm_adjust_band,
        crop_type_source=crop_type_source,
        crop_type_remap=crop_type_remap,
    ))


def test_Image_init_default_parameters():
    m = ptjpl.Image(default_image())
    assert m.ea_source == 'NLDAS'
    assert m.LWin_source == 'NLDAS'
    assert m.rs_source == 'NLDAS'
    assert m.ta_source == 'NLDAS'
    assert m.windspeed_source == 'NLDAS'
    assert m.topt_source == 'projects/openet/assets/ptjpl/ancillary/Topt_from_max_convolved'
    assert m.faparmax_source == 'projects/openet/assets/ptjpl/ancillary/fAPARmax'
    # assert m.et_reference_source is None
    # assert m.et_reference_band is None
    # assert m.et_reference_factor is None
    # assert m.et_reference_factor is None
    # assert m.latitude is None
    # assert m.longitude is None
    assert m.floor_Topt is True


# Todo: Break these up into separate functions?
def test_Image_init_calculated_properties():
    m = ptjpl.Image(default_image())
    assert utils.getinfo(m._time_start) == SCENE_TIME
    # assert m._scene_id.getInfo() == SCENE_ID
    # assert m._wrs2_tile.getInfo() == 'p{}r{}'.format(
    #     SCENE_ID.split('_')[1][:3], SCENE_ID.split('_')[1][3:])


def test_Image_init_date_properties():
    m = ptjpl.Image(default_image())
    assert utils.getinfo(m._date)['value'] == SCENE_TIME
    assert utils.getinfo(m._year) == int(SCENE_DATE.split('-')[0])
    assert utils.getinfo(m._month) == int(SCENE_DATE.split('-')[1])
    assert utils.getinfo(m._start_date)['value'] == utils.millis(SCENE_0UTC_DT)
    assert utils.getinfo(m._end_date)['value'] == utils.millis(
        SCENE_0UTC_DT + datetime.timedelta(days=1))
    assert utils.getinfo(m._doy) == SCENE_DOY


# CGM - scene_id is not currently being used in the model
# def test_Image_init_scene_id_property():
#     """Test that the system:index from a merged collection is parsed"""
#     input_img = default_image()
#     m = ptjpl.Image(input_img.set('system:index', '1_2_' + SCENE_ID))
#     assert utils.getinfo(m._scene_id) == SCENE_ID


@pytest.mark.parametrize('variable', ['albedo', 'emissivity', 'LST', 'NDVI'])
def test_Image_init_variable_properties(variable):
    """Test the band name and if properties are set on the variable images"""
    output = utils.getinfo(getattr(default_image_obj(), variable))
    # This assumes the band name is the lowercase of the variable name
    assert output['bands'][0]['id'] == variable.lower()
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'band, time, interp_flag, xy, expected',
    [
        ['shortwave_radiation', SCENE_TIME, False, TEST_POINT, 933.5350],
        ['shortwave_radiation', SCENE_TIME, True, TEST_POINT, 900.3892],
        ['temperature', SCENE_TIME, False, TEST_POINT, 41.7300],
        ['temperature', SCENE_TIME, True, TEST_POINT, 41.4162],
        ['specific_humidity', SCENE_TIME, False, TEST_POINT, 0.0062153],
        ['specific_humidity', SCENE_TIME, True, TEST_POINT, 0.0062153],
        # Other NLDAS bands
        # ['longwave_radiation', SCENE_TIME, False, TEST_POINT, 0],
        # ['longwave_radiation', SCENE_TIME, True, TEST_POINT, 0],
        # ['wind_u', SCENE_TIME, False, TEST_POINT, 0],
        # ['wind_u', SCENE_TIME, True, TEST_POINT, 0],
        # ['wind_v', SCENE_TIME, False, TEST_POINT, 0],
        # ['wind_v', SCENE_TIME, True, TEST_POINT, 0],
    ]
)
def test_Image_nldas_interpolate(band, time, interp_flag, xy, expected, tol=0.0001):
    output = utils.point_image_value(
        ptjpl.Image.nldas_interpolate(band, ee.Date(time), interp_flag), xy
    )
    assert abs(output[band] - expected) <= tol


@pytest.mark.parametrize(
    'ea_source, xy, expected',
    [
        ['NLDAS', TEST_POINT, 1012.1051],
        # Check string/float constant values
        ['1000', TEST_POINT, 1000],
        [1000, TEST_POINT, 1000],
    ]
)
def test_Image_ea_sources(ea_source, xy, expected, tol=0.01):
    """Test getting Ea values for a single date at a real point"""
    m = default_image_obj(ea_source=ea_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.ea), xy)
    assert abs(output['ea'] - expected) <= tol


def test_Image_ea_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(ea_source='').ea)


@pytest.mark.parametrize(
    'LWin_source, xy, expected',
    [
        ['NLDAS', TEST_POINT, 441.6057],
        # Check string/float constant values
        ['440', TEST_POINT, 440],
        [440, TEST_POINT, 440],
    ]
)
def test_Image_LWin_sources(LWin_source, xy, expected, tol=0.01):
    """Test getting LWin values for a single date at a real point"""
    m = default_image_obj(LWin_source=LWin_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.LWin), xy)
    assert abs(output['LWin'] - expected) <= tol


def test_Image_LWin_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(LWin_source='').LWin)


@pytest.mark.parametrize(
    'rs_source, xy, expected',
    [
        ['NLDAS', TEST_POINT, 903.4149],
        # Check string/float constant values
        ['900.0', TEST_POINT, 900.0],
        [900.0, TEST_POINT, 900.0],
    ]
)
def test_Image_rs_sources(rs_source, xy, expected, tol=0.0001):
    """Test getting Ta values for a single date at a real point"""
    m = default_image_obj(rs_source=rs_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.rs), xy)
    assert abs(output['rs'] - expected) <= tol


def test_Image_rs_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(rs_source='').rs)


@pytest.mark.parametrize(
    'ta_source, xy, expected',
    [
        ['NLDAS', TEST_POINT, 314.3867],
        # Check string/float constant values
        ['315.0', TEST_POINT, 315.0],
        [315.0, TEST_POINT, 315.0],
    ]
)
def test_Image_ta_sources(ta_source, xy, expected, tol=0.0001):
    """Test getting Ta values for a single date at a real point"""
    m = default_image_obj(ta_source=ta_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.ta), xy)
    assert abs(output['Ta'] - expected) <= tol


def test_Image_ta_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(ta_source='').ta)


@pytest.mark.parametrize(
    'windspeed_source, xy, expected',
    [
        ['NLDAS', TEST_POINT, 3.0524198],
        # Check string/float constant values
        ['3.0', TEST_POINT, 3.0],
        [3.0, TEST_POINT, 3.0],
    ]
)
def test_Image_windspeed_sources(windspeed_source, xy, expected, tol=0.0001):
    """Test getting Ta values for a single date at a real point"""
    m = default_image_obj(windspeed_source=windspeed_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.U), xy)
    assert abs(output['U'] - expected) <= tol


def test_Image_windspeed_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(windspeed_source='').U)


@pytest.mark.parametrize(
    'topt_source, xy, expected',
    [
        # Test using asset ID string
        ['projects/openet/ptjpl/ancillary/Topt', TEST_POINT, 19.7467],
        ['projects/openet/ptjpl/ancillary/Topt_from_max', TEST_POINT, 15.7395],
        ['projects/openet/ptjpl/ancillary/Topt_from_max_convolved', TEST_POINT, 15.7395],
        # # Test using ee.Image object
        # [ee.Image('projects/openet/ptjpl/ancillary/Topt'), TEST_POINT, 16.79],
        # Check string/float constant values
        ['40.0', TEST_POINT, 40.0],
        [40.0, TEST_POINT, 40.0],
    ]
)
def test_Image_topt_sources(topt_source, xy, expected, tol=0.001):
    """Test getting Topt values for a single date at a real point"""
    m = default_image_obj(topt_source=topt_source, floor_Topt=False)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.Topt), xy)
    assert abs(output['Topt'] - expected) <= tol


def test_Image_topt_floor(topt_source=15, ta_source=40, tol=0.001):
    m = default_image_obj(topt_source=topt_source, ta_source=ta_source+273.15, floor_Topt=True)
    output = utils.point_image_value(ee.Image(m.Topt), TEST_POINT)
    assert abs(output['Topt'] - ta_source) <= tol


def test_Image_topt_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(topt_source='').Topt)


@pytest.mark.parametrize(
    'faparmax_source, xy, expected',
    [
        # Check string/float constant values
        ['0.6', TEST_POINT, 0.6],
        [0.6, TEST_POINT, 0.6],
        # Check using asset ID string
        ['projects/openet/ptjpl/ancillary/fAPARmax', TEST_POINT, 0.3340],
        # # Check using keywords
        # ['ASSET', [-120.113, 36.336], 0.6144],
        # # Check using ee.Image object
        # [ee.Image('projects/openet/ptjpl/ancillary/fAPARmax'), TEST_POINT, 0.6144],
    ]
)
def test_Image_faparmax_sources(faparmax_source, xy, expected, tol=0.0001):
    """Test getting fAPARmax values for a single date at a real point"""
    m = default_image_obj(faparmax_source=faparmax_source)
    # Uncomment to check values for other dates
    # m._start_date = ee.Date(start_date)
    # m._end_date =  m._start_date.advance(1, 'day')
    output = utils.point_image_value(ee.Image(m.fAPARmax), xy)
    assert abs(output['fAPARmax'] - expected) <= tol


def test_Image_faparmax_sources_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(faparmax_source='').fAPARmax)


# TODO: Add tests for all the added properties


# We probably don't need to test the properties and band names for all of them
def test_Image_SWin_properties():
    """"""
    output = utils.getinfo(default_image_obj().SWin)
    assert output['bands'][0]['id'] == 'SWin'


def test_Image_SWin_value():
    """Test that a non-zero value is returned for the default inputs"""
    output = utils.constant_image_value(ee.Image(default_image_obj().SWin))
    assert output['SWin'] > 0


def test_Image_LE_defaults(expected=517.674, tol=0.001):
    output = utils.constant_image_value(ee.Image(default_image_obj().LE))
    assert abs(output['LE'] - expected) <= tol


def test_Image_PET_defaults(expected=718.188, tol=0.001):
    output = utils.constant_image_value(ee.Image(default_image_obj().PET))
    assert abs(output['PET'] - expected) <= tol


def test_Image_ESI_defaults(expected=0.721, tol=0.001):
    output = utils.constant_image_value(ee.Image(default_image_obj().ESI))
    assert abs(output['ESI'] - expected) <= tol


def test_Image_et_properties():
    """Test band name and if properties are set on the image"""
    output = utils.getinfo(default_image_obj().et)
    assert output['bands'][0]['id'] == 'et'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_et_defaults(expected=6.128, tol=0.001):
    output = utils.constant_image_value(ee.Image(default_image_obj().et))
    assert abs(output['et'] - expected) <= tol


# CGM - Test if there are any conditions that should return nodata
# @pytest.mark.parametrize(
#     'albedo, emissivity, lst, ndvi, expected',
#     [
#         [0.2, 0.99, 300, 0.80, None],
#     ]
# )
# def test_Image_et_nodata(albedo, emissivity, lst, ndvi, expected):
#     output_img = default_image_obj(albedo=albedo, emissivity=emissivity, lst=lst, ndvi=ndvi)
#     output = utils.constant_image_value(ee.Image(output_img.et))
#     assert output['et'] is None


def test_Image_et_reference_properties():
    """Test if properties are set on the reference ET image"""
    output = utils.getinfo(default_image_obj().et_reference)
    assert output['bands'][0]['id'] == 'et_reference'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


@pytest.mark.parametrize(
    'source, band, factor, xy, expected',
    [
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 1, TEST_POINT, 12.9],
        ['IDAHO_EPSCOR/GRIDMET', 'etr', 0.85, TEST_POINT, 12.9 * 0.85],
        ['projects/openet/assets/reference_et/california/cimis/daily/v1',
         'etr', 1, TEST_POINT, 11.7893],
        ['projects/openet/reference_et/california/cimis/daily/v1',
         'etr', 1, TEST_POINT, 11.7893],
        ['projects/earthengine-legacy/assets/projects/openet/reference_et/california/cimis/daily/v1',
         'etr', 1, TEST_POINT, 11.7893],
        [10, 'FOO', 1, TEST_POINT, 10.0],
        [10, 'FOO', 0.85, TEST_POINT, 8.5],
    ]
)
def test_Image_et_reference_sources(source, band, factor, xy, expected, tol=0.001):
    """Test getting reference ET values for a single date at a real point"""
    output = utils.point_image_value(default_image_obj(
        et_reference_source=source, et_reference_band=band,
        et_reference_factor=factor).et_reference, xy)
    assert abs(output['et_reference'] - expected) <= tol


# # DEADBEEF - Current implementation does not use etr_source for computing etr
# def test_Image_etr_values(etr_source=15, etr_factor=0.85, tol=0.0001):
#     output = utils.constant_image_value(default_image_obj(
#         etr_source=etr_source, etr_factor=etr_factor).etr)
#     assert abs(output['constant'] - etr_source * etr_factor) <= tol


def test_Image_et_fraction_properties():
    """Test if properties are set on the ET fraction image"""
    output = utils.getinfo(default_image_obj().et_fraction)
    assert output['bands'][0]['id'] == 'et_fraction'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


# def test_Image_et_fraction_defaults(expected=0.721, tol=0.001):
#     output = utils.constant_image_value(ee.Image(default_image_obj().et_fraction))
#     assert abs(output['et_fraction'] - expected) <= tol


# # DEADBEEF - Current implementation does not use etr_source for computing etr
# def test_Image_etf_values(etr_source=15, etr_factor=0.85, expected=13.335, tol=0.0001):
#     output = utils.constant_image_value(default_image_obj(
#         etr_source=etr_source, etr_factor=etr_factor).etf)
#     assert abs(output['etf'] - expected / (etr_source * etr_factor)) <= tol


# CGM - Can't check mask, time, and calculate() until ET is working
def test_Image_mask_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().mask)
    assert output['bands'][0]['id'] == 'mask'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_mask_values():
    assert utils.constant_image_value(default_image_obj().mask)['mask'] == 1


def test_Image_time_properties():
    """Test if properties are set on the time image"""
    output = utils.getinfo(default_image_obj().time)
    assert output['bands'][0]['id'] == 'time'
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_time_values():
    """The time band should have the 0 UTC time in it for interpolation"""
    assert utils.constant_image_value(
        default_image_obj().time)['time'] == utils.millis(SCENE_0UTC_DT)


def test_Image_calculate_properties():
    """Test if properties are set on the output image"""
    output = utils.getinfo(default_image_obj().calculate(['ndvi']))
    assert output['properties']['system:index'] == SCENE_ID
    assert output['properties']['system:time_start'] == SCENE_TIME
    assert output['properties']['image_id'] == COLL_ID + SCENE_ID


def test_Image_calculate_variables_default():
    output = utils.getinfo(default_image_obj().calculate())
    assert set([x['id'] for x in output['bands']]) == {'et'}


def test_Image_calculate_variables_custom():
    variables = {'ndvi'}
    output = utils.getinfo(default_image_obj().calculate(variables))
    assert set([x['id'] for x in output['bands']]) == variables


def test_Image_calculate_variables_all():
    variables = {'et', 'et_fraction', 'et_reference', 'mask', 'ndvi', 'time'}
    output = utils.getinfo(default_image_obj().calculate(variables=list(variables)))
    assert set([x['id'] for x in output['bands']]) == variables


def test_Image_calculate_values():
    """Test if the calculate method returns values"""
    output_img = default_image_obj().calculate(['et'])
    # output_img = default_image_obj().calculate(['et', 'et_reference', 'et_fraction'])
    assert utils.constant_image_value(output_img.select(['et']))['et'] > 0
    # assert utils.constant_image_value(output_img.select(['et_reference']))['et_reference'] > 0
    # assert utils.constant_image_value(output_img.select(['et_fraction']))['et_fraction'] > 0


def test_Image_calculate_variables_valueerror():
    """Test if calculate method raises a valueerror for invalid variables"""
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj().calculate(['FOO']))


def test_Image_from_landsat_c2_sr_default_image():
    """Test that the classmethod is returning a class object"""
    output = ptjpl.Image.from_landsat_c2_sr(
        ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828')
        # 'LANDSAT/LC08/C02/T1_L2/LC08_038031_20130828'
    )
    assert type(output) == type(default_image_obj())


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LT04/C02/T1_L2/LT04_044033_19830812',
        'LANDSAT/LT05/C02/T1_L2/LT05_044033_20110716',
        'LANDSAT/LE07/C02/T1_L2/LE07_044033_20170708',
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_landsat_c2_sr_landsat_image(image_id):
    """Test instantiating the class from a real Landsat images"""
    output = utils.getinfo(ptjpl.Image.from_landsat_c2_sr(ee.Image(image_id)).NDVI)
    assert output['properties']['system:index'] == image_id.split('/')[-1]


def test_Image_from_landsat_c2_sr_exception():
    """Test instantiating the class for an invalid image ID"""
    with pytest.raises(Exception):
        # Intentionally using .getInfo() since utils.getinfo() might catch the exception
        ptjpl.Image.from_landsat_c2_sr(ee.Image('FOO')).ndvi.getInfo()


def test_Image_from_landsat_c2_sr_scaling():
    """Test if Landsat SR images images are being scaled"""
    sr_img = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
    input_img = (
        ee.Image.constant([10909, 10909, 10909, 10909, 10909, 10909, 44177.6, 21824])
        .rename(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'])
        .set({'SPACECRAFT_ID': ee.String(sr_img.get('SPACECRAFT_ID')),
              'system:id': ee.String(sr_img.get('system:id')),
              'system:index': ee.String(sr_img.get('system:index')),
              'system:time_start': ee.Number(sr_img.get('system:time_start'))})
    )

    # LST correction and cloud score masking do not work with a constant image
    #   and must be explicitly set to False
    output = utils.constant_image_value(ptjpl.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False,
        cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}).albedo)
    assert abs(output['albedo'] - 0.1) <= 0.01

    output = utils.constant_image_value(ptjpl.Image.from_landsat_c2_sr(
        input_img, c2_lst_correct=False,
        cloudmask_args={'cloud_score_flag': False, 'filter_flag': False}).LST)
    assert abs(output['lst'] - 300) <= 0.1


def test_Image_from_landsat_c2_sr_cloud_mask_args():
    """Test if the cloud_mask_args parameter can be set (not if it works)"""
    output = ptjpl.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        cloudmask_args={'snow_flag': True, 'cirrus_flag': True}
    )
    assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_cloud_score_mask_arg():
    """Test if the cloud_score_flag parameter can be set in cloudmask_args"""
    output = ptjpl.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        cloudmask_args={'cloud_score_flag': True, 'cloud_score_pct': 50})
    assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_arg():
    """Test if the c2_lst_correct parameter can be set (not if it works)"""
    output = ptjpl.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702', c2_lst_correct=True)
    assert type(output) == type(default_image_obj())


def test_Image_from_landsat_c2_sr_c2_lst_correct_fill():
    """Test if the c2_lst_correct fills the LST holes in Nebraska"""
    image_id = 'LANDSAT/LC08/C02/T1_L2/LC08_031034_20160702'
    xy = [-102.08284, 37.81728]

    # CGM - Is the uncorrected test needed?
    uncorrected = utils.point_image_value(
        ptjpl.Image.from_landsat_c2_sr(image_id, c2_lst_correct=False).LST, xy)
    assert uncorrected['lst'] is None
    corrected = utils.point_image_value(
        ptjpl.Image.from_landsat_c2_sr(image_id, c2_lst_correct=True).LST, xy)
    assert corrected['lst'] > 0
    # # Exact test values copied from openet-core
    # assert abs(corrected['lst'] - 306.83) <= 0.25


@pytest.mark.parametrize(
    'image_id',
    [
        'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716',
        'LANDSAT/LC09/C02/T1_L2/LC09_044033_20220127',
    ]
)
def test_Image_from_image_id(image_id):
    """Test instantiating the class using the from_image_id method"""
    output = utils.getinfo(ptjpl.Image.from_image_id(image_id).NDVI)
    assert output['properties']['system:index'] == image_id.split('/')[-1]
    assert output['properties']['image_id'] == image_id


def test_Image_from_method_kwargs():
    """Test that the init parameters can be passed through the helper methods"""
    assert ptjpl.Image.from_landsat_c2_sr(
        'LANDSAT/LC08/C02/T1_L2/LC08_042035_20150713',
        ea_source='FOO').ea_source == 'FOO'


# TODO: Move the following tests to the correct spot in order of functions
# CGM - Copied from SIMS test_b_model.py
def test_Model_crop_type_source_exception():
    with pytest.raises(ValueError):
        # Intentionally using .getInfo() since utils.getinfo() might catch the exception
        default_image_obj(crop_type_source='FOO').crop_type.getInfo()


def test_Model_crop_type_constant_value():
    output = utils.constant_image_value(default_image_obj(crop_type_source=10).crop_type)
    assert output['crop_type'] == 10


@pytest.mark.parametrize(
    'year, expected',
    [
        [2007, 2008],
        [2008, 2008],
        [2016, 2016],
        # [2022, 2022],
        # [2023, 2022],
    ]
)
def test_Model_crop_type_source_cdl_collection(year, expected):
    """Test that the CDL collection is filtered to a single year and is limited
    to years with data
    """
    image = default_image_obj(crop_type_source='USDA/NASS/CDL')
    image._year = year
    output = utils.getinfo(image.crop_type)
    assert output['properties']['id'] == f'USDA/NASS/CDL/{expected}'


def test_Model_crop_type_source_cdl_image():
    output = utils.getinfo(default_image_obj(crop_type_source='USDA/NASS/CDL/2008').crop_type)
    assert output['properties']['id'] == 'USDA/NASS/CDL/2008'


def test_Model_crop_type_source_cdl_image_exception():
    """Requesting a CDL image that doesn't exist should raise an EE exception"""
    with pytest.raises(Exception):
        # Intentionally using .getInfo() since utils.getinfo() might catch the exception
        default_image_obj(crop_type_source='USDA/NASS/CDL/2099').crop_type.getInfo()


@pytest.mark.parametrize(
    'crop_type_source',
    [
        'projects/openet/assets/crop_type/v2023a',
        'projects/openet/assets/crop_type/v2021a',
        # 'projects/openet/crop_type/v2021a',
        # 'projects/earthengine-legacy/assets/projects/openet/crop_type/v2021a',
    ]
)
def test_Model_crop_type_source_openet_crop_type(crop_type_source):
    output = utils.getinfo(default_image_obj(crop_type_source=crop_type_source).crop_type)
    expected = crop_type_source.replace('projects/earthengine-legacy/assets/', '')
    assert output['properties']['id'] == expected


@pytest.mark.parametrize(
    'crop_type_source, xy, expected',
    [
        ['USDA/NASS/CDL', TEST_POINT, 36],
        ['USDA/NASS/CDL/2016', TEST_POINT, 36],
        ['projects/openet/assets/crop_type/v2023a', TEST_POINT, 47],
        ['projects/openet/assets/crop_type/v2021a', TEST_POINT, 47],
        # ['projects/openet/crop_type/v2021a', TEST_POINT, 47],
        # ['projects/earthengine-legacy/assets/projects/openet/crop_type/v2021a', TEST_POINT, 47],
    ]
)
def test_Model_crop_type_values(crop_type_source, xy, expected):
    output = utils.point_image_value(
        default_image_obj(crop_type_source=crop_type_source).crop_type, xy
    )
    assert output['crop_type'] == expected


@pytest.mark.parametrize(
    'crop_type_source, expected',
    [
        [36, 1],   # alfalfa
        [87, 1],   # wetland
        [190, 1],  # wetland
        [195, 1],  # wetland
        [152, 0],  # shrubland
    ]
)
def test_Model_crop_mask_constant_value(crop_type_source, expected):
    output = utils.constant_image_value(
        default_image_obj(crop_type_source=crop_type_source).crop_mask
    )
    assert output['crop_mask'] == expected


def test_Model_crop_pm_adjust_source_exception():
    with pytest.raises(ValueError):
        utils.getinfo(default_image_obj(crop_pm_adjust_source='FOO').crop_pm_adjust)


def test_Model_crop_pm_adjust_constant_value():
    output = utils.constant_image_value(
        default_image_obj(crop_pm_adjust_source=2, crop_type_source=36).crop_pm_adjust
    )
    assert output['crop_pm_adjust'] == 2


@pytest.mark.parametrize(
    'crop_pm_adjust_source, xy, expected',
    [
        ['projects/openet/assets/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
         TEST_POINT, 1.3024652077503645],
        # ['projects/openet/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
        #  TEST_POINT, 1.3024652077503645],
    ]
)
def test_Model_crop_pm_adjust_values(crop_pm_adjust_source, xy, expected, tol=0.001):
    output = utils.point_image_value(
        default_image_obj(crop_pm_adjust_source=crop_pm_adjust_source).crop_pm_adjust, xy
    )
    assert abs(output['crop_pm_adjust'] - expected) <= tol


# def test_Model_crop_pm_adjust_band():
#     output = utils.constant_image_value(default_image_obj(
#         crop_pm_adjust_source=2, crop_type_source=36).crop_pm_adjust)
#     assert output['crop_pm_adjust'] == 2


@pytest.mark.parametrize(
    'crop_pm_adjust_flag, crop_pm_adjust_source, xy, expected',
    [
        [True, 2, TEST_POINT, 2 * 6.128485756837151],
        [False, 2, TEST_POINT, 6.128485756837151],
        [True, 'projects/openet/assets/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
         TEST_POINT, 6.128485756837151 * 1.302],
        [False, 'projects/openet/assets/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
         TEST_POINT, 6.128485756837151],
        # [True, 'projects/openet/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
        #  TEST_POINT, 6.128485756837151 * 1.302],
        # [False, 'projects/openet/ptjpl/ancillary/alpha/gridmet_1980-2020_dgs',
        #  TEST_POINT, 6.128485756837151],
    ]
)
def test_Model_et_crop_pm_adjust_flag(crop_pm_adjust_flag, crop_pm_adjust_source,
                                      xy, expected, tol=0.01):
    output = utils.point_image_value(default_image_obj(
        crop_pm_adjust_flag=crop_pm_adjust_flag,
        crop_pm_adjust_source=crop_pm_adjust_source).et, xy)
    assert abs(output['et'] - expected) <= tol
