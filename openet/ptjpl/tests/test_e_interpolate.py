# import pprint

import ee
import pytest

import openet.ptjpl.interpolate as interpolate
import openet.ptjpl.utils as utils


def scene_coll(variables, et_fraction=0.4, et=5, ndvi=0.6):
    """Return a generic scene collection to test scene interpolation functions

    Parameters
    ----------
    variables : list
        The variables to return in the collection
    et_fraction : float
    et : float
    ndvi : float

    Returns
    -------
    ee.ImageCollection

    """
    img = (
        ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
        .select(['SR_B3']).double().multiply(0)
    )
    mask = img.add(1).updateMask(1).uint8()

    # The "date" is used for the time band since it needs to be the 0 UTC time
    # The "time" is advanced to match the typical Landsat overpass time
    date1 = ee.Number(ee.Date.fromYMD(2017, 7, 8).millis())
    date2 = ee.Number(ee.Date.fromYMD(2017, 7, 16).millis())
    date3 = ee.Number(ee.Date.fromYMD(2017, 7, 24).millis())
    time1 = ee.Number(ee.Date.fromYMD(2017, 7, 8).advance(18, 'hours').millis())
    time2 = ee.Number(ee.Date.fromYMD(2017, 7, 16).advance(18, 'hours').millis())
    time3 = ee.Number(ee.Date.fromYMD(2017, 7, 24).advance(18, 'hours').millis())

    # Mask and time bands currently get added on to the scene collection
    #   and images are unscaled just before interpolating in the export tool
    scene_img = (
        ee.Image([img.add(et_fraction), img.add(et), img.add(ndvi), mask])
        .rename(['et_fraction', 'et', 'ndvi', 'mask'])
    )
    scene_coll = ee.ImageCollection.fromImages([
        scene_img.addBands([img.add(date1).rename('time')])
            .set({'system:index': 'LE07_044033_20170708',
                  'system:time_start': time1,
                  'scale_factor': 1.0}),
        scene_img.addBands([img.add(date2).rename('time')])
            .set({'system:index': 'LC08_044033_20170716',
                  'system:time_start': time2,
                  'scale_factor': 1.0}),
        scene_img.addBands([img.add(date3).rename('time')])
            .set({'system:index': 'LE07_044033_20170724',
                  'system:time_start': time3,
                  'scale_factor': 1.0}),
    ])
    return scene_coll.select(variables)


def test_from_scene_et_fraction_t_interval_daily_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_fraction(
        scene_coll(['et_fraction', 'ndvi', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'ndvi'],
        interp_args={'interp_method': 'linear', 'interp_days': 32},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='daily')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['ndvi']['2017-07-10'] - 0.6) <= tol
    assert abs(output['et_fraction']['2017-07-10'] - 0.4) <= tol
    assert abs(output['et_reference']['2017-07-10'] - 8.0) <= tol
    assert abs(output['et']['2017-07-10'] - (8.0 * 0.4)) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 0.4) <= tol
    assert abs(output['et_fraction']['2017-07-31'] - 0.4) <= tol
    assert '2017-08-01' not in output['et_fraction'].keys()
    # assert output['count']['2017-07-01'] == 3


def test_from_scene_et_fraction_t_interval_monthly_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_fraction(
        scene_coll(['et_fraction', 'ndvi', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'ndvi', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='monthly')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['ndvi']['2017-07-01'] - 0.6) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 0.4) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.5) <= tol
    assert abs(output['et']['2017-07-01'] - (236.5 * 0.4)) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_fraction_t_interval_custom_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_fraction(
        scene_coll(['et_fraction', 'ndvi', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'ndvi', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='custom')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['ndvi']['2017-07-01'] - 0.6) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 0.4) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.5) <= tol
    assert abs(output['et']['2017-07-01'] - (236.5 * 0.4)) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_actual_t_interval_daily_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     # 'interp_factor': 1.0,
                     'interp_resample': 'nearest'},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='daily')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et_fraction']['2017-07-10'] - 0.5970309972763062) <= tol
    assert abs(output['et_reference']['2017-07-10'] - 8) <= tol
    assert abs(output['et']['2017-07-10'] - 4.776247978210449) <= tol
    assert abs(output['et']['2017-07-01'] - 3.988095283508301) <= tol
    assert abs(output['et']['2017-07-31'] - 5.0) <= tol
    assert '2017-08-01' not in output['et'].keys()


def test_from_scene_et_actual_t_interval_monthly_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     # 'interp_factor': 1.0,
                     'interp_resample': 'nearest'},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='monthly')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et']['2017-07-01'] - 145.9705047607422) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.5) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 145.9705047607422 / 236.5) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_actual_t_interval_custom_values(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     # 'interp_factor': 1.0,
                     'interp_resample': 'nearest'},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'nearest'},
        t_interval='custom')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et']['2017-07-01'] - 145.9705047607422) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.5) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 145.9705047607422 / 236.5) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_actual_t_interval_monthly_et_reference_factor(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     # 'interp_factor': 1.0,
                     'interp_resample': 'nearest'},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 0.5,
                    'et_reference_resample': 'nearest'},
        t_interval='monthly')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et']['2017-07-01'] - 145.9705047607422) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.5 * 0.5) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 145.970505 / 236.5 / 0.5) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_actual_t_interval_monthly_et_reference_resample(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask']),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction', 'count'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     'interp_factor': 1.0,
                     'interp_resample': 'bilinear'},
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_factor': 1.0,
                    'et_reference_resample': 'bilinear'},
        t_interval='monthly')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et']['2017-07-01'] - 145.86253356933594) <= tol
    assert abs(output['et_reference']['2017-07-01'] - 236.0560913) <= tol
    assert abs(output['et_fraction']['2017-07-01'] - 145.8625336 / 236.0560913) <= tol
    assert output['count']['2017-07-01'] == 3


def test_from_scene_et_actual_t_interval_daily_et_fraction_max(tol=0.0001):
    output_coll = interpolate.from_scene_et_actual(
        scene_coll(['et', 'time', 'mask'], et=100),
        start_date='2017-07-01', end_date='2017-08-01',
        variables=['et', 'et_reference', 'et_fraction'],
        interp_args={'interp_method': 'linear', 'interp_days': 32,
                     'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                     'interp_band': 'eto',
                     # 'interp_factor': 1.0,
                     'interp_resample': 'nearest',
                     'et_fraction_min': 0.0,
                     'et_fraction_max': 1.4,
                     },
        model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                    'et_reference_band': 'eto',
                    'et_reference_resample': 'nearest',
                    'et_reference_factor': 1.0,
                    },
        t_interval='daily')

    TEST_POINT = (-121.5265, 38.7399)
    output = utils.point_coll_value(output_coll, TEST_POINT, scale=30)
    assert abs(output['et_fraction']['2017-07-10'] - 1.4) <= tol


def test_from_scene_et_fraction_t_interval_bad_value():
    # Function should raise a ValueError if t_interval is not supported
    with pytest.raises(ValueError):
        interpolate.from_scene_et_fraction(
            scene_coll(['et', 'time', 'mask']),
            start_date='2017-07-01', end_date='2017-08-01', variables=['et'],
            interp_args={'interp_method': 'linear', 'interp_days': 32},
            model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                        'et_reference_band': 'etr',
                        'et_reference_factor': 0.5,
                        'et_reference_resample': 'nearest'},
            t_interval='deadbeef',
        )


def test_from_scene_et_fraction_t_interval_no_value():
    # Function should raise an Exception if t_interval is not set
    with pytest.raises(TypeError):
        interpolate.from_scene_et_fraction(
            scene_coll(['et', 'time', 'mask']),
            start_date='2017-07-01', end_date='2017-08-01',
            variables=['et', 'et_reference', 'et_fraction', 'count'],
            interp_args={'interp_method': 'linear', 'interp_days': 32},
            model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                        'et_reference_band': 'etr',
                        'et_reference_factor': 0.5,
                        'et_reference_resample': 'nearest'},
        )


def test_from_scene_et_actual_t_interval_bad_value():
    # Function should raise a ValueError if t_interval is not supported
    with pytest.raises(ValueError):
        interpolate.from_scene_et_actual(
            scene_coll(['et', 'time', 'mask']),
            start_date='2017-07-01', end_date='2017-08-01', variables=['et'],
            interp_args={'interp_method': 'linear', 'interp_days': 32,
                         'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                         'interp_band': 'etr', 'interp_resample': 'nearest'},
            model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                        'et_reference_band': 'etr',
                        'et_reference_factor': 1.0,
                        'et_reference_resample': 'nearest'},
            t_interval='deadbeef',
        )


def test_from_scene_et_actual_t_interval_no_value():
    # Function should raise an Exception if t_interval is not set
    with pytest.raises(TypeError):
        interpolate.from_scene_et_actual(
            scene_coll(['et', 'time', 'mask']),
            start_date='2017-07-01', end_date='2017-08-01', variables=['et'],
            interp_args={'interp_method': 'linear', 'interp_days': 32,
                         'interp_source': 'IDAHO_EPSCOR/GRIDMET',
                         'interp_band': 'etr',
                         'interp_resample': 'nearest'},
            model_args={'et_reference_source': 'IDAHO_EPSCOR/GRIDMET',
                        'et_reference_band': 'etr',
                        'et_reference_factor': 1.0,
                        'et_reference_resample': 'nearest'},
        )
