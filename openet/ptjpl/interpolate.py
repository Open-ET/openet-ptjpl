import datetime
import logging

from dateutil.relativedelta import *
import ee
import openet.core.interpolate
# TODO: import utils from openet.core
# import openet.core.utils as utils

from . import utils


def from_scene_et_fraction(scene_coll, start_date, end_date, variables,
                           interp_args, model_args, t_interval,
                           use_joins=False):
    """

    Parameters
    ----------
    scene_coll : ee.ImageCollection
        Non-daily ET fraction images that will be interpolated.
    start_date : str
        ISO format start date.
    end_date : str
        ISO format end date (exclusive, passed directly to .filterDate()).
    variables : list
        List of variables that will be returned in the Image Collection.
    interp_args : dict
        Parameters from the INTERPOLATE section of the INI file.
        interp_method : {'linear}, optional
            Interpolation method.  The default is 'linear'.
        interp_days : int, str, optional
            Number of extra days before the start date and after the end date
            to include in the interpolation calculation. The default is 32.
    model_args : dict
        Parameters from the MODEL section of the INI file.
    t_interval : {'daily', 'monthly', 'annual', 'custom'}
        Time interval over which to interpolate and aggregate values
        The 'custom' interval will aggregate all days within the start and end
        dates into an image collection with a single image.
    use_joins : bool, optional
        If True, use joins to link the target and source collections.

    Returns
    -------
    ee.ImageCollection

    Raises
    ------
    ValueError

    """
    # Get interp_method
    if 'interp_method' in interp_args.keys():
        interp_method = interp_args['interp_method']
    else:
        interp_method = 'linear'
        logging.debug('interp_method was not set, default to "linear"')

    # Get interp_days
    if 'interp_days' in interp_args.keys():
        interp_days = interp_args['interp_days']
    else:
        interp_days = 32
        logging.debug('interp_days was not set, default to 32')

    # Check that the input parameters are valid
    if t_interval.lower() not in ['daily', 'monthly', 'annual', 'custom']:
        raise ValueError(f'unsupported t_interval: {t_interval}')
    elif interp_method.lower() not in ['linear']:
        raise ValueError(f'unsupported interp_method: {interp_method}')

    if ((type(interp_days) is str or type(interp_days) is float) and
            utils.is_number(interp_days)):
        interp_days = int(interp_days)
    elif not type(interp_days) is int:
        raise TypeError('interp_days must be an integer')
    elif interp_days <= 0:
        raise ValueError('interp_days must be a positive integer')

    if not variables:
        raise ValueError('variables parameter must be set')

    # Adjust start/end dates based on t_interval
    # Increase the date range to fully include the time interval
    start_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    if t_interval.lower() == 'annual':
        start_dt = datetime.datetime(start_dt.year, 1, 1)
        # Covert end date to inclusive, flatten to beginning of year,
        # then add a year which will make it exclusive
        end_dt -= relativedelta(days=+1)
        end_dt = datetime.datetime(end_dt.year, 1, 1)
        end_dt += relativedelta(years=+1)
    elif t_interval.lower() == 'monthly':
        start_dt = datetime.datetime(start_dt.year, start_dt.month, 1)
        end_dt -= relativedelta(days=+1)
        end_dt = datetime.datetime(end_dt.year, end_dt.month, 1)
        end_dt += relativedelta(months=+1)
    start_date = start_dt.strftime('%Y-%m-%d')
    end_date = end_dt.strftime('%Y-%m-%d')

    # The start/end date for the interpolation include more days
    # (+/- interp_days) than are included in the ETr collection
    interp_start_dt = start_dt - datetime.timedelta(days=interp_days)
    interp_end_dt = end_dt + datetime.timedelta(days=interp_days)
    interp_start_date = interp_start_dt.date().isoformat()
    interp_end_date = interp_end_dt.date().isoformat()

    # Get reference ET parameters
    if 'et_reference_source' in model_args.keys():
        et_reference_source = model_args['et_reference_source']
    else:
        raise ValueError('et_reference_source was not set')

    if 'et_reference_band' in model_args.keys():
        et_reference_band = model_args['et_reference_band']
    else:
        raise ValueError('et_reference_band was not set')

    if 'et_reference_factor' in model_args.keys():
        et_reference_factor = model_args['et_reference_factor']
    else:
        et_reference_factor = 1.0
        logging.debug('et_reference_factor was not set, default to 1.0')
        # raise ValueError('et_reference_factor was not set')

    # CGM - Resampling is not working correctly so commenting out for now
    # if 'et_reference_resample' in model_args.keys():
    #     et_reference_resample = model_args['et_reference_resample']
    # else:
    #     et_reference_resample = 'nearest'
    #     logging.debug('et_reference_resample was not set, default to nearest')
    #     # raise ValueError('et_reference_resample was not set')

    if type(et_reference_source) is str:
        # Assume a string source is an single image collection ID
        #   not an list of collection IDs or ee.ImageCollection
        daily_et_ref_coll = ee.ImageCollection(et_reference_source) \
            .filterDate(start_date, end_date) \
            .select([et_reference_band], ['et_reference'])
    # elif isinstance(et_reference_source, computedobject.ComputedObject):
    #     # Interpret computed objects as image collections
    #     daily_et_reference_coll = et_reference_source \
    #         .filterDate(self.start_date, self.end_date) \
    #         .select([et_reference_band])
    else:
        raise ValueError(f'unsupported et_reference_source: {et_reference_source}')

    # Scale reference ET images (if necessary)
    # CGM - Resampling is not working correctly so not including for now
    if (et_reference_factor and et_reference_factor != 1):
        def et_reference_adjust(input_img):
            return input_img.multiply(et_reference_factor) \
                .set({'system:time_start': input_img.get('system:time_start')})
        daily_et_ref_coll = daily_et_ref_coll.map(et_reference_adjust)

    # Initialize variable list to only variables that can be interpolated
    interp_vars = ['et_fraction', 'ndvi']
    interp_vars = list(set(interp_vars) & set(variables))

    # To return ET, the ETf must be interpolated
    if 'et' in variables and 'et_fraction' not in interp_vars:
        interp_vars.append('et_fraction')

    # With the current interpolate.daily() function,
    #   something has to be interpolated in order to return et_reference
    if 'et_reference' in variables and 'et_fraction' not in interp_vars:
        interp_vars.append('et_fraction')

    # The time band is always needed for interpolation
    interp_vars.append('time')

    def interpolate_prep(img):
        """Prep WRS2 scene images for interpolation

        "Unscale" the images using the "scale_factor" property and
            convert to double.
        Add a mask and time band to each image in the scene_coll since
            interpolator is assuming time and mask bands exist.
        The interpolation could be modified to get the mask from the
            time band instead of setting it here.

        """
        mask_img = img.select(['et_fraction'], ['mask']) \
            .multiply(0).add(1).updateMask(1).uint8()
        time_img = img.select(['et_fraction'], ['time']) \
            .double().multiply(0) \
            .add(utils.date_to_time_0utc(ee.Date(img.get('system:time_start'))))

        # TODO: This should probably check if the scale_factor property
        # TODO:   exists and is not 1
        return img.select(interp_vars) \
            .double().multiply(ee.Number(img.get('scale_factor'))) \
            .addBands([mask_img, time_img]) \
            .set({'system:time_start': ee.Number(img.get('system:time_start'))})
        #     .set({'image_id': ee.String(img.get('system:index'))})

    # Filter scene collection to the interpolation range
    # This may not be needed since scene_coll was built to this range
    scene_coll = scene_coll \
        .filterDate(interp_start_date, interp_end_date) \
        .map(interpolate_prep)

    # For count, compute the composite/mosaic image for the mask band only
    if 'count' in variables:
        aggregate_coll = openet.core.interpolate.aggregate_to_daily(
            image_coll=scene_coll.select(['mask']),
            start_date=start_date, end_date=end_date)
        # The following is needed because the aggregate collection can be
        #   empty if there are no scenes in the target date range but there
        #   are scenes in the interpolation date range.
        # Without this the count image will not be built but the other
        #   bands will be which causes a non-homogeneous image collection.
        aggregate_coll = aggregate_coll.merge(
            ee.Image.constant(0).rename(['mask'])
                .set({'system:time_start': ee.Date(start_date).millis()}))

    # Interpolate to a daily time step
    # NOTE: the daily function is not computing ET (ETf x ETr)
    #   but is returning the target (ETr) band
    daily_coll = openet.core.interpolate.daily(
        target_coll=daily_et_ref_coll,
        source_coll=scene_coll.select(interp_vars),
        interp_method=interp_method, interp_days=interp_days,
        use_joins=use_joins,
        compute_product=False
    )

    # Compute ET from ETf and ETr (if necessary)
    # The check for et_fraction is needed since it is back computed from ET and ETr
    # if 'et' in variables or 'et_fraction' in variables:
    def compute_et(img):
        """This function assumes ETr and ETf are present"""
        et_img = img.select(['et_fraction'])\
            .multiply(img.select(['et_reference']))
        return img.addBands(et_img.double().rename('et'))
    daily_coll = daily_coll.map(compute_et)

    def aggregate_image(agg_start_date, agg_end_date, date_format):
        """Aggregate the daily images within the target date range

        Parameters
        ----------
        agg_start_date: ee.Date, str
            Start date (inclusive).
        agg_end_date : ee.Date, str
            End date (exclusive).
        date_format : str
            Date format for system:index (uses EE JODA format).

        Returns
        -------
        ee.Image

        Notes
        -----
        Since this function takes multiple inputs it is being called
        for each time interval by separate mappable functions

        """
        if 'et' in variables or 'et_fraction' in variables:
            et_img = daily_coll.filterDate(agg_start_date, agg_end_date) \
                .select(['et']).sum()
        if 'et_reference' in variables or 'et_fraction' in variables:
            # Get reference ET from the original collection instead of from
            #   the interpolated daily collection.
            # et_reference_img = daily_coll.select(['et_reference']) \
            et_reference_img = daily_et_ref_coll \
                .filterDate(agg_start_date, agg_end_date) \
                .sum()

        image_list = []
        if 'et' in variables:
            image_list.append(et_img.float())
        if 'et_reference' in variables:
            image_list.append(et_reference_img.float())
        if 'et_fraction' in variables:
            # Compute average et fraction over the aggregation period
            image_list.append(
                et_img.divide(et_reference_img).rename(
                    ['et_fraction']).float())
        if 'ndvi' in variables:
            # Compute average ndvi over the aggregation period
            ndvi_img = daily_coll \
                .filterDate(agg_start_date, agg_end_date) \
                .mean().select(['ndvi']).float()
            image_list.append(ndvi_img)
        if 'count' in variables:
            count_img = aggregate_coll \
                .filterDate(agg_start_date, agg_end_date) \
                .select(['mask']).sum().rename('count').uint8()
            image_list.append(count_img)

        return ee.Image(image_list) \
            .set({
            'system:index': ee.Date(agg_start_date).format(date_format),
            'system:time_start': ee.Date(agg_start_date).millis()})
        #     .set(interp_properties)\

    # Combine input, interpolated, and derived values
    if t_interval.lower() == 'daily':
        def agg_daily(daily_img):
            # CGM - Double check that this time_start is a 0 UTC time.
            # It should be since it is coming from the interpolate source
            #   collection, but what if source is GRIDMET (+6 UTC)?
            agg_start_date = ee.Date(daily_img.get('system:time_start'))
            # CGM - This calls .sum() on collections with only one image
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'day'),
                date_format='YYYYMMdd')

        return ee.ImageCollection(daily_coll.map(agg_daily))

    elif t_interval.lower() == 'monthly':
        def month_gen(iter_start_dt, iter_end_dt):
            iter_dt = iter_start_dt
            # Conditional is "less than" because end date is exclusive
            while iter_dt < iter_end_dt:
                yield iter_dt.strftime('%Y-%m-%d')
                iter_dt += relativedelta(months=+1)

        month_list = ee.List(list(month_gen(start_dt, end_dt)))

        def agg_monthly(agg_start_date):
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'month'),
                date_format='YYYYMM')

        return ee.ImageCollection(month_list.map(agg_monthly))

    elif t_interval.lower() == 'annual':
        def year_gen(iter_start_dt, iter_end_dt):
            iter_dt = iter_start_dt
            while iter_dt < iter_end_dt:
                yield iter_dt.strftime('%Y-%m-%d')
                iter_dt += relativedelta(years=+1)

        year_list = ee.List(list(year_gen(start_dt, end_dt)))

        def agg_annual(agg_start_date):
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'year'),
                date_format='YYYY')

        return ee.ImageCollection(year_list.map(agg_annual))

    elif t_interval.lower() == 'custom':
        # Returning an ImageCollection to be consistent
        return ee.ImageCollection(aggregate_image(
            agg_start_date=start_date, agg_end_date=end_date,
            date_format='YYYYMMdd'))


def from_scene_et_actual(scene_coll, start_date, end_date, variables,
                         interp_args, model_args, t_interval,
                         use_joins=False):
    """

    Parameters
    ----------
    scene_coll : ee.ImageCollection
        Non-daily 'et' images that will be interpolated.
    start_date : str
        ISO format start date.
    end_date : str
        ISO format end date (exclusive, passed directly to .filterDate()).
    variables : list
        List of variables that will be returned in the Image Collection.
    interp_args : dict
        Parameters from the INTERPOLATE section of the INI file.
        # TODO: Look into a better format for showing the options
        interp_method : {'linear}, optional
            Interpolation method.  The default is 'linear'.
        interp_days : int, str, optional
            Number of extra days before the start date and after the end date
            to include in the interpolation calculation. The default is 32.
        use_joins : bool, optional
            If True, use joins to link the target and source collections.
            If False, the source collection will be filtered for each target image.
    model_args : dict
        Parameters from the MODEL section of the INI file.
    t_interval : {'daily', 'monthly', 'annual', 'custom'}
        Time interval over which to interpolate and aggregate values
        The 'custom' interval will aggregate all days within the start and end
        dates into an image collection with a single image.
    use_joins : bool, optional
        If True, use joins to link the target and source collections.
        If False, the source collection will be filtered for each target image.

    Returns
    -------
    ee.ImageCollection

    Raises
    ------
    ValueError

    """
    # Get interp_method
    if 'interp_method' in interp_args.keys():
        interp_method = interp_args['interp_method']
    else:
        interp_method = 'linear'
        logging.debug('interp_method was not set, default to "linear"')

    # Get interp_days
    if 'interp_days' in interp_args.keys():
        interp_days = interp_args['interp_days']
    else:
        interp_days = 32
        logging.debug('interp_days was not set, default to 32')

    # Override the use_joins function parameter if it was set in interp_args
    # Eventually the separate function parameter will be deprecated and removed
    if 'use_joins' in interp_args.keys():
        use_joins = interp_args['use_joins']
    # else:
    #     use_joins = False
    #     logging.debug('use_joins was not set, default to False')

    # Check that the input parameters are valid
    if t_interval.lower() not in ['daily', 'monthly', 'annual', 'custom']:
        raise ValueError(f'unsupported t_interval: {t_interval}')
    elif interp_method.lower() not in ['linear']:
        raise ValueError(f'unsupported interp_method: {interp_method}')

    if ((type(interp_days) is str or type(interp_days) is float) and
            utils.is_number(interp_days)):
        interp_days = int(interp_days)
    elif not type(interp_days) is int:
        raise TypeError('interp_days must be an integer')
    elif interp_days <= 0:
        raise ValueError('interp_days must be a positive integer')

    if not variables:
        raise ValueError('variables parameter must be set')

    # Adjust start/end dates based on t_interval
    # Increase the date range to fully include the time interval
    start_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_dt = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    if t_interval.lower() == 'annual':
        start_dt = datetime.datetime(start_dt.year, 1, 1)
        # Covert end date to inclusive, flatten to beginning of year,
        # then add a year which will make it exclusive
        end_dt -= relativedelta(days=+1)
        end_dt = datetime.datetime(end_dt.year, 1, 1)
        end_dt += relativedelta(years=+1)
    elif t_interval.lower() == 'monthly':
        start_dt = datetime.datetime(start_dt.year, start_dt.month, 1)
        end_dt -= relativedelta(days=+1)
        end_dt = datetime.datetime(end_dt.year, end_dt.month, 1)
        end_dt += relativedelta(months=+1)
    start_date = start_dt.strftime('%Y-%m-%d')
    end_date = end_dt.strftime('%Y-%m-%d')

    # The start/end date for the interpolation include more days
    # (+/- interp_days) than are included in the ETr collection
    interp_start_dt = start_dt - datetime.timedelta(days=interp_days)
    interp_end_dt = end_dt + datetime.timedelta(days=interp_days)
    interp_start_date = interp_start_dt.date().isoformat()
    interp_end_date = interp_end_dt.date().isoformat()

    # Get reference ET collection
    if 'et_reference' in variables or 'et_fraction' in variables:
        if 'et_reference_source' not in model_args.keys():
            raise ValueError('et_reference_source was not set')

        if 'et_reference_band' not in model_args.keys():
            raise ValueError('et_reference_band was not set')

        # TODO: Check if model_args can be modified instead of making new variables
        if 'et_reference_factor' in model_args.keys():
            et_reference_factor = model_args['et_reference_factor']
        else:
            et_reference_factor = 1.0
            logging.debug('et_reference_factor was not set, default to 1.0')
            # raise ValueError('et_reference_factor was not set')

        # CGM - Resampling is not working correctly so commenting out for now
        # if 'et_reference_resample' in model_args.keys():
        #     et_reference_resample = model_args['et_reference_resample']
        # else:
        #     et_reference_resample = 'nearest'
        #     logging.debug('et_reference_resample was not set, default to nearest')
        #     # raise ValueError('et_reference_resample was not set')

        # Assume a string source is an single image collection ID
        #   not an list of collection IDs or ee.ImageCollection
        daily_et_ref_coll_id = model_args['et_reference_source']
        daily_et_ref_coll = ee.ImageCollection(daily_et_ref_coll_id) \
            .filterDate(start_date, end_date) \
            .select([model_args['et_reference_band']], ['et_reference'])

        # Scale reference ET images (if necessary)
        # CGM - Resampling is not working correctly so not including for now
        if (et_reference_factor and et_reference_factor != 1):
            def et_reference_adjust(input_img):
                return input_img.multiply(et_reference_factor) \
                    .copyProperties(input_img) \
                    .set({'system:time_start': input_img.get('system:time_start')})
            daily_et_ref_coll = daily_et_ref_coll.map(et_reference_adjust)

    # Get the interpolation collection
    # TODO: Add logic to try falling back on et_reference parameters
    #   in model_args (if they are available)
    if 'interp_source' not in interp_args.keys():
        raise ValueError('interp_source was not set')
    if 'interp_band' not in interp_args.keys():
        raise ValueError('interp_band was not set')

    # CGM - Resampling is not working correctly so commenting out for now
    # if 'interp_resample' not in interp_args.keys():
    #     interp_args['interp_resample'] = 'nearest'
    #     logging.debug('interp_resample was not set, defaulting to nearest')
    #     # raise ValueError('interp_resample was not set')

    # CGM - Factor is not currently being applied so commenting out for now
    # if 'interp_factor' not in interp_args.keys():
    #     interp_args['interp_factor'] = 1.0
    #     logging.debug('interp_factor was not set, defaulting to 1.0')
    #     # raise ValueError('interp_factor was not set')

    # Target collection needs to be filtered to the same date range as the
    #   scene collection in order to normalize the scenes.
    # It will be filtered again to the start/end when it is sent into
    #   interpolate.daily()
    daily_target_coll = ee.ImageCollection(interp_args['interp_source']) \
        .filterDate(interp_start_date, interp_end_date) \
        .select([interp_args['interp_band']])

    # For count, compute the composite/mosaic image for the mask band only
    if 'count' in variables:
        aggregate_coll = openet.core.interpolate.aggregate_to_daily(
            image_coll=scene_coll.select(['mask']),
            start_date=start_date, end_date=end_date)
        # The following is needed because the aggregate collection can be
        #   empty if there are no scenes in the target date range but there
        #   are scenes in the interpolation date range.
        # Without this the count image will not be built but the other
        #   bands will be which causes a non-homogeneous image collection.
        aggregate_coll = aggregate_coll.merge(
            ee.Image.constant(0).rename(['mask'])
                .set({'system:time_start': ee.Date(start_date).millis()}))

    # It might be more efficient to join the target collection to the scenes
    def normalize_et(img):
        img_date = ee.Date(img.get('system:time_start')) \
            .update(hour=0, minute=0, second=0)
        img_date = ee.Date(img_date.millis().divide(1000).floor().multiply(1000))
        target_img = ee.Image(daily_target_coll \
            .filterDate(img_date, img_date.advance(1, 'day')).first())

        # CGM - This is causing weird artifacts in the output images
        # if interp_args['interp_resample'].lower() in ['bilinear', 'bicubic']:
        #     target_img = target_img.resample(interp_args['interp_resample'])

        et_norm_img = img.select(['et']).divide(target_img).rename(['et_norm'])

        # Clamp the normalized ET image (et_fraction)
        if 'et_fraction_max' in interp_args.keys():
            et_norm_img = et_norm_img.min(float(interp_args['et_fraction_max']))
        if 'et_fraction_min' in interp_args.keys():
            et_norm_img = et_norm_img.max(float(interp_args['et_fraction_min']))
        # if ('et_fraction_min' in interp_args.keys() and
        #     'et_fraction_max' in interp_args.keys()):
        #     et_norm_img = et_norm_img.clamp(
        #         float(interp_args['et_fraction_min']),
        #         float(interp_args['et_fraction_max']))

        return img.addBands([et_norm_img.double(), target_img.rename(['norm'])])

    # The time band is always needed for interpolation
    scene_coll = scene_coll \
        .filterDate(interp_start_date, interp_end_date) \
        .select(['et'] + ['mask', 'time']) \
        .map(normalize_et)

    # # Join the target (normalization) image to the scene images
    # if use_joins:
    #     prev_filter = ee.Filter.And(
    #         ee.Filter.maxDifference(
    #             difference=(interp_days + 1) * 24 * 60 * 60 * 1000,
    #             leftField='system:time_start', rightField='system:time_start'),
    #         ee.Filter.greaterThan(leftField='system:time_start',
    #                               rightField='system:time_start')
    #     )
    #     scene_coll = ee.ImageCollection(
    #         ee.Join.saveFirst(matchKey='norm_img', ordering='system:time_start',
    #                           ascending=False)
    #             .apply(primary=scene_coll, secondary=target_coll,
    #                    condition=prev_filter)
    #     )

    # Interpolate to a daily time step
    daily_coll = openet.core.interpolate.daily(
        target_coll=daily_target_coll.filterDate(start_date, end_date),
        source_coll=scene_coll.select(['et_norm', 'time']),
        interp_method=interp_method, interp_days=interp_days,
        use_joins=use_joins,
        compute_product=True,
    )

    # The interpolate.daily() function is currently returning the product of
    # the source and target image named as "{source_band}_1".
    # This approach will not be valid if other bands are interpolated.
    daily_coll = daily_coll.select(['et_norm_1'], ['et'])

    # Convert normalized ET back to ET
    # This isn't needed if compute_product=True in daily() and band is renamed
    # The check for et_fraction is needed since it is back computed from ET and ETr
    # # if 'et' in variables or 'et_fraction' in variables:
    # def compute_et(img):
    #     """This function assumes ETr and ETf are present"""
    #     et_img = img.select(['et_norm']).multiply(
    #         img.select(['et_reference']))
    #     return img.addBands(et_img.double().rename('et'))
    # daily_coll = daily_coll.map(compute_et)


    def aggregate_image(agg_start_date, agg_end_date, date_format):
        """Aggregate the daily images within the target date range

        Parameters
        ----------
        agg_start_date: ee.Date, str
            Start date (inclusive).
        agg_end_date : ee.Date, str
            End date (exclusive).
        date_format : str
            Date format for system:index (uses EE JODA format).

        Returns
        -------
        ee.Image

        Notes
        -----
        Since this function takes multiple inputs it is being called
        for each time interval by separate mappable functions

        """
        if 'et' in variables or 'et_fraction' in variables:
            et_img = daily_coll.filterDate(agg_start_date, agg_end_date) \
                .select(['et']).sum()
        if 'et_reference' in variables or 'et_fraction' in variables:
            # Get the reference ET image from the reference ET collection,
            #   not the interpolated collection
            et_reference_img = daily_et_ref_coll \
                .filterDate(agg_start_date, agg_end_date) \
                .sum()

        image_list = []
        if 'et' in variables:
            image_list.append(et_img.float())
        if 'et_reference' in variables:
            image_list.append(et_reference_img.float())
        if 'et_fraction' in variables:
            # Compute average et fraction over the aggregation period
            image_list.append(
                et_img.divide(et_reference_img)
                    .rename(['et_fraction']).float())
        # if 'ndvi' in variables:
        #     # Compute average ndvi over the aggregation period
        #     ndvi_img = daily_coll \
        #         .filterDate(agg_start_date, agg_end_date) \
        #         .mean().select(['ndvi']).float()
        #     image_list.append(ndvi_img)
        if 'count' in variables:
            count_img = aggregate_coll \
                .filterDate(agg_start_date, agg_end_date) \
                .select(['mask']).sum().rename('count').uint8()
            image_list.append(count_img)

        return ee.Image(image_list) \
            .set({
            'system:index': ee.Date(agg_start_date).format(date_format),
            'system:time_start': ee.Date(agg_start_date).millis()})
        #     .set(interp_properties)\

    # Combine input, interpolated, and derived values
    if t_interval.lower() == 'daily':
        def agg_daily(daily_img):
            # CGM - Double check that this time_start is a 0 UTC time.
            # It should be since it is coming from the interpolate source
            #   collection, but what if source is GRIDMET (+6 UTC)?
            agg_start_date = ee.Date(daily_img.get('system:time_start'))
            # CGM - This calls .sum() on collections with only one image
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'day'),
                date_format='YYYYMMdd')

        return ee.ImageCollection(daily_coll.map(agg_daily))

    elif t_interval.lower() == 'monthly':
        def month_gen(iter_start_dt, iter_end_dt):
            iter_dt = iter_start_dt
            # Conditional is "less than" because end date is exclusive
            while iter_dt < iter_end_dt:
                yield iter_dt.strftime('%Y-%m-%d')
                iter_dt += relativedelta(months=+1)

        month_list = ee.List(list(month_gen(start_dt, end_dt)))

        def agg_monthly(agg_start_date):
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'month'),
                date_format='YYYYMM')

        return ee.ImageCollection(month_list.map(agg_monthly))

    elif t_interval.lower() == 'annual':
        def year_gen(iter_start_dt, iter_end_dt):
            iter_dt = iter_start_dt
            while iter_dt < iter_end_dt:
                yield iter_dt.strftime('%Y-%m-%d')
                iter_dt += relativedelta(years=+1)

        year_list = ee.List(list(year_gen(start_dt, end_dt)))

        def agg_annual(agg_start_date):
            return aggregate_image(
                agg_start_date=agg_start_date,
                agg_end_date=ee.Date(agg_start_date).advance(1, 'year'),
                date_format='YYYY')

        return ee.ImageCollection(year_list.map(agg_annual))

    elif t_interval.lower() == 'custom':
        # Returning an ImageCollection to be consistent
        return ee.ImageCollection(aggregate_image(
            agg_start_date=start_date, agg_end_date=end_date,
            date_format='YYYYMMdd'))
