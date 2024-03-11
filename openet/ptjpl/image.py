import ee
import openet.core.common

from openet.ptjpl.daylight_hours import sha_deg_from_doy_lat, sunrise_from_sha, daylight_from_sha
from . import landsat
from . import ptjpl
from . import utils

DOWNSAMPLE_METHOD = "bilinear"

STEFAN_BOLTZMAN_CONSTANT = 5.67036713e-8
# KRN = 0.6
KPAR = 0.5


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated

    https://stevenloria.com/lazy-properties/
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


class Image:
    """Earth Engine based PT-JPL Image"""

    _C2_LST_CORRECT = True  # Enable (True) C2 LST correction to recalculate LST

    def __init__(
            self, image,
            windspeed_source='NLDAS',
            ea_source='NLDAS',
            rs_source='NLDAS',
            LWin_source='NLDAS',
            ta_source='NLDAS',
            topt_source='projects/openet/assets/ptjpl/ancillary/Topt_from_max_convolved',
            faparmax_source='projects/openet/assets/ptjpl/ancillary/fAPARmax',
            latitude=None,
            longitude=None,
            floor_Topt=True,
            **kwargs
    ):
        """Construct a generic PT-JPL Image

        Parameters
        ----------
        image : ee.Image
            A "prepped" PT-JPL input image.
            Bands: 'albedo', 'ndvi', 'lst'
            Properties: 'system:index', 'system:time_start'.
        ea_source : {'NLDAS'}, optional
            Actual vapor pressure source keyword (the default is 'NLDAS').
        rs_source : {'NLDAS'}, optional
            Incoming shortwave solar radiation source keyword
            (the default is 'NLDAS').
        ta_source : {'NLDAS'}, optional
            Air temperature source keyword (the default is 'NLDAS').
        topt_source : str, optional
            Optimal temperature source.
        faparmax_source : str, optional
            Maximum fraction of absorbed PAR source.
        latitude : ee.Image, ee.Number, float, optional
            Latitude [deg].  If not set will default to ee.Image.pixelLonLat()
            mapped to the input image.
        longitude : ee.Image, ee.Number, float, optional
            Longitude [deg].  If not set will default to ee.Image.pixelLonLat()
            mapped to the input image.
        floor_Topt : bool, optional
            ? (the default is True).
        kwargs : dict, optional
            et_reference_source : str, float
                Reference ET source (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_band : str
                Reference ET band name (the default is None).
                Parameter is required if computing 'et_fraction' or 'et_reference'.
            et_reference_factor : float, None
                Reference ET scaling factor.  The default is None which is
                equivalent to 1.0 (or no scaling).
            et_reference_resample : {'nearest', 'bilinear', 'bicubic', None}
                Reference ET resampling.  The default is None which is
                equivalent to nearest neighbor resampling.

        """
        self.image = ee.Image(image)

        # Set as "lazy_property" below in order to return custom properties
        # self.albedo = self.image.select('albedo')
        # self.emissivity = self.image.select('emissivity')
        # self.lst = self.image.select('lst')
        # self.ndvi = self.image.select('ndvi')

        # Copy system properties
        self._id = self.image.get('system:id')
        self._index = self.image.get('system:index')
        self._time_start = self.image.get('system:time_start')
        self._properties = {
            'system:index': self._index,
            'system:time_start': self._time_start,
            'image_id': self._id,
        }

        # # Build SCENE_ID from the (possibly merged) system:index
        # scene_id = ee.List(ee.String(self._index).split('_')).slice(-3)
        # self._scene_id = ee.String(scene_id.get(0)).cat('_') \
        #     .cat(ee.String(scene_id.get(1))).cat('_') \
        #     .cat(ee.String(scene_id.get(2)))

        # Build WRS2_TILE from the scene_id
        # self._wrs2_tile = ee.String('p').cat(self._scene_id.slice(5, 8)) \
        #     .cat('r').cat(self._scene_id.slice(8, 11))

        # Set server side date/time properties using the 'system:time_start'
        self._date = ee.Date(self._time_start)
        self._year = ee.Number(self._date.get('year'))
        self._month = ee.Number(self._date.get('month'))
        self._start_date = ee.Date(utils.date_to_time_0utc(self._date))
        self._end_date = self._start_date.advance(1, 'day')
        self._doy = ee.Number(self._date.getRelative('day', 'year')).add(1).int()

        # Model input parameters
        self.ea_source = ea_source
        self.rs_source = rs_source
        self.LWin_source = LWin_source
        self.ta_source = ta_source
        self.windspeed_source = windspeed_source
        self.topt_source = topt_source
        self.faparmax_source = faparmax_source

        # PM Crop Adjust
        try:
            self.crop_pm_adjust_flag = kwargs['crop_pm_adjust_flag']
        except:
            self.crop_pm_adjust_flag = False

        # CGM - Should the other parameters only be read if flag is true?
        # CGM - Should the crop type source and remap have defaults?
        # if self.crop_pm_adjust_flag:
        try:
            self.crop_pm_adjust_source = kwargs['crop_pm_adjust_source']
        except:
            self.crop_pm_adjust_source = None
        try:
            self.crop_pm_adjust_band = kwargs['crop_pm_adjust_band']
        except:
            self.crop_pm_adjust_band = None
        try:
            self.crop_type_source = kwargs['crop_type_source']
        except:
            self.crop_type_source = None
            # self.crop_type_source = 'USDA/NASS/CDL'
        try:
            self.crop_type_remap = kwargs['crop_type_remap']
        except:
            self.crop_type_remap = None
            # self.crop_type_remap = 'CDL'

        # Reference ET parameters
        try:
            self.et_reference_source = kwargs['et_reference_source']
        except:
            self.et_reference_source = None
        try:
            self.et_reference_band = kwargs['et_reference_band']
        except:
            self.et_reference_band = None
        try:
            self.et_reference_factor = kwargs['et_reference_factor']
        except:
            self.et_reference_factor = None
        try:
            self.et_reference_resample = kwargs['et_reference_resample']
        except:
            self.et_reference_resample = None

        # Check reference ET parameters
        if self.et_reference_factor and not utils.is_number(self.et_reference_factor):
            raise ValueError('et_reference_factor must be a number')
        if self.et_reference_factor and self.et_reference_factor < 0:
            raise ValueError('et_reference_factor must be greater than zero')
        resample_methods = ['nearest', 'bilinear', 'bicubic']
        if (self.et_reference_resample and
                self.et_reference_resample.lower() not in resample_methods):
            raise ValueError('unsupported et_reference_resample method')

        # Allow latitude and longitude to be set as a number to help with testing
        if latitude is None:
            self.latitude = self.NDVI.multiply(0).add(
                ee.Image.pixelLonLat().select(['latitude']))
        elif utils.is_number(latitude):
            self.latitude = ee.Image.constant(latitude)
        elif isinstance(latitude, ee.computedobject.ComputedObject):
            self.latitude = latitude
        else:
            raise ValueError('invalid latitude parameter')

        if longitude is None:
            self.longitude = self.NDVI.multiply(0).add(
                ee.Image.pixelLonLat().select(['longitude']))
        elif utils.is_number(longitude):
            self.longitude = ee.Image.constant(longitude)
        elif isinstance(longitude, ee.computedobject.ComputedObject):
            self.longitude = longitude
        else:
            raise ValueError('invalid longitude parameter')

        self.time_start = self.image.get('system:time_start')
        self.date = ee.Date(self.time_start)
        self.floor_Topt = floor_Topt

    def calculate(self, variables=['et']):
        """Return a multiband image of calculated variables

        Parameters
        ----------
        variables : list

        Returns
        -------
        ee.Image

        """
        output_images = []
        for v in variables:
            if v.lower() == 'et':
                output_images.append(self.et.float())
            elif v.lower() == 'et_fraction':
                output_images.append(self.et_fraction.float())
            elif v.lower() == 'et_reference':
                output_images.append(self.et_reference.float())
            elif v.lower() == 'lst':
                output_images.append(self.LST.float())
            elif v.lower() == 'mask':
                output_images.append(self.mask)
            elif v.lower() == 'ndvi':
                output_images.append(self.NDVI.float())
            # elif v.lower() == 'qa':
            #     output_images.append(self.qa)
            elif v.lower() == 'quality':
                output_images.append(self.quality)
            elif v.lower() == 'time':
                output_images.append(self.time)
            else:
                raise ValueError(f'unsupported variable: {v}')

        return ee.Image(output_images).set(self._properties)

    @lazy_property
    def et(self):
        """Compute PT-JPL actual ET [mm day-1]"""
        return self.ET.rename(['et']).set(self._properties)

    @lazy_property
    def et_fraction(self):
        """Fraction of reference ET (equivalent to the Kc)"""
        return self.ET.divide(self.et_reference).rename(['et_fraction']).set(self._properties)

    @lazy_property
    def et_reference(self):
        """Reference ET for the image date"""
        if utils.is_number(self.et_reference_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.et_reference_source)
            et_reference_img = ee.Image.constant(self.et_reference_source)
        elif type(self.et_reference_source) is str:
            # Assume a string source is an image collection ID (not an image ID)
            et_reference_coll = (
                ee.ImageCollection(self.et_reference_source)
                .filterDate(self._start_date, self._end_date)
                .select([self.et_reference_band])
            )
            et_reference_img = ee.Image(et_reference_coll.first())
            if self.et_reference_resample in ['bilinear', 'bicubic']:
                et_reference_img = et_reference_img.resample(self.et_reference_resample)
        else:
            raise ValueError(f'unsupported et_reference_source: {self.et_reference_source}')

        if self.et_reference_factor:
            et_reference_img = et_reference_img.multiply(self.et_reference_factor)

        return self.NDVI.multiply(0).add(et_reference_img) \
            .rename(['et_reference']).set(self._properties)

    @lazy_property
    def mask(self):
        """Mask of all active pixels (based on the final et)"""
        return self.et.multiply(0).add(1).updateMask(1) \
            .rename(['mask']).set(self._properties).uint8()

    @lazy_property
    def quality(self):
        """Set quality to 1 for all active pixels (for now)"""
        return self.mask.rename(['quality']).set(self._properties)

    @lazy_property
    def time(self):
        """Return an image of the 0 UTC time (in milliseconds)"""
        return self.mask \
            .double().multiply(0).add(utils.date_to_time_0utc(self._date)) \
            .rename(['time']).set(self._properties)

    @lazy_property
    def hour_utc(self):
        """Hour of day image"""
        return ee.Number(self._date.getFraction('day')).multiply(24.0)

    @lazy_property
    def utc_offset(self):
        return self.longitude.divide(180.0).multiply(12.0)

    @lazy_property
    def hour_solar(self):
        return self.utc_offset.add(self.hour_utc)

    @lazy_property
    def SHA(self):
        return sha_deg_from_doy_lat(self._doy, self.latitude)

    @lazy_property
    def sunrise_hour(self):
        return sunrise_from_sha(self.SHA)

    @lazy_property
    def daylight_hours(self):
        return daylight_from_sha(self.SHA)

    @lazy_property
    def albedo(self):
        """Albedo"""
        return self.image.select(['albedo']).set(self._properties)

    @lazy_property
    def emissivity(self):
        """Emissivity"""
        return self.image.select(['emissivity']).set(self._properties)

    @lazy_property
    def LST(self):
        """Land surface temperature (LST)"""
        return self.image.select(['lst']).set(self._properties)

    @lazy_property
    def ST_K(self):
        """Surface temperature Kelvin"""
        return self.LST.rename(['ST_K']).set(self._properties)

    @lazy_property
    def ST_C(self):
        """Surface temperature Celsius"""
        return self.ST_K.subtract(273.15).rename(['ST_C']).set(self._properties)

    @lazy_property
    def NDVI(self):
        """Normalized difference vegetation index (NDVI)"""
        return self.image.select(['ndvi']).set(self._properties)

    @lazy_property
    def NDWI(self):
        """Normalized Difference Water Index (NDWI)"""
        return self.image.select(['ndwi']).set(self._properties)

    @lazy_property
    def MNDWI(self):
        """Modified Normalized Difference Water Index (NDWI)"""
        return self.image.select(['mndwi']).set(self._properties)

    @lazy_property
    def WRI(self):
        """Water Ratio Index (WRI)"""
        return self.image.select(['wri']).set(self._properties)

    @lazy_property
    def water_mask(self):
        """Water pixel identification"""
        return (
            self.NDWI.gt(0)
            .And(self.MNDWI.gt(0))
            .And(self.WRI.gt(1))
            .And(self.NDVI.lt(0))
            .rename(['water_mask']).set(self._properties)
        )

    @lazy_property
    def SAVI(self):
        """Soil adjusted vegetation index (SAVI)"""
        return ptjpl.SAVI(self.NDVI).rename(['SAVI']).set(self._properties)

    @lazy_property
    def fAPAR(self):
        """Fraction of absorbed photosynthetically active radiation"""
        return ptjpl.fAPAR(self.SAVI).rename(['fAPAR']).set(self._properties)

    @lazy_property
    def fIPAR(self):
        """Fraction of intercepted photosynthetically active radiation"""
        return ptjpl.fIPAR(self.NDVI).rename(['fIPAR']).set(self._properties)

    @lazy_property
    def LAI(self):
        """Leaf area index (LAI) image"""
        return ptjpl.LAI(self.NDVI, self.fIPAR).rename(['LAI']).set(self._properties)

    @lazy_property
    def SWin(self):
        """Instantaneous incoming shortwave radiation in watts per square meter"""
        SWin = self.rs

        return SWin.rename(['SWin']).set(self._properties)

    @lazy_property
    def SWout(self):
        """Instantaneous outgoing shortwave radiation in watts per square meter"""
        return ptjpl.SWout(self.SWin, self.albedo).rename(['SWout']).set(self._properties)

    @lazy_property
    def SWnet(self):
        """Instantaneous outgoing shortwave radiation in watts per square meter"""
        return ptjpl.SWnet(self.SWin, self.SWout).rename(['SWnet']).set(self._properties)

    # @lazy_property
    # def LWin(self):
    #     """Instantaneous incoming longwave radiation in watts per square meter"""
    #     return ptjpl.LWin(self.Ta_K, self.Ea_Pa) \
    #         .rename(['LWin']).set(self._properties)

    @lazy_property
    def LWin(self):
        """Incoming longwave radiation in watts per square meter

        Returns
        -------
        ee.Image

        Raises
        ------
        ValueError
            If `self.LWin_source` is not supported.

        """
        if utils.is_number(self.LWin_source):
            LWin = ee.Image.constant(float(self.LWin_source))
        elif isinstance(self.LWin_source, ee.computedobject.ComputedObject):
            LWin = ee.Image(self.LWin_source)
        elif self.LWin_source.upper() == 'NLDAS':
            LWin = self.nldas_interpolate('longwave_radiation', self._date)
        else:
            raise ValueError(f'Unsupported LWin source: {self.LWin_source}\n')

        return LWin.rename(['LWin']).resample(DOWNSAMPLE_METHOD)

    @lazy_property
    def LWout(self):
        """Instantaneous outgoing longwave radiation in watts per square meter"""
        return ptjpl.LWout(self.LST, self.emissivity).rename(['LWout']).set(self._properties)

    @lazy_property
    def LWnet(self):
        """Instantaneous net longwave radiation in watts per square meter"""
        return ptjpl.LWnet(self.LWin, self.LWout).rename(['LWnet']).set(self._properties)

    @lazy_property
    def Rn(self):
        """Instantaneous net radiation"""
        return ptjpl.Rn(self.SWin, self.SWout, self.LWin, self.LWout) \
            .rename(['Rn']).set(self._properties)

    @lazy_property
    def Rnd(self):
        """Daily average net radiation"""
        return ptjpl.Rnd(self.Rn, self.hour_solar, self.sunrise_hour, self.daylight_hours) \
            .rename(['Rnd']).set(self._properties)

    @lazy_property
    def Rns(self):
        """Net radiation to the soil in watts per square meter"""
        return ptjpl.Rns(self.Rn, self.LAI, self.water_mask).rename(['Rns']).set(self._properties)

    @lazy_property
    def Rnc(self):
        """Net radiation to the canopy in watts per square meter"""
        return ptjpl.Rnc(self.Rn, self.Rns, self.water_mask).rename(['Rnc']).set(self._properties)

    @lazy_property
    def U(self):
        """Windspeed"""
        if utils.is_number(self.windspeed_source):
            windspeed_img = ee.Image.constant(float(self.windspeed_source))
        elif isinstance(self.windspeed_source, ee.computedobject.ComputedObject):
            windspeed_img = self.windspeed_source
        # elif self.windspeed_source.upper() == 'CFSV2':
        #     # It would be more correct to compute the magnitude for each image,
        #     #   then compute the average.
        #     # Do we need daily, 6hr, or interpolated instantaneous data?
        #     windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H') \
        #         .select([
        #         'u-component_of_wind_height_above_ground',
        #         'v-component_of_wind_height_above_ground']) \
        #         .filterDate(self._date, self._date.advance(1, 'day'))
        #     windspeed_img = windspeed_coll.mean() \
        #         .expression('sqrt(b(0) ** 2 + b(1) ** 2)')
        elif self.windspeed_source.upper() == 'NLDAS':
            wind_u = self.nldas_interpolate('wind_u', self._date)
            wind_v = self.nldas_interpolate('wind_v', self._date)
            windspeed_img = wind_u.expression(
                "sqrt(wind_u ** 2 + wind_v ** 2)",
                {"wind_u": wind_u, "wind_v": wind_v}
            )
        else:
            raise ValueError(f'Invalid windspeed_source: {self.windspeed_source}\n')

        windspeed_img = windspeed_img.resample(DOWNSAMPLE_METHOD)
        windspeed_img = self.LST.multiply(0).add(windspeed_img)

        return windspeed_img.rename(['U']).set(self._properties)

    @lazy_property
    def SVP_kPa(self):
        """Saturation vapor pressure in kilopascal"""
        return ptjpl.SVP_kPa(self.Ta_C).rename(['SVP_kPa']).set(self._properties)

    @lazy_property
    def VPD_kPa(self):
        """Vapor pressure deficit in kilopascal"""
        return ptjpl.VPD_kPa(self.Ea_kPa, self.SVP_kPa).rename(['VPD_kPa']).set(self._properties)

    @lazy_property
    def RH(self):
        """Relative humidity"""
        return ptjpl.RH(self.Ea_kPa, self.SVP_kPa).rename(['RH']).set(self._properties)

    @lazy_property
    def Td_C(self):
        """Relative humidity"""
        return ptjpl.Td(self.Ta_C, self.RH).rename(['Td']).set(self._properties)

    @lazy_property
    def W(self):
        """Instantaneous water heat flux"""
        W = ptjpl.W(
            self.ST_C,
            self.Td_C,
            self.U,
            self.Rn,
            self.SWnet,
            self.water_mask
        )

        return W.rename(['W']).set(self._properties)

    @lazy_property
    def G(self):
        """Instantaneous soil heat flux from Rn and fractional vegetation cover"""
        G = ptjpl.G(self.Rn, self.fIPAR, self.Rns, self.W, self.water_mask)
        # G = G.where(self.water_mask, self.W)
        # G = G.updateMask(self.water_mask.Not())

        return G.rename(['G']).set(self._properties)

    @lazy_property
    def A(self):
        """Instantaneous available energy from the difference of net radiation and ground heat flux"""
        return self.Rn.subtract(self.G).rename(['A']).set(self._properties)

    @lazy_property
    def delta(self):
        """Slope of saturation vapor pressure to air temperature"""
        return ptjpl.delta(self.Ta_C).rename(['delta']).set(self._properties)

    @lazy_property
    def fwet(self):
        """Relative surface wetness"""
        return ptjpl.fwet(self.RH).rename(['fwet']).set(self._properties)

    @lazy_property
    def fg(self):
        """Green canopy fraction"""
        return ptjpl.fg(self.fAPAR, self.fIPAR).rename(['fg']).set(self._properties)

    @lazy_property
    def fM(self):
        """Plant moisture constraint"""
        return ptjpl.fM(self.fAPAR, self.fAPARmax).rename(['fM']).set(self._properties)

    @lazy_property
    def fSM(self):
        """Soil moisture constraint"""
        return ptjpl.fSM(self.RH, self.VPD_kPa).rename(['fSM']).set(self._properties)

    @lazy_property
    def fT(self):
        """Plant temperature constraint"""
        return ptjpl.fT(self.Ta_C, self.Topt).rename(['fT']).set(self._properties)

    @lazy_property
    def epsilon(self):
        """Plant temperature constraint"""
        return ptjpl.epsilon(self.delta).rename(['epsilon']).set(self._properties)

    @lazy_property
    def LEc(self):
        """Instantaneous canopy transpiration in watts per square meter"""
        return ptjpl.LEc(self.fwet, self.fg, self.fT, self.fM, self.epsilon, self.Rnc) \
            .rename(['LEc']).set(self._properties)

    @lazy_property
    def LEi(self):
        """Instantaneous interception evaporation in watts per square meter"""
        return ptjpl.LEi(self.fwet, self.epsilon, self.Rnc) \
            .rename(['LEi']).set(self._properties)

    @lazy_property
    def LEs(self):
        """Instantaneous soil evaporation in watts per square meter"""
        return ptjpl.LEs(self.fwet, self.fSM, self.epsilon, self.Rns, self.G) \
            .rename(['LEs']).set(self._properties)

    @lazy_property
    def LE(self):
        """Instantaneous evapotranspiration in watts per square meter"""
        return ptjpl.LE(self.LEc, self.LEi, self.LEs, self.PET, self.water_mask) \
            .rename(['LE']).set(self._properties)

    @lazy_property
    def EF(self):
        """Evaporative fraction"""
        return ptjpl.EF(self.LE, self.Rn, self.G).rename(['EF']).set(self._properties)

    @lazy_property
    def LEd(self):
        """Daily latent heat flux"""
        return ptjpl.LEd(self.EF, self.Rnd).rename(['LEd']).set(self._properties)

    @lazy_property
    def PET(self):
        """Potential instantaneous latent heat flux in watts per square meter"""
        return ptjpl.PET(self.epsilon, self.Rn, self.G).rename(['PET']).set(self._properties)

    @lazy_property
    def ET(self):
        """Actual evapotranspiration in mm/day"""
        et_img = ptjpl.ET(self.LEd, self.daylight_hours)

        if self.crop_pm_adjust_flag:
            et_img = et_img.multiply(self.crop_pm_adjust)

        return et_img

        # return ptjpl.ET(self.LEd, self.daylight_hours) \
        #     .rename(['ET']).set(self._properties)

    @lazy_property
    def ESI(self):
        """Ratio of instantaneous latent heat flux to potential latent heat flux"""
        return ptjpl.ESI(self.LE, self.PET).rename(['ESI']).set(self._properties)

    @lazy_property
    def ea(self):
        """Actual vapor pressure [Pa]

        Returns
        -------
        ee.Image

        Raises
        ------
        ValueError
            If `self.ea_source` is not supported.

        """
        if utils.is_number(self.ea_source):
            ea_img = ee.Image.constant(float(self.ea_source))
        elif isinstance(self.ea_source, ee.computedobject.ComputedObject):
            ea_img = ee.Image(self.ea_source)
        elif self.ea_source.upper() == 'NLDAS':
            sph_img = self.nldas_interpolate('specific_humidity', self._date)

            # We could get pressure from NLDAS or compute from elevation?
            # NLDAS pressure is in Pa
            pair_img = self.nldas_interpolate('pressure', self._date)

            # # Elevation could be made into a class property
            # elev_img = ee.Image('USGS/SRTMGL1_003')
            # pair_img = elev_img.multiply(-0.0065).add(293.0).divide(293.0) \
            #     .pow(5.26).multiply(101.3).multiply(1000)

            # Compute Ea from sph (and pair) and convert kPa to Pa
            ea_img = sph_img.multiply(0.378).add(0.622).pow(-1) \
                .multiply(sph_img).multiply(pair_img)
        else:
            raise ValueError(f'Unsupported ea_source: {self.ea_source}\n')

        return ea_img.rename(['ea']).resample(DOWNSAMPLE_METHOD)

    @lazy_property
    def Ea_Pa(self):
        """Actual vapor pressure (Pascal)"""
        return self.ea.rename(['Ea_Pa']).set(self._properties)

    @lazy_property
    def Ea_kPa(self):
        """Actual vapor pressure (kilopascal)"""
        return self.Ea_Pa.divide(1000.0).rename(['Ea_kPa']).set(self._properties)

    @lazy_property
    def rs(self):
        """Incoming shortwave radiation in watts per square meter

        Returns
        -------
        ee.Image

        Raises
        ------
        ValueError
            If `self.rs_source` is not supported.

        """
        if utils.is_number(self.rs_source):
            rs_img = ee.Image.constant(float(self.rs_source))
        elif isinstance(self.rs_source, ee.computedobject.ComputedObject):
            rs_img = ee.Image(self.rs_source)
        elif self.rs_source.upper() == 'NLDAS':
            rs_img = self.nldas_interpolate('shortwave_radiation', self._date)
        else:
            raise ValueError(f'Unsupported rs_source: {self.rs_source}\n')

        return rs_img.rename(['rs']).resample(DOWNSAMPLE_METHOD)

    @lazy_property
    def ta(self):
        """Air temperature [K]

        Returns
        -------
        ee.Image

        Raises
        ------
        ValueError
            If `self.ta_source` is not supported.

        """
        if utils.is_number(self.ta_source):
            ta_img = ee.Image.constant(float(self.ta_source))
        elif isinstance(self.ta_source, ee.computedobject.ComputedObject):
            ta_img = ee.Image(self.ta_source)
        elif self.ta_source.upper() == 'NLDAS':
            ta_img = self.nldas_interpolate('temperature', self._date).add(273.15)
        else:
            raise ValueError(f'Unsupported ta_source: {self.ta_source}\n')

        return ta_img.rename(['Ta']).resample(DOWNSAMPLE_METHOD)

    @staticmethod
    def nldas_interpolate(band_name, interp_date, interp_flag=True):
        """
        Parameters
        ----------
        band_name : str
        interp_date : ee.Date
        interp_flag : bool

        Returns
        -------
        ee.Image

        Notes
        -----
        Interpolate rs hourly image at image time
        Hourly Rs is time average so time starts are 30 minutes early
        Move image time 30 minutes earlier to simplify filtering/interpolation

        """
        nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002').select([band_name])

        if interp_flag:
            # Interpolate
            # TODO: Check if NLDAS data is average over trailing hour, instantaneous, something else.
            #   The date may need to be shifted depending.
            # interp_date = interp_date.advance(-0.5, 'hour')
            prev_img = ee.Image(nldas_coll.filterDate(
                interp_date.advance(-1, 'hour'), interp_date).first())
            next_img = ee.Image(nldas_coll.filterDate(
                interp_date, interp_date.advance(1, 'hour')).first())
            prev_time = ee.Number(prev_img.get('system:time_start'))
            next_time = ee.Number(next_img.get('system:time_start'))
            interp_time = interp_date.millis().subtract(prev_time) \
                .divide(next_time.subtract(prev_time))
            output_img = next_img.subtract(prev_img).multiply(interp_time).add(prev_img)
        else:
            # Select the first NLDAS image after the image date
            output_coll = ee.ImageCollection(nldas_coll) \
                .filterDate(interp_date, interp_date.advance(1, 'hour'))
            # Select the first NLDAS image closest to the image date
            # output_coll = ee.ImageCollection(nldas_coll) \
            #     .filterDate(interp_date.advance(30, 'minute'),
            #                 interp_date.advance(30, 'minute'))
            output_img = ee.Image(output_coll.first())

        return output_img

    @lazy_property
    def Ta_K(self):
        """Near-surface air temperature in Kelvin"""
        return self.ta.rename(['Ta_K']).set(self._properties)

    @lazy_property
    def Ta_C(self):
        """Near-surface air temperature in Celsius"""
        Ta_C = self.Ta_K.subtract(273.15)
        Ta_C = Ta_C.resample(DOWNSAMPLE_METHOD)
        Ta_C = self.LST.multiply(0).add(Ta_C)

        return Ta_C.rename(['Ta_C']).set(self._properties)

    @lazy_property
    def Topt(self):
        """Optimum temperature [C]"""
        if utils.is_number(self.topt_source):
            # Interpret numbers as constant images
            topt_img = ee.Image.constant(float(self.topt_source))
        elif (type(self.topt_source) is str and
              (self.topt_source.lower().startswith('projects/') or
               self.topt_source.lower().startswith('users/'))):
            # Read the source as an image asset ID str
            topt_img = ee.Image(self.topt_source)
        # elif isinstance(self.topt_source, computedobject.ComputedObject):
        #     # Use the ee.Image() directly
        #     topt_img = ee.Image(self.topt_source)
        else:
            raise ValueError(f'unsupported topt_source: {self.topt_source}')

        topt_img = topt_img.resample(DOWNSAMPLE_METHOD)
        topt_img = self.Ta_C.multiply(0).add(topt_img)

        if self.floor_Topt:
            topt_img = topt_img.max(self.Ta_C)

        return topt_img.rename(['Topt']).set(self._properties)

    @lazy_property
    def fAPARmax(self):
        """Maximum fraction of absorbed photosynthetically activate radiation"""
        if utils.is_number(self.faparmax_source):
            # Interpret numbers as constant images
            fAPARmax_img = ee.Image.constant(float(self.faparmax_source))
        elif (type(self.faparmax_source) is str and
              (self.faparmax_source.lower().startswith('projects/') or
               self.faparmax_source.lower().startswith('users/'))):
            # Read the source as an image asset ID str
            fAPARmax_img = ee.Image(self.faparmax_source)
        # elif isinstance(self.faparmax_source, computedobject.ComputedObject):
        #     # Use the ee.Image() directly
        #     fAPARmax_img = ee.Image(self.fAPARmax_source)
        else:
            raise ValueError(f'unsupported fAPARmax_source: {self.faparmax_source}')

        return fAPARmax_img.rename(['fAPARmax']).resample(DOWNSAMPLE_METHOD)

    @classmethod
    def from_image_id(cls, image_id, **kwargs):
        """Constructs an PT-JPL Image instance from an image ID

        Parameters
        ----------
        image_id : str
            An earth engine image ID.
            (i.e. 'LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
        kwargs
            Keyword arguments to pass through to model init.

        Returns
        -------
        new instance of Image class

        """
        # DEADBEEF - Should the supported image collection IDs and helper
        # function mappings be set in a property or method of the Image class?
        collection_methods = {
            'LANDSAT/LT04/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LT05/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LE07/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC08/C02/T1_L2': 'from_landsat_c2_sr',
            'LANDSAT/LC09/C02/T1_L2': 'from_landsat_c2_sr',
        }

        try:
            method_name = collection_methods[image_id.rsplit('/', 1)[0]]
        except KeyError:
            raise ValueError(f'unsupported collection ID: {image_id}')
        except Exception as e:
            raise Exception(f'unhandled exception: {e}')

        method = getattr(Image, method_name)

        return method(ee.Image(image_id), **kwargs)

    @classmethod
    def from_landsat_c2_sr(cls, sr_image, cloudmask_args={}, **kwargs):
        """Returns a PT-JPL Image instance from a Landsat C02 level 2 (SR) image

        Parameters
        ----------
        sr_image : ee.Image, str
            A raw Landsat Collection 2 level 2 (SR) image or image ID.
        cloudmask_args : dict
            keyword arguments to pass through to cloud mask function
        kwargs : dict
            Keyword arguments to pass through to Image init function

        Returns
        -------
        Image

        Notes
        -----
        https://www.usgs.gov/faqs/how-do-i-use-a-scale-factor-landsat-level-2-science-products?qt-news_science_products=0#qt-news_science_products
        https://www.usgs.gov/core-science-systems/nli/landsat/landsat-collection-2-level-2-science-products

        """
        sr_image = ee.Image(sr_image)

        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(sr_image.get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_4': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_5': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_7': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL'],
            'LANDSAT_8': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
            'LANDSAT_9': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL'],
        })
        output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir', 'QA_PIXEL']

        prep_image = (
            sr_image.select(input_bands.get(spacecraft_id), output_bands)
            .multiply([0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.0000275, 0.00341802, 1])
            .add([-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 149.0, 0])
        )

        # Default the cloudmask flags to True if they were not
        # Eventually these will probably all default to True in openet.core
        if 'cirrus_flag' not in cloudmask_args.keys():
            cloudmask_args['cirrus_flag'] = True
        if 'dilate_flag' not in cloudmask_args.keys():
            cloudmask_args['dilate_flag'] = True
        if 'shadow_flag' not in cloudmask_args.keys():
            cloudmask_args['shadow_flag'] = True
        if 'snow_flag' not in cloudmask_args.keys():
            cloudmask_args['snow_flag'] = True
        if 'cloud_score_flag' not in cloudmask_args.keys():
            cloudmask_args['cloud_score_flag'] = False
        if 'cloud_score_pct' not in cloudmask_args.keys():
            cloudmask_args['cloud_score_pct'] = 100
        if 'filter_flag' not in cloudmask_args.keys():
            cloudmask_args['filter_flag'] = False
        # QA_RADSAT band will need to be added above if applying saturated masking
        # if 'saturated_flag' not in cloudmask_args.keys():
        #     cloudmask_args['saturated_flag'] = False

        cloud_mask = openet.core.common.landsat_c2_sr_cloud_mask(sr_image, **cloudmask_args)

        # Check if passing c2_lst_correct or soil_emis_coll_id arguments
        if "c2_lst_correct" in kwargs.keys():
            assert isinstance(kwargs['c2_lst_correct'], bool), "selection type must be a boolean"
            # Remove from kwargs since it is not a valid argument for Image init
            c2_lst_correct = kwargs.pop('c2_lst_correct')
        else:
            c2_lst_correct = cls._C2_LST_CORRECT

        if c2_lst_correct:
            lst = openet.core.common.landsat_c2_sr_lst_correct(sr_image, landsat.ndvi(prep_image))
        else:
            lst = prep_image.select(['tir'])

        # Build the input image
        input_image = ee.Image([
            landsat.albedo_metric(prep_image),
            landsat.emissivity_metric(prep_image),
            lst.rename(['lst']),
            # CGM - Don't compute LST since it is being provided
            # prep_image.select(['tir'], ['lst']),
            # landsat.lst(prep_image),
            landsat.ndvi(prep_image),
            landsat.ndwi(prep_image),
            landsat.mndwi(prep_image),
            landsat.wri(prep_image),
            # cloud_mask.Not(),
        ])

        # Apply the cloud mask and add properties
        input_image = (
            input_image.updateMask(cloud_mask)
            .set({'system:index': sr_image.get('system:index'),
                  'system:time_start': sr_image.get('system:time_start'),
                  'system:id': sr_image.get('system:id'),
            })
        )

        # Instantiate the class
        return cls(input_image, **kwargs)

    @lazy_property
    def crop_type(self):
        """Crop type

        Parameters
        ----------
        crop_type_source : int, str
            CDL image collection ID: 'USDA/NASS/CDL'
                Collection will be filtered to a single year that is closest
                to the Image year.
            CDL image ID for a specific year: 'USDA/NASS/CDL/2018'
            OpenET crop type image collection ID:
                'projects/openet/crop_type/v2021a'
                Collection will be mosaiced to a single image.
            Integer (will be converted to an EE constant image)
        _year : ee.Number
            Year is needed if crop_type_source is the CDL collection.

        Returns
        -------
        ee.Image

        Raises
        ------
        ValueError for unsupported crop_type_sources

        Notes
        -----
        Copied from _crop_type() in the OpenET SIMS model.py

        """
        properties = ee.Dictionary()

        if utils.is_number(self.crop_type_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.et_reference_source)
            crop_type_img = ee.Image.constant(self.crop_type_source).rename('crop_type')
            # properties = properties.set('id', 'constant')
        elif (type(self.crop_type_source) is str and
              self.crop_type_source.upper() == 'USDA/NASS/CDL'):
            # Use the CDL image closest to the image date
            year_min = ee.Number(2008)
            year_max = ee.Date(
                ee.ImageCollection('USDA/NASS/CDL')
                .limit(1, 'system:index', False).first()
                .get('system:time_start')
            ).get('year')
            # year_max = ee.Number.parse(ee.ImageCollection('USDA/NASS/CDL')\
            #     .limit(1, 'system:index', False).first().get('system:index'))
            cdl_year = ee.Number(self._year).min(year_max).max(year_min)
            cdl_coll = (
                ee.ImageCollection('USDA/NASS/CDL')
                .filterDate(ee.Date.fromYMD(cdl_year, 1, 1),
                            ee.Date.fromYMD(cdl_year.add(1), 1, 1))
                .select(['cropland'])
            )
            crop_type_img = ee.Image(cdl_coll.first())
            properties = properties.set('id', crop_type_img.get('system:id'))
        elif (type(self.crop_type_source) is str and
                self.crop_type_source.upper().startswith('USDA/NASS/CDL')):
            crop_type_img = ee.Image(self.crop_type_source).select(['cropland'])
            properties = properties.set('id', crop_type_img.get('system:id'))
        elif (type(self.crop_type_source) is str and
              (('projects/openet/crop_type' in self.crop_type_source.lower()) or
               ('projects/openet/assets/crop_type' in self.crop_type_source.lower()))):
            # Use the crop_type image closest to the image date
            crop_coll = ee.ImageCollection(self.crop_type_source)
            cdl_year = (
                ee.Number(self._year)
                .min(ee.Date(crop_coll.aggregate_max('system:time_start')).get('year'))
                .max(ee.Date(crop_coll.aggregate_min('system:time_start')).get('year'))
            )
            crop_type_coll = (
                crop_coll
                .filterDate(ee.Date.fromYMD(cdl_year, 1, 1),
                            ee.Date.fromYMD(cdl_year.add(1), 1, 1))
            )
            crop_type_img = crop_type_coll.mosaic()
            properties = properties.set('id', crop_type_coll.get('system:id'))
        # TODO: Add support for generic ee.Image and ee.ImageCollection sources
        # elif isinstance(self.crop_type_source, computedobject.ComputedObject):
        else:
            raise ValueError(f'unsupported crop_type_source: {self.crop_type_source}')

        return crop_type_img.rename(['crop_type']).set(properties)

    @lazy_property
    def crop_mask(self):
        # CGM - This is the crop list we are using to set when SIMS is included
        #   in the ensemble.
        # It is identical to the remap in ensemble_interpolate_asset_export.py
        #   except that wetlands (87, 190, 195) are being included.
        # TODO: This remap should probably be set from an external file
        return (
            self.crop_type.remap(
                list(range(0, 255)),
                [0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                 1, 1, 1, 0, 1, 1, 1, 1, 1, 0,
                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 # 100
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 # 200
                 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1])
            .rename('crop_mask').set(self._properties)
        )

    @lazy_property
    def crop_pm_adjust(self):
        """Crop Penman-Monteith Adjustment"""
        if utils.is_number(self.crop_pm_adjust_source):
            # Interpret numbers as constant images
            # CGM - Should we use the ee_types here instead?
            #   i.e. ee.ee_types.isNumber(self.crop_pm_adjust_source)
            crop_pm_adjust_img = ee.Image.constant(self.crop_pm_adjust_source)
        elif (type(self.crop_pm_adjust_source) is str and
              (self.crop_pm_adjust_source.startswith('projects/') or
               self.crop_pm_adjust_source.startswith('users/'))):
            # TODO: Copy all properties from source image to output
            # Assume a string source is an image ID
            crop_pm_adjust_img = ee.Image(self.crop_pm_adjust_source)
            if self.crop_pm_adjust_band:
                crop_pm_adjust_img = crop_pm_adjust_img.select([self.crop_pm_adjust_band])
            else:
                crop_pm_adjust_img = crop_pm_adjust_img.select([0])
            crop_pm_adjust_img = crop_pm_adjust_img.resample(DOWNSAMPLE_METHOD)
            # if DOWNSAMPLE_METHOD in ['bilinear', 'bicubic']:
            #     crop_pm_adjust_img = crop_pm_adjust_img.resample(DOWNSAMPLE_METHOD)
        else:
            raise ValueError(f'unsupported crop_pm_adjust_source: '
                             f'{self.crop_pm_adjust_source}')

        return (
            self.NDVI.multiply(0).add(1)
            .where(self.crop_mask, crop_pm_adjust_img)
            .rename(['crop_pm_adjust']).set(self._properties)
        )
