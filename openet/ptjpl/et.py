import ee

from openet.ptjpl import daily_integration
from openet.ptjpl import daylight_hours as daylight
from openet.ptjpl import meteorology
from openet.ptjpl import ptjpl
from openet.ptjpl import vegetation


# SI units watts per square meter per kelvin to the fourth
STEFAN_BOLTZMAN_CONSTANT = 5.67036713e-8

# Priestley-Taylor coefficient alpha
PT_ALPHA = 1.26
BETA = 1.0

# maximum proportion of soil heat flux out of net radiation to the soil
G_MAX_PROPORTION = 0.35

# psychrometric constant gamma in pascals over kelvin
# same as value for ventilated (Asmann type) psychrometers, with an air movement of some 5 m/s
# http://www.fao.org/docrep/x0490e/x0490e07.htm
PSYCHROMETRIC_GAMMA = 0.0662  # Pa/K

KRN = 0.6
KPAR = 0.5


# TODO: These could be moved into a Model class?
def et(
        LST,
        emissivity,
        NDVI,
        albedo,
        Ta_K,
        Ea_Pa,
        SWin,
        # LWin=0,  # CGM - Computed in the function, commenting out for now
        Topt,
        fAPARmax,
        datetime,
        cloud_mask=None,
        latitude=None,
        # hour,
        # sunrise_hour,
        # daylight_hours,
):
    """Evapotranspiration [mm?]

    Parameters
    ----------
    LST : ee.Image
        Land surface temperature [K].
    emissivity : ee.Image
        Emissivity.
    NDVI : ee.Image
        Normalized difference vegetation index.
    albedo : ee.Image
        Albedo.
    Ta_K : ee.Image
        Air temperature [K].
    Ea_Pa : ee.Image
        Vapor pressure [Pa].
    SWin : ee.Image
        Incoming shortwave radiation [?].
    # CGM - This is computed below, commenting out for now
    # LWin : ee.Image
    #     Incoming longwave radiation [?].
    Topt : ee.Image
        Optimum temperature [C].
    fAPARmax : ee.Image
        Maximum fraction of absorbed photosynthetically active radiation (PAR).
    cloud_mask : ee.Image
        Mask with clouds flagged as 1.
    datetime : ee.Date
        Image datetime.
    latitude : ee.Image
        Latitude [deg].  If not set will default to ee.Image.pixelLonLat().
    # CGM - Computed below based on datetime and latitude
    # hour :
    # sunrise_hour :
    # daylight_hours :

    Returns
    -------
    ee.Image

    References
    ----------


    """
    # PREPROCESS
    # CGM - Computing hour, sunrise_hour, and daylight_hours from datetime and
    # latitude instead of passing as parameters
    # TODO: Check if int() and double() are needed
    doy = ee.Number(datetime.getRelative('day', 'year')).int().add(1).double()
    if latitude is None:
        latitude = ee.Image.pixelLonLat().select(['latitude'])
    # CGM - Is hour integer, float, local, UTC, etc?
    hour = ee.Date(datetime).get('hour')
    #     .add(ee.Date(datetime).get('minute').divide(60))

    sha_deg = daylight.sha_deg_from_doy_lat(doy, latitude)
    sunrise_hour = daylight.sunrise_from_sha(sha_deg)
    daylight_hours = daylight.daylight_from_sha(sha_deg)

    # PT-JPL START

    # calculate outgoing shortwave from incoming shortwave and albedo
    SWout = SWin.multiply(albedo)

    # calculate atmospheric emissivity
    eta1 = Ea_Pa.divide(Ta_K).multiply(0.465)
    # eta1 = 0.465 * Ea_Pa / Ta_K

    atmospheric_emissivity = eta1.expression(
        '1 - (1 + eta1) * exp(-sqrt(1.2 + 3 * eta1))', {'eta1': eta1}
    )
    # atmospheric_emissivity = eta1.multiply(3).add(1.2).sqrt().multiply(-1).exp() \
    #     .multiply(eta1.add(1)).multiply(-1).add(1)
    # atmospheric_emissivity = (1 - (1 + eta1) * np.exp(-(1.2 + 3 * eta1) ** 0.5))

    # CGM - Should clouds be masked out in the input image?
    # consider cloud atmosphere as black-body
    atmospheric_emissivity = atmospheric_emissivity.where(cloud_mask, 1.0)
    # atmospheric_emissivity = np.where(cloud_mask, 1.0, atmospheric_emissivity)

    # calculate incoming longwave radiation
    LWin = Ta_K.expression(
        'atmospheric_emissivity * STEFAN_BOLTZMAN_CONSTANT * Ta_K ** 4',
        {
            'Ta_K': Ta_K,
            'atmospheric_emissivity': atmospheric_emissivity,
            'STEFAN_BOLTZMAN_CONSTANT': STEFAN_BOLTZMAN_CONSTANT
        }
    )
    # LWin = Ta_K.pow(4).multiply(atmospheric_emissivity).multiply(STEFAN_BOLTZMAN_CONSTANT)
    # LWin = atmospheric_emissivity * STEFAN_BOLTZMAN_CONSTANT * Ta_K ** 4

    # calculate outgoing longwave radiation
    LWout = LST.expression(
        'emissivity * STEFAN_BOLTZMAN_CONSTANT * LST ** 4',
        {
            'LST': LST,
            'emissivity': emissivity,
            'STEFAN_BOLTZMAN_CONSTANT': STEFAN_BOLTZMAN_CONSTANT
        }
    )
    # LWout = LST.multiply(4).multiply(emissivity).multiply(STEFAN_BOLTZMAN_CONSTANT)
    # LWout = emissivity * STEFAN_BOLTZMAN_CONSTANT * LST ** 4

    # calculate net shortwave radiation
    SWnet = SWin.subtract(SWout)

    # calculate net longwave radiation
    LWnet = LWin.subtract(LWout)

    # calculate instantaneous net radiation
    Rn = SWnet.add(LWnet)

    # integrate net radiation to daily value
    Rnd = daily_integration.daily_integration(Rn, hour, sunrise_hour, daylight_hours)

    # METEOROLOGICAL VARIABLES

    # convert temperatures from kelvin to celsius
    Ta_C = meteorology.kelvin_to_celsius(Ta_K)

    # scale water vapor pressure from Pa to kPa
    Ea_kPa = meteorology.pascal_to_kilopascal(Ea_Pa)

    # calculate meteorology
    # TODO - Split into separate functions
    SVP_kPa, VPD_kPa, RH = meteorology.meteorology(Ta_C, Ea_kPa)

    # calculate relative surface wetness from relative humidity
    fwet = meteorology.fwet_from_RH(RH)

    # calculate slope of saturation to vapor pressure curve Pa/K
    delta = meteorology.delta_from_Ta(Ta_C)

    # VEGETATION VALUES

    # calculate fAPAR from NDVI
    fAPAR = vegetation.fAPAR_from_ndvi(NDVI)

    # calculate fIPAR from NDVI
    fIPAR = vegetation.fIPAR_from_ndvi(NDVI)

    # calculate green canopy fraction (fg) from fAPAR and fIPAR,
    # constrained between zero and one
    fg = fAPAR.divide(fIPAR).clamp(0, 1)
    # fg = np.clip(fAPAR / fIPAR, 0, 1)

    # calculate plant moisture constraint (fM) from fraction of
    # photosynthetically active radiation, constrained between zero and one
    fM = fAPAR.divide(fAPARmax).clamp(0, 1)
    # fM = np.clip(fAPAR / fAPARmax, 0.0, 1.0)

    # calculate soil moisture constraint from mean relative humidity and vapor
    # pressure deficit, constrained between zero and one
    fSM = RH.pow(VPD_kPa.divide(BETA)).clamp(0, 1)
    # fSM = np.clip(RH ** (VPD_kPa / BETA), 0.0, 1.0)

    # calculate plant temperature constraint (fT) from optimal phenology
    fT = Ta_C.expression('exp(-(((Ta_C - Topt) / Topt) ** 2))', {'Ta_C': Ta_C, 'Topt': Topt})
    # fT = Ta_C.subtract(Topt).divide(Topt).pow(2).multiply(-1).exp()
    # fT = np.exp(-(((Ta_C - Topt) / Topt) ** 2))

    # calculate leaf area index
    LAI = NDVI.expression('-log(1 - fIPAR) * (1 / KPAR)', {'fIPAR': fIPAR, 'KPAR': KPAR})
    # LAI = fIPAR.multiply(-1).add(1).log().multiply(-1).divide(KPAR)
    # LAI = -np.log(1 - fIPAR) * (1 / KPAR)

    # calculate delta / (delta + gamma)
    epsilon = delta.expression(
        'delta / (delta + PSYCHROMETRIC_GAMMA)',
        {'delta': delta, 'PSYCHROMETRIC_GAMMA': PSYCHROMETRIC_GAMMA}
    )
    # epsilon = delta.add(PSYCHROMETRIC_GAMMA).pow(-1).multiply(delta)
    # epsilon = delta / (delta + PSYCHROMETRIC_GAMMA)

    # SOIL EVAPORATION

    # calculate net radiation of the soil from leaf area index
    Rns = ptjpl.Rns(Rn, LAI)

    # calculate instantaneous soil heat flux from Rn and fractional vegetation cover
    G = ptjpl.G(Rn, fIPAR, Rns)

    # calculate soil evaporation (LEs) from relative surface wetness,
    # soil moisture constraint, priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), net radiation of the soil,
    # and soil heat flux
    LEs = Rns.expression(
        '(fwet + fSM * (1 - fwet)) * PT_ALPHA * epsilon * (Rns - G)',
        {'fwet': fwet, 'fSM': fSM, 'PT_ALPHA': PT_ALPHA, 'epsilon': epsilon, 'Rns': Rns, 'G': G}
    )
    # LEs = Rns.subtract(G).multiply(PT_ALPHA).multiply(epsilon) \
    #     .multiply(fwet.multiply(-1).add(1).multiply(fsm).add(fwet))
    # LEs = (fwet + fSM * (1 - fwet)) * PT_ALPHA * epsilon * (Rns - G)

    # TODO: TODO: Convert this section if necessary
    # CGM - Why would LEs be missing or nodata?  Is this step necessary?
    # replace missing soil evaporation with zero
    # LEs = LEs.updateMask(0)
    # LEs = np.where(np.isnan(LEs), 0, LEs)

    # floor soil evaporation at zero
    LEs = LEs.max(0)

    # CANOPY TRANSPIRATION

    # calculate net radiation of the canopy from net radiation of the soil
    Rnc = Rn.subtract(Rns)

    # calculate canopy transpiration (LEc) from priestley taylor, relative
    # surface wetness, green canopy fraction, plant temperature constraint,
    # plant moisture constraint, epsilon = delta / (delta + gamma),
    # and net radiation of the canopy
    LEc = Rnc.expression(
        'PT_ALPHA * (1 - fwet) * fg * fT * fM * epsilon * Rnc',
        {
            'PT_ALPHA': PT_ALPHA,
            'fwet': fwet,
            'fg': fg,
            'fT': fT,
            'fM': fM,
            'epsilon': epsilon,
            'Rnc': Rnc
        }
    )
    # LEc =  fwet.multiply(-1).add(1).multiply(PT_ALPHA).multiply(epsilon) \
    #     .multiply(fg).multiply(fT).multiply(fM).multiply(Rnc)
    # LEc = PT_ALPHA * (1 - fwet) * fg * fT * fM * epsilon * Rnc

    # # TODO: Convert this section if necessary
    # # CGM - Why would LEs be missing or nodata?  Is this step necessary?
    # # replace missing canopy transpiration with zero
    # LEc = LEc.updateMask(0)
    # LEc = np.where(np.isnan(LEc), 0, LEc)

    # floor canopy transpiration at zero
    LEc = LEc.max(0)

    # INTERCEPTION EVAPORATION

    # calculate interception evaporation (LEi) from relative surface wetness
    # and net radiation of the canopy
    LEi = fwet.multiply(PT_ALPHA).multiply(epsilon).multiply(Rnc)

    # # TODO: Convert this section if necessary
    # # CGM - Why would LEs be missing or nodata?  Is this step necessary?
    # # replace missing interception evaporation with zero
    # LEi = LEi.updateMask(0)
    # LEi = np.where(np.isnan(LEi), 0, LEi)

    # floor interception evaporation at zero
    LEi = LEi.max(0)

    # COMBINED EVAPOTRANSPIRATION

    # combine soil evaporation (LEs), canopy transpiration (LEc), and
    # interception evaporation (LEi) into instantaneous evapotranspiration (LE)
    LE = LEs.add(LEc).add(LEi)

    # # TODO: TODO: Convert this section if necessary
    # # CGM - This might happen automatically depending on the calculations above
    # # mask LE with LST
    # LE = np.where(np.isnan(NDVI), np.nan, LE)

    # # TODO: Convert this section if necessary
    # # CGM - What is the purpose of this section?
    # #   Is the int16 cast important?
    # # mask latent heat flux to the soil
    # LEs = np.where(np.isnan(LE), np.nan, LEs)
    # LEs = np.where(np.isnan(LEs), 0, LEs)
    # # CGM - Is the int16 cast important?
    # LEs = np.int16(LEs / LE * 100)
    #
    # # mask latent heat flux to the canopy
    # LEc = np.where(np.isnan(LE), np.nan, LEc)
    # LEc = np.where(np.isnan(LEc), 0, LEc)
    # LEc = np.int16(LEc / LE * 100)
    #
    # # mask interception
    # LEi = np.where(np.isnan(LE), np.nan, LEi)
    # LEi = np.where(np.isnan(LEi), 0, LEi)
    # LEi = np.int16(LEi / LE * 100)

    # if LE is clipped to Rn, should partitions be scaled down?
    # wouldn't make more sense to clip LE to PET?
    # CGM - Can't use .clamp() with an image so use .max().min() combo
    LE = LE.max(0).min(Rn)
    # LE = np.clip(LE, 0, Rn)

    # DAILY EVAPOTRANSPIRATION

    # calculate evaporative fraction (EF) from evapotranspiration, net
    # radiation, and soil heat flux
    EF = LE.divide(Rn.subtract(G)).where(Rn.eq(G), 0)
    # EF = np.where(Rn - G == 0, 0, LE / (Rn - G))

    # calculate daily latent heat flux from evaporative fraction and daily net
    # radiation with minimum of zero
    LE_daily = EF.multiply(Rnd).max(0)
    # LE_daily = np.clip(EF * Rnd, 0, None)

    # calculate daytime daily evapotranspiration in kilograms equivalent to millimeters
    ET_daily_kg = daily_integration.calculate_vapor(LE_daily, daylight_hours)

    return ET_daily_kg
