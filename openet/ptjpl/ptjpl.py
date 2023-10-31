"""
OpenET implementation of PT-JPL from LST
Gregory Halverson, Jet Propulsion Laboratory
"""
import ee

from openet.ptjpl.daily_integration import daily_integration, calculate_vapor

__author__ = "Gregory Halverson"

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


def SAVI(NDVI):
    """Return soil adjusted vegetation index (SAVI) image"""
    SAVI = NDVI.multiply(0.45).add(0.132)

    return SAVI


def fAPAR(SAVI):
    """Return fraction of absorbed photosynthetically active radiation"""
    fAPAR = SAVI.multiply(1.3632).subtract(0.048).clamp(0, 1)

    return fAPAR


def fIPAR(NDVI):
    """Return fraction of intercepted photosynthetically active radiation"""
    fIPAR = NDVI.subtract(0.05).clamp(0, 1)

    return fIPAR


def LAI(NDVI, fIPAR):
    """Return leaf area index (LAI) image"""

    # calculate leaf area index
    LAI = NDVI.expression('-log(1 - fIPAR) * (1 / KPAR)', {'fIPAR': fIPAR, 'KPAR': KPAR})

    return LAI


def SWout(SWin, albedo):
    """
    Parameters
    ----------
    :param SWin: ee.Image incoming shortwave radiation in watts per square meter
    :param albedo: ee.Image land-surface albedo

    Returns
    -------
    :return: ee.Image Instantaneous outgoing shortwave radiation in watts per square meter
    """
    SWout = SWin.multiply(albedo)

    return SWout


def LWin(Ta_K, Ea_Pa):
    """
    Parameters
    ----------
    :param Ta_K: ee.Image near-surface air temperature in Kelvin
    :param Ea_Pa: ee.Image near-surface vapor pressure in Pascal

    Returns
    -------
    :return: ee.Image Instantaneous incoming longwave radiation in watts per square meter
    """
    # calculate atmospheric emissivity
    eta1 = Ea_Pa.divide(Ta_K).multiply(0.465)
    atmospheric_emissivity = eta1.expression(
        '1 - (1 + eta1) * exp(-sqrt(1.2 + 3 * eta1))',
        {'eta1': eta1}
    )

    # calculate incoming longwave radiation
    LWin = Ta_K.expression(
        'atmospheric_emissivity * STEFAN_BOLTZMAN_CONSTANT * Ta_K ** 4',
        {
            'Ta_K': Ta_K,
            'atmospheric_emissivity': atmospheric_emissivity,
            'STEFAN_BOLTZMAN_CONSTANT': STEFAN_BOLTZMAN_CONSTANT
        }
    )

    return LWin


def LWout(LST, emissivity):
    """
    Parameters
    ----------
    :param LST: ee.Image land-surface temperature
    :param emissivity: ee.Image land-surface emissivity

    Returns
    -------
    :return: ee.Image Instantaneous outgoing longwave radiation in watts per square meter
    """
    # calculate outgoing longwave radiation
    LWout = LST.expression(
        'emissivity * STEFAN_BOLTZMAN_CONSTANT * LST ** 4',
        {
            'LST': LST,
            'emissivity': emissivity,
            'STEFAN_BOLTZMAN_CONSTANT': STEFAN_BOLTZMAN_CONSTANT
        }
    )

    return LWout


def SWnet(SWin, SWout):
    """
    Parameters
    ----------
    :param SWin: ee.Image incoming shortwave radiation in watts per square meter
    :param SWout: ee.Image outgoing shortwave radiation in watts per square meter

    Returns
    -------
    :return: ee.Image Net shortwave in watts per square meter
    """
    # calculate net shortwave radiation
    SWnet = SWin.subtract(SWout)

    return SWnet


def LWnet(LWin, LWout):
    """
    Parameters
    ----------
    :param LWin: ee.Image incoming longwave radiation in watts per square meter
    :param LWout: ee.Image outgoing longwave radiation in watts per square meter

    Returns
    -------
    :return: ee.Image New longwave in watts per square meter
    """
    # calculate net longwave radiation
    LWnet = LWin.subtract(LWout)

    return LWnet


def Rn(SWin, SWout, LWin, LWout):
    """
    Parameters
    ----------
    :param SWin: ee.Image incoming shortwave radiation in watts per square meter
    :param SWout: ee.Image outgoing shortwave radiation in watts per square meter
    :param LWin: ee.Image incoming longwave radiation in watts per square meter
    :param LWout: ee.Image outgoing longwave radiation in watts per square meter

    Returns
    -------
    :return: ee.Image Instantaneous net radiation in watts per square meter
    """
    # calculate net shortwave radiation
    SWnet = SWin.subtract(SWout)

    # calculate net longwave radiation
    LWnet = LWin.subtract(LWout)

    # calculate instantaneous net radiation
    Rn = SWnet.add(LWnet)

    return Rn


def Rnd(Rn, hour_of_day, sunrise_hour, daylight_hours):
    return daily_integration(Rn, hour_of_day, sunrise_hour, daylight_hours)


def Rns(Rn, LAI, water_mask):
    """Net radiation of the soil from leaf area index
    Parameters
    ----------
    :param Rn: ee.Image
    :param LAI: ee.Image
    :param water_mask: ee.Image

    Returns
    -------
    ee.Image
    """
    Rns = Rn.expression('Rn * exp(-KRN * LAI)', {'Rn': Rn, 'KRN': KRN, 'LAI': LAI})

    Rns = Rns.updateMask(water_mask.Not())

    return Rns


def Rnc(Rn, Rns, water_mask):
    """Net radiation of the soil from leaf area index
    Parameters
    ----------
    :param Rn: ee.Image
    :param Rns: ee.Image

    Returns
    -------
    ee.Image
    """
    # calculate net radiation of the canopy from net radiation of the soil
    Rnc = Rn.subtract(Rns)
    Rnc = Rnc.updateMask(water_mask.Not())

    return Rnc


def W(WST, Td_C, U, SWnet, Rn, water_mask):
    """
    Water heat flux method from AquaSEBS
    http://www.mdpi.com/2072-4292/8/7/583
    """
    W_MAX_PROPORTION = 1

    Tn = WST.subtract(Td_C).multiply(0.5)

    eta = WST.expression(
        "0.35 + 0.015 * water_surface_temperature + 0.0012 * (Tn ** 2)",
        {
            "water_surface_temperature": WST,
            "Tn": Tn
        }
    )

    S = U.multiply(3.3)

    beta = WST.expression(
        "4.5 + 0.05 * water_surface_temperature + (eta + 0.47) * S",
        {
            "water_surface_temperature": WST,
            "eta": eta,
            "S": S
        }
    )

    Te = SWnet.divide(beta).add(Td_C)
    W = Te.subtract(WST).multiply(beta)
    Wmax = Rn.multiply(W_MAX_PROPORTION)
    W = W.min(Wmax)
    W = W.updateMask(water_mask)

    return W


def G(Rn, fIPAR, Rns, W, water_mask=None):
    """Instantaneous soil heat flux from Rn and fractional vegetation cover

    Parameters
    ----------
    :param Rn: ee.Image instantaneous net radiation
    :param fIPAR: ee.Image fraction of intercepted photosynthetically active radiation
    :param Rns: ee.Image net radiation to the soil in watts per square meter
    :param W: ee.Image water heat flux in watts per square meter

    Returns
    -------
    :return: ee.Image soil heat flux in watts per square meter

    """
    G = Rn.expression(
        'Rn * (0.05 + (1 - fIPAR) * 0.265)',
        {
            'Rn': Rn,
            'fIPAR': fIPAR.clamp(0, 1)
        }
    )

    # floor soil heat flux at zero
    G = G.max(0)

    # constrain soil heat flux to maximum allowed proportion of net radiation to the soil
    G = G.min(Rns.multiply(G_MAX_PROPORTION))

    # fill in water heat flux
    # G = G.where(water_mask, W)
    # G = G.where(water_mask.Not(), W)
    G = Rn.multiply(0).add(G.unmask(0).add(W.unmask(0)))

    return G


def SVP_kPa(Ta_C):
    # calculate saturation vapor pressure in kPa from air temperature in celsius
    SVP = Ta_C.multiply(17.27).divide(Ta_C.add(237.7)).exp().multiply(0.611)

    # floor saturation vapor pressure at 1
    SVP = SVP.max(1)

    return SVP


def VPD_kPa(Ea_kPa, SVP_kPa):
    # calculate vapor pressure deficit from water vapor pressure
    VPD = SVP_kPa.subtract(Ea_kPa)

    # CGM - I think the .gte(0) is correct since you are telling updateMask()
    #   which pixels are valid
    # TODO: Check if the updateMask() can can be skipped (and then change tests from images to numbers)
    # lower bound of vapor pressure deficit is zero, negative values replaced with nodata
    VPD = VPD.updateMask(VPD.gte(0))
    # VPD = np.where(VPD < 0, np.nan, VPD)

    return VPD


def RH(Ea_kPa, SVP_kPa):
    """
    Calculates saturation vapor pressure, vapor pressure deficit, and relative
    humidity from air temperature and water vapor pressure.
    :param Ea_kPa: water vapor pressure in kilopascals
    :param SVP_kPa: saturation vapor pressure in kilopascals
    :return: relative humidity as a proportion between 0 and 1
    """
    # calculate relative humidity from water vapor pressure and saturation vapor pressure
    RH = Ea_kPa.divide(SVP_kPa)

    # upper bound of relative humidity is one, results higher than one are capped at one
    RH = RH.min(1)

    return RH


def Td(Ta_C, RH):
    Td = RH.expression("Ta_C - ((100 - RH * 100) / 5.0)", {"Ta_C": Ta_C, "RH": RH})

    return Td


def fwet(RH):
    """
    Calculates relative surface wetness from relative humidity.
    :param RH: relative humidity as a proportion between 0 and 1
    :return: relative surface wetness as a proportion between 0 and 1
    """
    fwet = RH.pow(4)

    return fwet


# def delta(Ta_C):
#     """
#     Calculates slope of saturation vapor pressure to air temperature.
#     :param Ta_C: air temperature in celsius
#     :return: slope of the vapor pressure curve in kilopascals/kelvin
#     """
#     delta = Ta_C.multiply(17.27).divide(Ta_C.add(237.7)).exp() \
#         .multiply(0.6108 * 4098).divide(Ta_C.add(237.3).pow(2.0))
#
#     return delta

def delta(Ta_C):
    """
    Calculates slope of saturation vapor pressure to air temperature.
    :param Ta_C: air temperature in celsius
    :return: slope of the vapor pressure curve in kilopascals/kelvin
    """
    delta = Ta_C.expression(
        '4098 * (0.6108 * exp(17.27 * Ta_C / (237.7 + Ta_C))) / (Ta_C + 237.3) ** 2',
        {'Ta_C': Ta_C}
    )

    return delta


def fg(fAPAR, fIPAR):
    """green canopy fraction"""
    fg = fAPAR.divide(fIPAR).clamp(0, 1)

    return fg


def fM(fAPAR, fAPARmax):
    """plant moisture constraint"""
    fM = fAPAR.divide(fAPARmax).clamp(0, 1)

    return fM


def fSM(RH, VPD_kPa):
    """plant moisture constraint"""
    fSM = RH.pow(VPD_kPa.divide(BETA)).clamp(0, 1)

    return fSM


def fT(Ta_C, Topt):
    """plant temperature constraint"""
    fT = Ta_C.expression('exp(-(((Ta_C - Topt) / Topt) ** 2))', {'Ta_C': Ta_C, 'Topt': Topt})

    return fT


def epsilon(delta):
    epsilon = delta.expression(
        'delta / (delta + PSYCHROMETRIC_GAMMA)',
        {'delta': delta, 'PSYCHROMETRIC_GAMMA': PSYCHROMETRIC_GAMMA}
    )

    return epsilon


def LEc(fwet, fg, fT, fM, epsilon, Rnc):
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

    return LEc


def LEs(fwet, fSM, epsilon, Rns, G):
    # calculate soil evaporation (LEs) from relative surface wetness,
    # soil moisture constraint, priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), net radiation of the soil,
    # and soil heat flux
    LEs = Rns.expression(
        '(fwet + fSM * (1 - fwet)) * PT_ALPHA * epsilon * (Rns - G)',
        {
            'fwet': fwet,
            'fSM': fSM,
            'PT_ALPHA': PT_ALPHA, 'epsilon': epsilon,
            'Rns': Rns, 'G': G
        })

    # floor soil evaporation at zero
    LEs = LEs.max(0)

    return LEs


def LEi(fwet, epsilon, Rnc):
    # INTERCEPTION EVAPORATION

    # calculate interception evaporation (LEi) from relative surface wetness
    # and net radiation of the canopy
    LEi = fwet.multiply(PT_ALPHA).multiply(epsilon).multiply(Rnc)

    # floor interception evaporation at zero
    LEi = LEi.max(0)

    return LEi


def LE(LEc, LEi, LEs, PET, water_mask):
    return PET.where(water_mask.Not(), LEc.add(LEi).add(LEs))
    # LE = LEc.add(LEi).add(LEs)
    # LE = LE.where(water_mask, PET)
    #
    # return LE


def EF(LE, Rn, G):
    return LE.expression('LE / (Rn - G)', {'LE': LE, 'Rn': Rn, 'G': G})


def LEd(EF, Rnd):
    return EF.multiply(Rnd)


def PET(epsilon, Rn, G):
    PET = Rn.expression(
        'PT_ALPHA * epsilon * (Rn - G)',
        {'PT_ALPHA': PT_ALPHA, 'epsilon': epsilon, 'Rn': Rn, 'G': G}
    )

    return PET


def ET(LEd, daylight_hours):
    return calculate_vapor(LEd, daylight_hours)


def ESI(LE, PET):
    return LE.divide(PET)
