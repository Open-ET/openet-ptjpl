def savi_from_ndvi(ndvi):
    """
    Linearly calculates Soil-Adjusted Vegetation Index from LST.
    :param ndvi: normalized difference vegetation index clipped between 0 and 1
    :return: soil-adjusted vegetation index
    """
    return ndvi.multiply(0.45).add(0.132)
    # return ndvi * 0.45 + 0.132


def fAPAR_from_savi(savi):
    """
    Linearly calculates fraction of absorbed photosynthetically active radiation from SAVI.
    :param savi: soil adjusted vegetation index
    :return: fraction of absorbed photosynthetically active radiation
    """
    # CGM - Double check that the correct operation is ".subtract(0.048)"
    return savi.multiply(1.3632).subtract(0.048).clamp(0, 1)
    # return np.clip(savi * 1.3632 + -0.048, 0, 1)


def fAPAR_from_ndvi(ndvi):
    """
    Calculates fraction of absorbed photosynthetically active radiation from LST.
    :param ndvi: LST clipped between 0 and 1
    :return:
    """
    return fAPAR_from_savi(savi_from_ndvi(ndvi))


def fIPAR_from_ndvi(ndvi):
    """
    Calculate fraction of intercepted photosynthetically active radiation from LST
    :param ndvi: LST clipped between 0 and 1
    :return: fraction of intercepted photosynthetically active radiation
    """
    # CGM - The clamp on +1 seems unecessary since NDVI can't be >
    return ndvi.subtract(0.05).clamp(0, 1)
    # return np.clip(ndvi - 0.05, 0, 1)
