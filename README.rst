===============
OpenET - PT-JPL
===============

|version| |build|

**WARNING: This code is in development, is being provided without support, and is subject to change at any time without notification**

This repository provides an Earth Engine Python API based implementation of the PT-JPL model for computing evapotranspiration (ET).

Through ecophysiological constraint functions, PT-JPL retrieves actual ET by reducing potential ET starting with the Priestley-Taylor equation (PriestleyTaylor1972_). A series of ecophysiological scalar functions, based on atmospheric vapor pressure deficit, relative humidity, and vegetation indices simultaneously reduce potential ET to actual ET, and partition total ET into three sources for canopy transpiration, soil evaporation, and interception evaporation (Fisher2008_). PT-JPL is run globally and continuously in space and time with no need for calibration or site-specific parameters.

Model Design
============

The primary component of the PT-JPL model is the Image() class.  The Image class can be used to compute a single ET image from a single input image.  The Image class should generally be instantiated from an Earth Engine Landsat image using the collection specific methods listed below.  ET image collections can be built by computing ET in a function that is mapped over a collection of input images.  Please see the `Example Notebooks`_ for more details.

Input Collections
=================

PT-JPL can currently be computed for Landsat Collection 2 Level 2 (SR/ST) images  images from the following Earth Engine image collections:

 * LANDSAT/LT05/C02/T1_L2
 * LANDSAT/LE07/C02/T1_L2
 * LANDSAT/LC08/C02/T1_L2
 * LANDSAT/LC09/C02/T1_L2

Landsat Collection 2 SR/ST Input Image
--------------------------------------

To instantiate the class for a Landsat Collection 2 SR/ST image, use the Image.from_landsat_c2_sr method.

The input Landsat image must have the following bands and properties:

=================  ======================================
SPACECRAFT_ID      Band Names
=================  ======================================
LANDSAT_5          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, QA_PIXEL
LANDSAT_7          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B7, ST_B6, QA_PIXEL
LANDSAT_8          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, QA_PIXEL
LANDSAT_9          SR_B1, SR_B2, SR_B3, SR_B4, SR_B5, SR_B6, SR_B7, ST_B10, QA_PIXEL
=================  ======================================

=================  =============================================
Property           Description
=================  =============================================
system:index       - Landsat Scene ID
                   - Must be in the Earth Engine format (e.g. LC08_044033_20170716)
                   - Used to lookup the scene specific c-factor
system:time_start  Image datetime in milliseconds since 1970
SPACECRAFT_ID      - Used to determine which Landsat type
                   - Must be: LANDSAT_5, LANDSAT_7, LANDSAT_8, or LANDSAT_9
=================  =============================================

Model Output
------------

The primary output of the PT-JPL model is the actual ET (ETa) in millimeters.

Example
-------

.. code-block:: python

    import openet.ptjpl as ptjpl

    landsat_img = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044033_20170716')
    et_actual = ptjpl.Image.from_landsat_c2_sr(landsat_img).et

Example Notebooks
=================


Meteorology
===========

The `GEE NLDAS hourly image collection <https://developers.google.com/earth-engine/datasets/catalog/NASA_NLDAS_FORA0125_H002>`__ is the default and only currently supported meteorology source for the hourly air temperature, vapor pressure, windspeed, and incoming short and longwave solar radiation.

Ancillary Datasets
==================

Optimal temperature (Topt)
--------------------------


Maximum fraction of absorbed photosynthetically active radiation (fAPARmax)
---------------------------------------------------------------------------


Installation
============

The OpenET PT-JPL python module can be installed via pip:

.. code-block:: console

    pip install openet-ptjpl

Dependencies
============

 * `earthengine-api <https://github.com/google/earthengine-api>`__
 * `openet-core <https://github.com/Open-ET/openet-core-beta>`__

OpenET Namespace Package
========================

Each OpenET model is stored in the "openet" folder (namespace).  The model can then be imported as a "dot" submodule of the main openet module.

.. code-block:: console

    import openet.ptjpl as model

Development and Testing
=======================

Please see the `CONTRIBUTING.rst <CONTRIBUTING.rst>`__.

References
==========

.. _references:

.. [Fisher2008]
 | Fisher, J., K. Tu, and D. Baldocchi (2008). Global estimates of the land-atmosphere water flux based on monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites, *Remote Sensing of Environment*, 112(3), 901-919.
 | `https://doi.org/10.1016/j.rse.2007.06.025 <https://doi.org/10.1016/j.rse.2007.06.025>`__
.. [PriestleyTaylor1972]
 | Priestley, C. and R. Taylor (1972). On the assessment of surface heat flux and evaporation using large scale parameters. *Monthly Weather Review*, 100, 81–92.
 | `https://doi.org/10.1175/1520-0493%281972%29100%3C0081%3AOTAOSH%3E2.3.CO%3B2 <https://doi.org/10.1175/1520-0493%281972%29100%3C0081%3AOTAOSH%3E2.3.CO%3B2>`__

.. |build| image:: https://github.com/Open-ET/openet-ptjpl/actions/workflows/tests.yml/badge.svg
   :alt: Build status
   :target: https://github.com/Open-ET/openet-ptjpl
.. |version| image:: https://badge.fury.io/py/openet-ptjpl.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/openet-ptjpl
