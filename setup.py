import codecs
import os
import re

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Single sourcing code from here:
#   https://packaging.python.org/guides/single-sourcing-package-version/
here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError('Unable to find version string.')

model_name = 'PTJPL'
version = find_version('openet', model_name.lower(), '__init__.py')

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=f'openet-{model_name.lower()}',
    version=version,
    description=f'Earth Engine based {model_name} Model',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    license='Apache',
    author='Gregory Halverson',
    author_email='Gregory.H.Halverson@jpl.nasa.gov',
    url=f'https://github.com/Open-ET/openet-{model_name.lower()}',
    download_url=f'https://github.com/Open-ET/openet-{model_name.lower()}/'
                 f'archive/v{version}.tar.gz',
    install_requires=['earthengine-api', 'openet-core', 'python-dateutil'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=[f'openet.{model_name.lower()}'],
    keywords=f'{model_name} OpenET Evapotranspiration Earth Engine',
    classifiers = [
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    zip_safe=False,
)
