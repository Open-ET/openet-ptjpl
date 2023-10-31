# try:
#     from importlib import metadata
# except ImportError:  # for Python<3.8
#     import importlib_metadata as metadata

from .collection import Collection
from .image import Image
from . import interpolate

MODEL_NAME = 'PTJPL'

# # __version__ = metadata.version(__package__ or __name__)
# __version__ = metadata.version(__package__.replace('.', '-') or __name__.replace('.', '-'))
# # __version__ = metadata.version('openet-ptjpl')
