import os
from setuptools import setup
import sys

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def conf_path(name):
  if sys.prefix == '/usr':
    conf_path = os.path.join('/etc', name)
  else:
    conf_path = os.path.join(sys.prefix, 'etc', name)
  return conf_path

setup(
    name = "fm128_radar_knmi",
    version = "0.0.1",
    author = "Ronald van Haren",
    author_email = "r.vanharen@esciencecenter.nl",
    description = ("A python library to download and convert KNMI "
                   "radar files to WRFDA fm128_radar ascii"),
    license = "Apache 2.0",
    keywords = "WRFDA netCDF WRF radar",
    url = "https://github.com/ERA-URBAN/fm128_radar_knmi",
    packages=['fm128_radar_knmi'],
    data_files=[(os.path.join(conf_path('fm128_radar_knmi')),
                 ['fm128_radar_knmi/radar.config'])],
    scripts=['fm128_radar_knmi/scripts/fm128_radar_knmi'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved ::Apache Software License",
    ],
)
