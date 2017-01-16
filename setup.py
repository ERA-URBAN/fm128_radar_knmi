import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "fm128-radar-knmi",
    version = "0.0.1",
    author = "Ronald van Haren",
    author_email = "r.vanharen@esciencecenter.nl",
    description = ("A python library to download and convert KNMI "
                   "radar files to WRFDA fm128_radar ascii"),
    license = "Apache 2.0",
    keywords = "WRFDA netCDF WRF radar",
    url = "https://github.com/ERA-URBAN/fm128_radar_knmi",
    packages=['fm128_radar_knmi'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved ::Apache Software License",
    ],
)
