#from distutils.core import setup
from setuptools import setup#, find_packages

with open("README.md", "r") as fh:
    ## skip the HTML, which doesn't work on PyPI
    long_description = "".join(fh.readlines()[0:])

setup(
    name='mie_dust',
    version='0.1dev1',
    author='Everett Schlawin, Kayli Glidic',
    packages=['tshirt','tshirt.pipeline',
              'tshirt.pipeline.instrument_specific'],
    url="https://github.com/eas342/tshirt",
    description="A package to analyze time series data, especially for exoplanets",
    include_package_data=True,
    install_requires=[
        "numpy>=1.15",
        "scipy>=1.1.0",
        "astropy>=2.0",
        "miepython",
    ],
    long_description=long_description,
    long_description_content_type='text/markdown'
)
