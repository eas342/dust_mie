#from distutils.core import setup
from setuptools import setup#, find_packages

with open("README.md", "r") as fh:
    ## skip the HTML, which doesn't work on PyPI
    long_description = "".join(fh.readlines()[0:])

setup(
    name='mie_dust',
    version='0.1.dev1',
    author='Everett Schlawin',
    packages=['mie_dust'],
    url="https://github.com/eas342/mie_dust",
    description="A package to evaluate dust coefficients for a variety of compositions",
    include_package_data=True,
    install_requires=[
        "numpy>=1.15",
        "scipy>=1.1.0",
        "astropy>=2.0",
        "joblib",
        "miepython",
    ],
    long_description=long_description,
    long_description_content_type='text/markdown'
)
