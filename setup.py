#!/usr/bin/env python
import setuptools
from setuptools import find_packages
from distutils.core import setup
from distutils.core import Extension
from distutils import log
import re, os

packages = find_packages(exclude=('tests', 'doc'))
provides = ['taurex_petitrad', ]

requires = []

install_requires = ['taurex', ]

entry_points = {'taurex.plugins': 'petitrad = taurex_petitrad'}

with open("README.md", "r") as fh:
    long_description = fh.read()

version = "0.1.0-alpha"

setup(name='taurex_petitrad',
      author="Ahmed Faris Al-Refaie",
      author_email="ahmed.al-refaie.12@ucl.ac.uk",
      license="BSD",
      description='petitRADTRANS plugin for TauREx-3 ',
      packages=packages,
      long_description=long_description,
      long_description_content_type='text/markdown',
      entry_points=entry_points,
      keywords=['exoplanet',
                'chemistry'
                'taurex',
                'plugin',
                'taurex3',
                'petitradtrans',
                'atmosphere',
                'atmospheric'],
      provides=provides,
      requires=requires,
      install_requires=install_requires)