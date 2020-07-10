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

setup(name='taurex_petitrad',
      author="Ahmed Faris Al-Refaie",
      author_email="ahmed.al-refaie.12@ucl.ac.uk",
      license="BSD",
      description='petitRADTRANS plugin for TauREx-3 ',
      packages=packages,
      
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      install_requires=install_requires)