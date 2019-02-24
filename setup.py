#!/usr/bin/env python

from distutils.core import setup

setup(name='SUPPA',
      version='2.3',
      description='SUPPA2: Fast quantification of differential splicing',
      url='https://github.com/comprna/SUPPA',
      scripts=['suppa.py'],
      author="comprna",
      packages=['suppa'],
      install_requires=['scipy', 'numpy', 'pandas', 'statsmodels', 'sklearn'],
      python_requires='>=3',
)
