# -*- coding: utf-8 -*-
#
#  setup.py
#  xman
#
#  Created by Alan D Snow, 2017.
#  BSD 3-Clause

from setuptools import setup, find_packages

requires = [
    'descartes',
    'gazar',
    'matplotlib',
    'numpy',
    'shapely',
]

setup(name='xman',
      version='0.0.1',
      description="Cross section (X-Section) extractor and Manning's "
                  "roughness (N) solver.",
      # long_description='',
      author='Alan D. Snow',
      author_email='alansnow21@gmail.com',
      url='https://github.com/snowman2/xman',
      license='BSD 3-Clause',
      keywords='cross section, mannning, hydraulics, hydrology',
      packages=find_packages(),
      classifiers=[
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      install_requires=requires,
      extras_require={
          'tests': [
              'coveralls',
              'flake8',
              'pytest',
              'pytest-cov',
              'pylint',
          ],
          'docs': [
              'mock',
              'sphinx',
              'sphinxcontrib-napoleon',
          ]
      },
      )
