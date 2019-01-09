#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 3 14:14:06 2019
Setup script for ProTiler
@author: Wei He
@email: whe3@mdanderson.org
"""

from setuptools import setup
import re

def main():
    version = re.search('^__version__\s*=\s*"(.*)"',open('bin/protiler').read(), re.M).group(1)
    setup(name='ProTiler',
          version=version,
          author='Wei He, Liang Zhang, Oscar Villarreal',
          author_email='whe3@mdanderson.org',
          description="Call HS regions from CRISPR tiling screen data and predict HS region from common protein features",
          include_package_data=True,
          packages=['ProTiler'],
          package_dir={'ProTiler':'ProTiler'},
          package_data={'ProTiler':['StaticFiles/*']},
          url='http://mageck.sourceforge.net',
          scripts=['bin/protiler'],
          install_requires=[
                  'numpy',
                  'pandas',
                  'matplotlib',
                  'sklearn',
                  'argparse',
                  'seaborn'],
          zip_safe = True
        )
if __name__ == '__main__':
    main()
