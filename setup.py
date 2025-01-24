# SPDX-License-Identifier: GPL-3.0-or-later
from setuptools import setup

setup(
    name='vagusNerve',
    version='0.1',    
    description='Semi-analytic model of the evoked compound action potential',
    url='https://github.com/tharayil/vagusNerve',
    author='Blue Brain Project, EPFL',
    license='GPL-3.0',
    packages=['vagusNerve'],
    install_requires=[
    'scikit-learn',
    'scipy',
    'numpy',
    'pandas',
    'notebook',
    'ipython',
    'matplotlib',
    'ipympl',
    'quantities',
    'openpyxl',
    'pytest',
    'pytest-cov'
    ]
   )
