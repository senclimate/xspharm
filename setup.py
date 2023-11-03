#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy', 'xarray', 'pyspharm>=1.0.9']

test_requirements = [ ]

setup(
    author="Sen Zhao",
    author_email='zhaos2016@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="xarray interface to spherical harmonic transform",
    entry_points={
        'console_scripts': [
            'xspharm=xspharm.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='xspharm',
    name='xspharm',
    packages=find_packages(include=['xspharm', 'xspharm.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/senclimate/xspharm',
    version='0.1.0',
    zip_safe=False,
)
