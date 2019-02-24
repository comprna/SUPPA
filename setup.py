#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='SUPPA',
    packages=find_packages(),
    scripts=['suppa2'],
    version='2.3',
    description='A tool to study splicing across multiple conditions at high speed and accuracy.',
    author='GP Alamancos',
    author_email='eduardo.eyras@upf.edu',
    license='MIT',
    url='https://github.com/comprna/SUPPA',
    download_url='https://github.com/comprna/SUPPA/archive/v2.3.tar.gz',
    keywords=['alternative', 'splicing', 'analysis', 'transcriptomics'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.5'],
    install_requires=['scipy>=0.15.1',
                      'numpy>=1.11.0',
                      'pandas>=0.18.0',
                      'statsmodels>=0.6.1',
                      'scikit-learn>=0.16.1'],
    python_requires='>=3',
)
