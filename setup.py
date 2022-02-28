#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = '0'

setup(
    name='counting_molecules',
    version=version,
    author='Jack Hellerstedt and Tyler Hennen',
    author_email='hellerstedt.jack@gmail.com',
    description='Library to count molecules in STM images',
    long_description=long_description,
    url='https://github.com/thennen/counting-molecules',
    project_urls = {
        "Bug Tracker": "https://github.com/thennen/counting-molecules/issues"
    },
    license='MIT',
    packages=find_packages(),
    install_requires=['matplotlib', 'matplotlib_scalebar', 'nanonispy', 'mahotas', 'scipy', 'scikit-image', 'sklearn'],
)