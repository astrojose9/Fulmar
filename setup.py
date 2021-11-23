# -*- coding: utf-8 -*-
from os import path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Pull fulmar version from single source of truth file
try:  # Python 2
    execfile(path.join("fulmar", 'version.py'))
except:  # Python 3
    exec(open(path.join("fulmar", 'version.py')).read())


# If Python3: Add "README.md" to setup.
# Useful for PyPI (pip install fulmar-astro). Useless for users using Python2
try:

    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = ' '
except NameError:
    long_description = ' '


setup(name='fulmar-astro',
      version=FULMAR_VERSION,
      description='Follow-Up Lightcurves Multitool Assisting Radial velocities',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/astrojose9/fulmar',
      author='José Rodrigues',
      author_email='jose.rodrigues@astro.up.pt',
      license='MIT',
      classifiers=['Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Development Status :: 5 - Production/Stable',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python'],
      packages=['fulmar'],
      include_package_data=True,
      package_data={'': ['*.csv', '*.json']},
      install_requires=[
          'arviz',
          'astropy==4.3.1',
          "celerite2",
          "corner",
          'exoplanet',
          'lightkurve>=2',
          'matplotlib',
          'numpy<1.21,>=1.17',
          'pymc3',
          'pymc3-ext',
          'transitleastsquares>=1.0.31'
      ]
      )
