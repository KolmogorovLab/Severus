import os
import sys

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from severus.__version__ import __version__

with open('requirements.txt') as f:
    install_requires = f.read().splitlines()

setup(name='severus',
      version=__version__,
      description='A tool for somatic structural variant calling using long reads',
      url='https://github.com/KolmogorovLab/Severus',
      author='Ayse Keskus',
      author_email = 'aysegokce.keskus@nih.gov',
      license='BSD-3-Clause',
      packages=['severus'],
      package_data={'severus': ['vntrs/*']},
      install_requires=install_requires,
      python_requires='>=3.8',
      entry_points={'console_scripts': ['severus = severus.main:main']},
      )
