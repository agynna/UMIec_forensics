#!/usr/bin/env python3
from setuptools import setup

install_requires = ["fdstools"]

setup(name='umierrorcorrect',
      description='UMIerrorcorrect Forensics',
      #long_description = open('README.md').read(),
      #url='http://github.com/',
      author='Froste Svensson',
      author_email='froste.svensson@gmail.com',
      license='mit',
      package_data={'umierrorcorrect': ['README.md']
                   },
      include_package_data=True,
      install_requires=install_requires,
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['umierrorcorrect_forensics/run_tssv.py'],
      zip_safe=False)
