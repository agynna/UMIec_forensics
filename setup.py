#!/usr/bin/env python3
from setuptools import setup
from setuptools import find_packages

pack = find_packages() 

install_requires = ["https://github.com/tobbeost/umierrorcorrect/archive/v0.19.tar.gz","fdstools"]

exec(open('umierrorcorrect/version.py').read())

setup(name='umierrorcorrect',
      version=__version__,
      description='UMI error correct',
      long_description = open('README.md').read(),
      url='http://github.com/',
      author='Froste Svensson',
      author_email='froste.svensson@gmail.com',
      packages=pack,
      license='mit',
      package_data={'umierrorcorrect': ['README.md']
                   },
      include_package_data=True,
      install_requires=install_requires,
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['umierrorcorrect_forensics/run_tssv.py'],
      zip_safe=False)
