#!/usr/bin/env python3
from setuptools import setup

install_requires = ["fdstools>=2",
                    "matplotlib",
                    "pysam",
                    "umierrorcorrect>=0.19"]

plotting_requires = ["matplotlib",
                     "pandas",
                     "seaborn"]

setup(name='umierrorcorrect_forensics',
      description='UMIerrorcorrect Forensics',
      #long_description = open('README.md').read(),
      #url='http://github.com/',
      author='Froste Svensson & Arvid H Gynnå',
      author_email='froste.svensson@gmail.com',
      license='mit',
      package_data={'umierrorcorrect_forensics': ['README.md']
                   },
      include_package_data=True,
      install_requires=install_requires,
      extras_requires={'diagnostics' : plotting_requires}
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['umierrorcorrect_forensics/run_umierrorcorrect_forensics.py'],
      zip_safe=False)
