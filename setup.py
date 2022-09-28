#!/usr/bin/env python3
from setuptools import setup
from setuptools.command.install import install

install_requires = ["fdstools>=2.0.1",
                    "umierrorcorrect==0.22",
                    "pysam>=0.19.1"]

setup(name='umierrorcorrect_forensics',
      description='UMIerrorcorrect Forensics',
      #long_description = open('README.md').read(),
      #url='http://github.com/',
      author='Froste Svensson & Arvid H Gynnå',
      author_email='froste.svensson@gmail.com',
      license='mit',
      package_data={'umierrorcorrect_forensics': ['README.md']},
      version="0.1",
      include_package_data=True,
      install_requires=install_requires,
      #dependency_links=["https://github.com/tobbeost/umierrorcorrect/archive/v0.19.tar.gz"],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages=["umierrorcorrect_forensics",
                "umierrorcorrect_forensics.tools"],
      #cmdclass={'install': CustomInstall},
      scripts=["umierrorcorrect_forensics/run_umierrorcorrect_forensics.py",
               "umierrorcorrect_forensics/run_flash.py",
               "umierrorcorrect_forensics/run_tssv.py",
               "umierrorcorrect_forensics/run_fdstools.py",
               "umierrorcorrect_forensics/fastq2sam.py",
               "umierrorcorrect_forensics/run_umifilter.py",
               "umierrorcorrect_forensics/convert_fastq2bam.py",
               "umierrorcorrect_forensics/convert_bam2fastq.py",
               "umierrorcorrect_forensics/tools/uncollapse_reads.py",
               "umierrorcorrect_forensics/tools/downsample.py",
               "umierrorcorrect_forensics/tools/barcode_diversity.py",
               "umierrorcorrect_forensics/tools/get_stutters.py"],
      zip_safe=False)
