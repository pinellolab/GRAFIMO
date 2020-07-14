"""

Build GRAFIMO

"""

from setuptools import setup, find_packages, Extension
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as sd
from distutils.command.build_ext import build_ext as be
import sys


if sys.version_info[:2] < (3,6): # python 3.7 is required
    """
    Check Python version
    It must be >= 3.6
    """
    
    sys.stderr.write("Pyhton >= 3.6 is required to run GRAFIMO\n")
    sys.exit(1)

# read README.md    
encoding_arg={'encoding':'utf-8'} if sys.version_info[0] >= 3 else dict()
readme='README.md'
with open(readme, **encoding_arg) as infile:
    long_description=infile.read()

# Cython code build
CYTHON_V_REQUIRED='0.28' # minimum Cython version required

def check_cython():
    """
        Check Cython version
    """

    try:
        from Cython import __version__ as cyv

    except ImportError:
        sys.stderr.write("Cython not found on your machine. Install Cython >= "+
                            str(CYTHON_V_REQUIRED))
        sys.exit(1)

    if LooseVersion(cyv) < LooseVersion(CYTHON_V_REQUIRED):
        sys.stderr.write("Found Cython v" + str(cyv)+" . Cython v"+
                            str(CYTHON_V_REQUIRED)+" is required")
        sys.exit(1)

extensions=[
    Extension('motif_processing', sources=['src/grafimo/motif_processing.pyx']),
]

class BuildExt(be):

    def run(self):

        check_cython()
        from Cython.Build import cythonize
        self.extensions=cythonize(self.extensions)

        super().run()

class SDist(sd):

    def run(self):

        check_cython()
        from Cython.Build import cythonize
        cythonize(extensions)
        super().run()

# definition of setup()
setup(
      name='grafimo',
      version='1.0.1',
      author='Manuel Tognon',
      author_email='manu.tognon@gmail.com',
      url='https://github.com/InfOmics/GRAFIMO',
      description='GRAph-based Find Indivividual Motif Occurrences',
      long_description=long_description,
      license='MIT', 
      cmdclass={'build_ext': BuildExt, 'sdist': SDist},
      ext_modules=extensions,
      packages=find_packages('src'),
      package_dir={'':'src'},
      entry_points={'console_scripts':['grafimo = grafimo.__main__:main']},
      install_requires=[
              'pandas~=0.24.2',
              'numpy~=1.16.4',
              'statsmodels~=0.10.0',
              'numba~=0.47',
              ],
      extras_require={
          'dev': ['Cython']
      },
      python_requires='>=3.6',
      classifiers=[
        "Development Status :: 1 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOSMojave",
        "Operating System :: Linux",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    )
    
