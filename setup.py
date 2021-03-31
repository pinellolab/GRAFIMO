"""

Build GRAFIMO

"""

from setuptools import setup, find_packages, Extension
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as sd
from distutils.command.build_ext import build_ext as be

import sys

try:  # sphinx could not be available
    from sphinx.setup_command import BuildDoc
except:
    errmsg = "\n\nPlease install \"sphinx\" before installing GRAFIMO.\n"
    sys.stderr.write(errmsg)
    sys.exit(3)


if sys.version_info[:2] < (3,6): # python 3.7 is required
    """
    Check Python version
    It must be >= 3.6
    """
    
    sys.stderr.write("Pyhton >= 3.6 is required to run GRAFIMO\n")
    sys.exit(1)

# read README file
encoding_arg={'encoding': 'utf-8'} if sys.version_info[0] >= 3 else dict()
readmefile = 'README.md'
with open(readmefile, **encoding_arg) as infile:
    long_description = infile.read()

# Cython code build
CYTHON_V_REQUIRED='0.28'  # minimum Cython version required

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

name = "GRAFIMO"
version = '1.1.4'

# definition of setup()
setup(
      name='grafimo',
      version=version,
      author='Manuel Tognon',
      author_email='manu.tognon@gmail.com',
      url="https://github.com/pinellolab/GRAFIMO",
      description='GRAph-based Finding of Indivividual Motif Occurrences',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT', 
      cmdclass={'build_ext': BuildExt, 'sdist': SDist, 'build_sphinx': BuildDoc},
      command_options={
          'build_sphinx': {
              'project': ('setup.py', name),
              'version': ('setup.py', version),
              'source_dir': ('setup.py', 'docs') 

          }
      },
      ext_modules=extensions,
      packages=find_packages('src'),
      package_dir={'':'src'},
      entry_points={'console_scripts':['grafimo = grafimo.__main__:main']},
      install_requires=[
              'pandas>=0.24.2',
              'numpy>=1.16.4',
              'statsmodels>=0.11.0',
              'numba>=0.47',
              'sphinx>=3.5.2',
              'colorama'
              ],
      extras_require={
          'dev': ['Cython']
      },
      python_requires='>=3.6',
      classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    )
    
