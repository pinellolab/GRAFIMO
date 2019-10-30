"""

Build GRAFIMO

"""

from setuptools import setup, find_packages
import sys


if sys.version_info[:2] < (3,7): # python 3.7 is required
    """
    Check Python version
    It must be >= 3.7
    """
    
    sys.stderr.write("Pyhton >= 3.7 is required to run GRAFIMO\n")
    sys.exit(1)

# read README.md    
encoding_arg={'encoding':'utf-8'} if sys.version_info[0] >= 3 else dict()
readme='README.md'
with open(readme, **encoding_arg) as infile:
    long_description=infile.read()
    

# definition of setup()
setup(
      name='grafimo',
      version='0.6',
      author='Manuel Tognon',
      author_email='manu.tognon@gmail.com',
      url='https://github.com/InfOmics/GRAFIMO',
      description='GRAph-based Find Indivividual Motif Occurrences',
      long_description=long_description,
      license='MIT', 
      #include_package_data = True,
      #package_data={'grafimo': ['./*']}, 
      packages=find_packages('src'),
      package_dir={'':'src'},
      entry_points={'console_scripts':['grafimo = grafimo.__main__:main']},
      install_requires=[
              'pandas~=0.24.2',
              'numpy~=1.16.4',
              ],
      python_requires='>=3.7',
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
    


