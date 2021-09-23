#from distutils.core import setup
# from setuptools import setup
#import numpy

try:
    from Cython.Build import cythonize
except:
    import pip
    pip.main(['install', 'cython'])
    from Cython.Build import cythonize

try:
    import numpy
except:
    import pip
    pip.main(['install', 'numpy'])
    import numpy

from numpy.distutils.core import Extension, setup
#from numpy.distutils.core import setup

description = """
SONiCS performs dense forward simulations of the PCR of Short Tandem Repeats
from capture experiments, calculates the likelihood of generating the provided
read support (reads per allele) out of the final PCR pool and determines
the most probable genotype based on the log likelihood distributions
from all such simulations.
"""

# TODO: put on pip, then this is necessary
# extmodules = [
#    Extension("pymc_extracted",
#              sources=['pymc_extracted.f'])#,
#]

setup(
        ext_modules=cythonize('sonics.pyx', language_level=3),
        include_paths=[numpy.get_include()],
        name="sonics",
        author="Kasia Kedzierska",
        author_email="kasia@well.ox.ac.uk",
        url="git@github.com:kzkedzierska/sonics.git",
        license="LICENCE",
        #setup_requires=['numpy', 'cython'],
        install_requires=['cython', 'numpy', 'scipy', 'pandas'],
        version="0.2.0",
        description=description
)
