from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
        ext_modules=cythonize("sonics.pyx"),
        include_dirs=[numpy.get_include()],
        name="sonics",
        version="0.2.0",
        description=description
)

description = """
SONiCS performs dense forward simulations of the PCR of Short Tandem Repeats
from capture experiments, calculates the likelihood of generating the provided
read support (reads per allele) out of the final PCR pool and determines
the most probable genotype based on the log likelihood distributions
from all such simulations.
"""
