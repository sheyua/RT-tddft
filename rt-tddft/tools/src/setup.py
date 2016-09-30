from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(Extension(
	'tdpost',
	sources=[
		'tdpost.pyx',
		'C_Rho.cpp'#, 'C_Vks.cpp'
	],
	include_dirs=[
		numpy.get_include()
	],
	extra_compile_args=['-O3'],
	language='c++',
)))
