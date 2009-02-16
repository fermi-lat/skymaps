# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-01-16-03
import os
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

progEnv.Tool('skymapsLib')
test_healpix = progEnv.Program('test_skymaps', listFiles(['src/test/*.cxx']))

libEnv.Tool('skymapsLib', depsOnly = 1)
libEnv.Append(CPPPATH = ['#/healpix/','#/healpix/src'])
skymapsLib = libEnv.SharedLibrary('skymaps', listFiles(['src/*.cxx']))

swigEnv.Tool('skymapsLib')
swigEnv.Append(CPPPATH = ['#/healpix/','#/healpix/src'])
swigEnv.Replace(SHLIBPREFIX = '_')
#swigEnv.Replace(SHLIBSUFFIX = '.pyd')
swigEnv.Append(RPATH = swigEnv['LIBDIR'])
pySkymapsLib = swigEnv.SharedLibrary('skymaps','python/swig_setup.i')

progEnv.Tool('registerObjects', 
        package = 'skymaps', 
	includes = listFiles(['skymaps/*.h']),
	libraries = [skymapsLib],
	python = ['python/skymaps.py',pySkymapsLib],		
	)

