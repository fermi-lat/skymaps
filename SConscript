# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-01-23-08
import os
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

progEnv.Tool('skymapsLib')
test_skymaps = progEnv.Program('test_skymaps', listFiles(['src/test/*.cxx']))

libEnv.Tool('skymapsLib', depsOnly = 1)
#libEnv.Append(CPPPATH = ['#/healpix/','#/healpix/src'])
skymapsLib = libEnv.SharedLibrary('skymaps', listFiles(['src/*.cxx']))


swigEnv.Tool('skymapsLib')
pySkymapsLib = swigEnv.SwigLibrary('_skymaps','python/swig_setup.i')

progEnv.Tool('registerTargets',
             package = 'skymaps',
             includes = listFiles(['skymaps/*.h']),
             libraryCxts = [[skymapsLib, libEnv]],
             swigLibraryCxts = [[pySkymapsLib, swigEnv]],
             testAppCxts = [[test_skymaps, progEnv]],
             data = ['data/LivetimeCubeTemplate'],
             python = ['python/skymaps.py']
             )

