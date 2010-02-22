# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-02-01-01
import os
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

if baseEnv['PLATFORM'] == "win32":
    libEnv.Tool('skymapsLib', depsOnly = 1)

progEnv.Tool('skymapsLib')
test_skymaps = progEnv.Program('test_skymaps', listFiles(['src/test/*.cxx']))

skymapsLib = libEnv.SharedLibrary('skymaps', listFiles(['src/*.cxx']))

swigEnv.Tool('skymapsLib')
pySkymapsLib = swigEnv.SwigLibrary('_skymaps','src/swig_setup.i')

progEnv.Tool('registerTargets',
             package = 'skymaps',
             includes = listFiles(['skymaps/*.h']),
             libraryCxts = [[skymapsLib, libEnv]],
             swigLibraryCxts = [[pySkymapsLib, swigEnv]],
             testAppCxts = [[test_skymaps, progEnv]],
             data = ['data/LivetimeCubeTemplate', 'data/aeff_P6_v1_diff_front.fits'],
             python = ['src/skymaps.py']
             )

