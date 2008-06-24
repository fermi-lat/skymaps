# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-01-02-01
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

progEnv.Tool('skymapsLib')
test_healpix = progEnv.Program('test_skymaps', listFiles(['src/test/*.cxx']))

libEnv.Tool('skymapsLib', depsOnly = 1)
progEnv.Tool('registerObjects', 
        package = 'skymaps', 
	includes = listFiles(['skymaps/*.h']),
	libraries = [libEnv.SharedLibrary('skymaps', listFiles(['src/*.cxx']))],		
	)

