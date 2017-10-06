# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-02-17-00
import os
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
package= 'skymaps'
libname = package+'Lib'
testname = 'test_'+package

libEnv.Tool('addLinkDeps', package=package, toBuild='shared')

progEnv.Tool(libname)
testapp = progEnv.Program(testname, listFiles(['src/test/*.cxx']))

if libEnv.get('CONTAINERNAME','') != 'ScienceTools_User':
    lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))
else:
    lib = libEnv.SharedLibrary(package, listFiles(['src/B*.cxx',
                                                   'src/Compos*.cxx',
                                                   'src/Convolution.cxx',
                                                   'src/DiffuseFunction.cxx',
                                                   'src/E*.cxx','src/G*.cxx',
                                                   'src/H*.cxx', 'src/I*.cxx',
                                                   'src/L*.cxx', 'src/P*.cxx',
                                                   'src/S*.cxx','src/W*.cxx']))

# SWIG
swigEnv = baseEnv.Clone()
#Get the path for the numpy includes
numpy_path = os.path.join(baseEnv['pythonSitePath'], 'numpy', 'core', 'include')
swigEnv.AppendUnique(CPPPATH=numpy_path)
swigEnv.Tool(libname)
pyLib = swigEnv.SwigLibrary('_'+package,'src/swig_setup.i')

progEnv.Tool('registerTargets',
             package = package,
             includes = listFiles([package+'/*.h']),
             libraryCxts = [[lib, libEnv]],
             ##  testAppCxts = [[testapp, progEnv]],
             data = ['data/LivetimeCubeTemplate', 'data/aeff_P6_v1_diff_front.fits'],
             swigLibraryCxts = [[pyLib, swigEnv]],
             python = ['src/skymaps.py']
             )

