# -*- python -*-
# @file SConscript
# @brief scons build specifications for skymaps
#
# $Header$
# Authors: T. Burnett <tburnett@u.washington.edu>, M.Roth <mar0@u.washington.edu>
# Version: skymaps-02-10-02
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

lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))
if baseEnv['PLATFORM']=='win32':
    # Add a post-build step to embed the manifest using mt.exe
    # The number at the end of the line indicates the file type (1: EXE; 2:DLL).
    libEnv.AddPostAction(lib, 'mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2')

# SWIG
swigEnv = baseEnv.Clone()
#Get the path for the numpy includes
import numpy
numpy_path = os.path.join(numpy.__path__[0],'core','include')
swigEnv.AppendUnique(CPPPATH=numpy_path)
swigEnv.Tool(libname)
pyLib = swigEnv.SwigLibrary('_'+package,'src/swig_setup.i')
if baseEnv['PLATFORM']=='win32':
    libEnv.AddPostAction(pyLib, 'mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2')

progEnv.Tool('registerTargets',
             package = package,
             includes = listFiles([package+'/*.h']),
             libraryCxts = [[lib, libEnv]],
             testAppCxts = [[testapp, progEnv]],
             data = ['data/LivetimeCubeTemplate', 'data/aeff_P6_v1_diff_front.fits'],
             swigLibraryCxts = [[pyLib, swigEnv]],
             python = ['src/skymaps.py']
             )

