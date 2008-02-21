#$Id$
import glob,os

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

if libEnv['PLATFORM'] == "win32":
	libEnv.AppendUnique(CPPFLAGS = "/wd4068")
	if libEnv['MSVS_VERSION'] == "8.0":
		libEnv.AppendUnique(CPPFLAGS = "/wd4812")

skymapsSharedLib = libEnv.SharedLibrary('skymaps', listFiles(['src/*.cxx']))

progEnv.Tool('skymapsLib')
test_healpix = progEnv.Program('test_skymaps', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'skymaps', 
	libraries = [skymapsSharedLib],		
	includes = listFiles(['skymaps/*.h']))

