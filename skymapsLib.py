# @file skymapsLib.py
# @brief scons package dependencies for skymaps
#
#$Header$
def generate(env, **kw):
	if not kw.get('depsOnly',0):
		env.Tool('addLibrary', library=['skymaps'])
        depends = 'facilities tip st_facilities healpix embed_python'.split()
        for pack in depends: env.Tool(pack+'Lib')

	# why these?
	#env.Tool('addLibrary', library=env['clhepLibs'])
	#env.Tool('addLibrary', library=env['cfitsioLibs'])

def exists(env):
	return 1
