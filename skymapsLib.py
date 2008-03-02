# @file skymapsLib.py
# @brief scons package dependencies for skymaps
#
#$Header$
def generate(env, **kw):
	env.Tool('addLibrary', library=['skymaps'], package = 'skymaps')
        depends = 'facilities tip'.split()
        for pack in depends: env.Tool(pack+'Lib')

	# why these?
	#env.Tool('addLibrary', library=env['clhepLibs'])
	#env.Tool('addLibrary', library=env['cfitsioLibs'])

def exists(env):
	return 1
