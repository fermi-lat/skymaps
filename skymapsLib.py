# @file skymapsLib.py
# @brief scons package dependencies
#
#$Header$
def generate(env, **kw):
	env.Tool('addLibrary', library=['astro'], package = 'astro')
        depends = 'facilities tip'.split()
        for pack in depends: env.Tool(pack+'Lib')

	# why these?
	#env.Tool('addLibrary', library=env['clhepLibs'])
	#env.Tool('addLibrary', library=env['cfitsioLibs'])

def exists(env):
	return 1
