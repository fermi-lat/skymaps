#$Id$
def generate(env, **kw):
    env.Tool('addLibrary', library = 'healpix')
    env.Tool('astroLib')
    
def exists(env):
    return 1;
