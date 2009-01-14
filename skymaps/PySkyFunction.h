/** @file PySkyFunction.h
@brief declare class PySkyFunction

$Header$

*/

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#include <string>

#ifndef skymaps_PySkyFunction_h
#define skymaps_PySkyFunction_h

namespace embed_python { class Module; }

// this needed to avoid include of Python.h, since PyObject is a struct
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif


namespace skymaps {


    /*** @class PySkyFunction
    @brief define a SkyFunction from a python module

    Note that a vector of the three components is passed to the python function,
    which may convert it back to a SkyDir thus:

        def __call__(self, v):
           return self.exposure(SkyDir(Hep3Vector(v[0],v[1],v[2])), [self.energy])


    */
    class PySkyFunction : public astro::SkyFunction{
    public:

        /** @brief ctor
        @param modulename name of python module 
        @param functionname name of a python function that takes arguments x,y,z, returns a function
        
        */
        PySkyFunction(std::string modulename, std::string functionname);


        /** @brief ctor
        @param callable a callable python object
        */
        PySkyFunction(PyObject* callable);

        ///@brief implement SkyFunction interface
        ///@param dir direction in sky
        ///@return value from the python function
        double operator()(const astro::SkyDir& dir)const;

        double average(const astro::SkyDir& dir, double angle, double tolerance)const;
        double level_ave(const astro::SkyDir& dir, double angle, int level) const;


    private:
        embed_python::Module* m_module;
        PyObject* m_fun;

    };


}

#endif

