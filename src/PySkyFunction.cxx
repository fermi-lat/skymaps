/** @file PySkyFunction.cxx
@brief implement class PySkyFunction 

$Header$
*/

#include "skymaps/PySkyFunction.h"
#include "embed_python/Module.h"
#ifdef _DEBUG
#undef _DEBUG /* Link with python24.lib and not python24_d.lib */ 
#endif
#include <Python.h>

#include <stdexcept>

using namespace skymaps;
using namespace embed_python;


PySkyFunction::PySkyFunction(std::string modulename, std::string functionname)
: m_module(&embed_python::Module("" , modulename))
, m_fun(embed_python::Module("" , modulename).attribute(functionname) )
{

}

PySkyFunction::PySkyFunction(PyObject* callable)
: m_fun(callable)
{
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw std::runtime_error("PySkyFunction: error in ctor"); 
    }
}


double PySkyFunction::operator()(const astro::SkyDir& dir)const
{
    // setup the argument: a list with components
    Hep3Vector vdir( dir());
    PyObject* list = PyList_New(3);
    PyList_SetItem( list, 0, PyFloat_FromDouble(vdir[0]) );
    PyList_SetItem( list, 1, PyFloat_FromDouble(vdir[1]) );
    PyList_SetItem( list, 2, PyFloat_FromDouble(vdir[2]) );
    PyObject* args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, list);

    // call the function
    PyObject* retobj  = PyObject_CallObject(m_fun, args);
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw std::runtime_error("PySkyFunction: error calling function"); 
    }
     Py_DECREF(list);
     Py_DECREF(args);

    // get its returned value, convert to a float
    double ret(PyFloat_AsDouble(retobj));
    Py_DECREF(retobj);
    return ret;
}
