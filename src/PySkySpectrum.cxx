/** @file PySkySpectrum.cxx
@brief implement class PySkySpectrum

$Header$
*/

#include "skymaps/PySkySpectrum.h"
#ifdef _DEBUG
#undef _DEBUG /* Link with python24.lib and not python24_d.lib */ 
#endif
#include <Python.h>

#include <stdexcept>

using namespace skymaps;

PySkySpectrum::PySkySpectrum(PyObject* value,PyObject * integral)
: m_value(value)
, m_integral(integral)
{
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw std::runtime_error("PySkySpectrum: error in ctor"); 
    }
}


double PySkySpectrum::value(const astro::SkyDir& dir, double e)const
{
    // setup the argument: a list with components
    CLHEP::Hep3Vector vdir( dir());
    PyObject* list = PyList_New(3);
    PyList_SetItem( list, 0, PyFloat_FromDouble(vdir[0]) );
    PyList_SetItem( list, 1, PyFloat_FromDouble(vdir[1]) );
    PyList_SetItem( list, 2, PyFloat_FromDouble(vdir[2]) );

    PyObject* args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, list);
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(e));

    // call the function
    PyObject* retobj  = PyObject_CallObject(m_value, args);
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw std::runtime_error("PySkySpectrum: error calling value function"); 
    }
     Py_DECREF(list);
     Py_DECREF(args);

    // get its returned value, convert to a float
    double ret(PyFloat_AsDouble(retobj));
    Py_DECREF(retobj);
    return ret;
}

double PySkySpectrum::integral(const astro::SkyDir& dir, double a, double b)const
{
    // setup the argument: a list with components
    CLHEP::Hep3Vector vdir( dir());
    PyObject* list = PyList_New(3);
    PyList_SetItem( list, 0, PyFloat_FromDouble(vdir[0]) );
    PyList_SetItem( list, 1, PyFloat_FromDouble(vdir[1]) );
    PyList_SetItem( list, 2, PyFloat_FromDouble(vdir[2]) );

    PyObject* args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, list);
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(a));
    PyTuple_SetItem(args, 2, PyFloat_FromDouble(b));

    // call the function
    PyObject* retobj  = PyObject_CallObject(m_integral, args);
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw std::runtime_error("PySkySpectrum: error calling value function"); 
    }
     Py_DECREF(list);
     Py_DECREF(args);

    // get its returned value, convert to a float
    double ret(PyFloat_AsDouble(retobj));
    Py_DECREF(retobj);
    return ret;
}
