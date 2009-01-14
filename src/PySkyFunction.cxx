/** @file PySkyFunction.cxx
@brief implement class PySkyFunction 

$Header$
*/

#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"
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
//copied from SkySpectrum.cxx by M. Kerr 14 January 2009
double PySkyFunction::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    using astro::SkyDir;
    using healpix::Healpix;

    static std::map<double, int> width_map;
    static bool map_built(false);

    int level, min_level = 6, max_level = 13;
    double result(0.0);

    // Get value for one point at center
    double previous = (*this) (dir);
    if (tolerance >= 0.5)  // If tolerance is higher than this, just return value at center.
        return previous;

    /* Build map of pixel widths by level, if not done yet.  Store twice the healpixel
       width for easy comparison. */
    if (!map_built)
    {
        width_map.clear();
        for (level = min_level; level <= max_level; ++level)
        {
            int nside(1 << level);
            int npix(12 * nside * nside);
            double width = sqrt(4 * M_PI / npix);  // Width of healpixel in radians
            width_map[2 * width] = level;
        }
        map_built = true;
    }

    // Use map to determine starting pixel level
    std::map<double, int>::iterator it = width_map.lower_bound(angle);
    if (it == width_map.end() || (it->first > angle && it != width_map.begin()))
        --it;
    level = it->second;

    // Get value for starting pixel level
    result = level_ave(dir, angle, level);

    // Iterate until result changes less than tolerance
    for(level += 1 ; fabs(result/previous -1.) > tolerance && level < max_level; ++ level)
    {
        previous = result;
        result = level_ave(dir, angle, level);
    }

    return result;

}

// Calculate average for a given level
double PySkyFunction::level_ave(const astro::SkyDir& dir, double angle, int level) const
{   

    int nside(1 << level);
    std::vector<int> v;
    healpix::Healpix hpx(nside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);
    hpx.query_disc(dir, angle, v); 
    double av(0);

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        healpix::HealPixel hp(*it, level);
        av += (*this) (hp());
    }

    return av/v.size();
}