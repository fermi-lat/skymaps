/** @file PySkySpectrum.h
@brief declare class PySkySpectrum

$Header$

*/

#include "astro/SkyDir.h"
#include "skymaps/SkySpectrum.h"

#include <string>

#ifndef skymaps_PySkySpectrum_h
#define skymaps_PySkySpectrum_h

// this needed to avoid include of Python.h, since PyObject is a struct
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif


namespace skymaps {
    /*** @class PySkySpectrum
    @brief define a SkySpectrum from a python module

    Note that a three components list and a scaler 
    is passed to the python function 'value'
    which may convert it back to a SkyDir thus:

        def value(v,energy):
           dir=Hep3Vector(v[0],v[1],v[2]))

    A three component list and the lower and upper energy
    are passed to the python function 'integral'
    
        def integral(v,elow,ehigh):
           dir=Hep3Vector(v[0],v[1],v[2]))

    */
    class PySkySpectrum : public skymaps::SkySpectrum {
    public:

        /** @brief ctor
        @param value a callable python object which should implement 
        the value function for a SkySpectrum
        @param integral a callable python object which should implement
        the integral function for a SkySpectrum
        */
        PySkySpectrum(PyObject* value,PyObject * integral);

        ///@brief implement SkySpectrum function 'value'
        ///@param dir direction in sky
        ///@param e energy
        ///@return value of the python function 'value'
        double value(const astro::SkyDir& dir, double e)const;

        ///@brief implement SkySpectrum function 'integral'
        ///@param dir direction in sky
        ///@param a lower energy
        ///@param a upper energy
        ///@return value from the python function 'integral'
        double integral(const astro::SkyDir& dir, double a, double b)const;

    private:
        PyObject* m_value;
        PyObject* m_integral;

    };


}

#endif

