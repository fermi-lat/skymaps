/** @file Band.h
@brief declare class Band

$Header$

*/


#include "astro/SkyDir.h"

#ifndef skymaps_Band_h
#define skymaps_Band_h

namespace skymaps {

    /*** @class Band
        @brief encapsulate concept of an energy band

    */
    class Band {
    public:
        /// default
        Band(): m_nside(-1), m_event_class(0){}

        Band(int nside, int event_class, double emin,double emax,
            double sigma, double gamma)
            : m_nside(nside)
            , m_event_class(event_class)
            , m_emin(emin)
            , m_emax(emax)
            , m_sigma(sigma)
            , m_gamma(gamma)
        {}
        astro::SkyDir dir(unsigned int index)const;
        unsigned int index(const astro::SkyDir& dir)const;

        /// @brief for identity,sorting: assume nside and event class is unique
        operator int()const{return (m_event_class ) | m_nside<<16;}

        int nside()const { return m_nside; }
        int event_class()const{return m_event_class; }
        double emin()const{return m_emin;}
        double emax()const{return m_emax;}
        double sigma()const{return m_sigma;}
        double gamma()const{return m_gamma;}
        double pixel_area()const;
    private:
        int m_nside;
        int m_event_class;
        double m_emin, m_emax;
        double m_sigma, m_gamma;
    };

}

#endif

