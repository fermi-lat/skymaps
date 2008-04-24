/** @file Band.h
@brief declare class Band

$Header$

*/


#include "astro/SkyDir.h"

#ifndef skymaps_Band_h
#define skymaps_Band_h

#include <vector>

namespace skymaps {

    /*** @class Band
        @brief encapsulate concept of an energy band

    */
    class Band {
    public:
        /// default
        Band(): m_nside(-1), m_event_class(0){}

        /** @brief ctor
            @param nside The HEALpix nside parameter. Total number of pixels will be 12*nside*nside
            @param event_class The LAT event class number, 0/1 for front/back currently
            @param emin,emax energy range
            @param sigma The sigma parameter, a scale factor for the PSF
            @param gamma the power law for the PSF

        */
        Band(int nside, int event_class, double emin,double emax,
            double sigma, double gamma)
            : m_nside(nside)
            , m_event_class(event_class)
            , m_emin(emin)
            , m_emax(emax)
            , m_sigma(sigma)
            , m_gamma(gamma)
        {}
        /// @brief direction for a pixel id
        astro::SkyDir dir(unsigned int index)const;

        /// @brief the pixel index from a direction
        unsigned int index(const astro::SkyDir& dir)const;

        /// @brief set a list of pixel id's within the radius about the direction
        void query_disk(const astro::SkyDir&dir, double radius, std::vector<int>& v)const;

        /// @brief the solid angle for this pixelization
        double pixelArea()const;

        /// @brief for identity,sorting: assume nside and event class is unique
        operator int()const{return m_event_class +10* m_nside;}

        int nside()const { return m_nside; }
        int event_class()const{return m_event_class; }
        double emin()const{return m_emin;}
        double emax()const{return m_emax;}
        double sigma()const{return m_sigma;}
        double gamma()const{return m_gamma;}
    private:
        int m_nside;
        int m_event_class;
        double m_emin, m_emax;
        double m_sigma, m_gamma;
    };

}

#endif

