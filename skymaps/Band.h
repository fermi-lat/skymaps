/** @file Band.h
@brief declare class Band

$Header$

*/

#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"

#ifndef skymaps_Band_h
#define skymaps_Band_h

#include <vector>
#include <map>

namespace healpix { class Healpix;}

namespace skymaps {

    /*** @class Band
        @brief encapsulate concept of an energy band and a collection of directions for photons in that band

        Implement SkyFunction to return no of entries in a bin in the given direction
        Implemnt  a map of pixel ids and contents.
        
    */
    class Band : public astro::SkyFunction{
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
            double sigma, double gamma);

        ///@brief implement SkyFunction interface
        ///@param dir direction in sky
        ///@return contents of pixel, if exists, otherwise zero
        double operator()(const astro::SkyDir& dir)const;

        ///@brief add an element by direction
        void add(const astro::SkyDir& dir, int count=1);

        ///@brief add an element, or to an existing element, by index and count
        void add(int index, int count);

        /// @brief direction for a pixel index
        astro::SkyDir dir( int index)const;

        /// @brief the pixel index from a direction
        int index(const astro::SkyDir& dir)const;

        /// @brief set a list of pixel ids and counts within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<astro::SkyDir,int> > & v)const;
 
        /// @brief set a list of pixel ids and counts within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<int,int> > & v)const;

        /// @brief the solid angle for this pixelization
        double pixelArea()const;

        /// @brief for identity,sorting: assume emin and event class is unique
        operator int()const{return m_event_class +10* static_cast<int>(m_emin+0.5);}

        ///@brief count of photons
        int photons()const; 
        int nside()const { return m_nside; }
        int event_class()const{return m_event_class; }
        double emin()const{return m_emin;}
        double emax()const{return m_emax;}
        double sigma()const{return m_sigma;}
        double gamma()const{return m_gamma;}

        typedef std::map<int,int> PixelMap;
        typedef PixelMap::iterator iterator;
        typedef PixelMap::const_iterator const_iterator;
        iterator begin(){return m_pixels.begin();}
        iterator end(){return m_pixels.begin();}
        const_iterator begin()const{return m_pixels.begin();}
        const_iterator end()const{return m_pixels.begin();}
        size_t size()const{return m_pixels.size();}

    private:
        PixelMap m_pixels; 
        int m_nside;
        int m_event_class;
        double m_emin, m_emax;
        double m_sigma, m_gamma;
        const healpix::Healpix* m_healpix; 
    };

}

#endif

