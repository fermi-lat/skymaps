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
        Implemnt a map of pixel ids and contents.
        
    */
    class Band : public astro::SkyFunction{
    public:

        /// @brief default ctor for pixel conversons
        Band(int nside=1);

        /** @brief ctor
            @param nside The HEALpix nside parameter. Total number of pixels will be 12*nside*nside
            @param event_class The LAT event class number, 0/1 for front/back currently
            @param emin,emax energy range
            @param sigma The sigma parameter, a scale factor for the PSF
            @param gamma the power law for the PSF
        */
        Band(int nside, int event_class, double emin,double emax,
            double sigma, double gamma);

        ///! copy ctor
        //Band(const Band& other);

        ///@brief implement SkyFunction interface
        ///@param dir direction in sky
        ///@return contents of pixel, if exists, otherwise zero
        double operator()(const astro::SkyDir& dir)const;

        ///@brief add an element by direction
        void add(const astro::SkyDir& dir, int count=1);

        ///@brief add an element, or to an existing element, by index and count
        void add(int index, int count);

        ///@brief add the contents of another Band (must have same parameters)
        void add(const Band& other);

        /// @brief direction for a pixel index
        astro::SkyDir dir( int index)const;

        /// @brief the pixel index from a direction
        int index(const astro::SkyDir& dir)const;

        /// @brief set a list of pixel ids and counts within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        /// @param v vector of SkyDir objects to be set
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<astro::SkyDir,int> > & v)const;
 
        /// @brief set a list of pixel ids and counts within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        /// @param v vector of pixel indices to be set
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<int,int> > & v)const;


        /// @brief fill a vector indeces of 7 or 8 neighbors of given index
        /// Note: requires nested indexing!
        void findNeighbors(int index, std::vector<int> &neighbors)const;


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

        void setSigma(double s){m_sigma=s;}
        void setGamma(double g){m_gamma=g;}

        typedef std::map<int,int> PixelMap;
        typedef PixelMap::iterator iterator;
        typedef PixelMap::const_iterator const_iterator;
        iterator begin(){return m_pixels.begin();}
        iterator end(){return m_pixels.begin();}
        const_iterator begin()const{return m_pixels.begin();}
        const_iterator end()const{return m_pixels.end();}
        size_t size()const{return m_pixels.size();}

    private:
        PixelMap m_pixels; 
        int m_nside;
        int m_event_class;
        double m_emin, m_emax;
        double m_sigma, m_gamma;
        const healpix::Healpix* m_healpix; 
    };

   /** @class WeightedSkyDir
         @brief a weighted SkyDir, used to describe pixels
    */
    class WeightedSkyDir : public astro::SkyDir {
    public:
        WeightedSkyDir(const astro::SkyDir& sdir=astro::SkyDir(), double weight=1): SkyDir(sdir), m_weight(weight){}
        double weight()const{return m_weight;}
    private:
        double m_weight;
    };


    /** @class WeightedSkyDirList
        @brief a vector of WeightedSKyDir objects

        */
    class WeightedSkyDirList : public astro::SkyFunction, public std::vector<WeightedSkyDir> {
    public:
        /** @brief ctor creates an object from a band
        */
        WeightedSkyDirList( const Band& band, const astro::SkyDir& sdir, double radius);

        /// @brief implement SkyFunction: return weight, or zero
        double operator()(const astro::SkyDir& sdir)const;

    private:
        const skymaps::Band& m_band; ///< the band used to make
    };


}

#endif

