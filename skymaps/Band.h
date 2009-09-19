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
            double sigma, double gamma, double sigma2=-1, double gamma2=-1, double frac2=0);

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
        /// @param include_empty [false] add empty pixels to the returned list
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<astro::SkyDir,int> > & v, bool include_empty=false)const;
 
        /// @brief set a list of pixel ids and counts within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        /// @param v vector of pixel indices to be set
        /// @param include_empty [false] add empty pixels to the returned list
        /// @return the number of photons 
        int query_disk(const astro::SkyDir&dir, double radius, 
            std::vector<std::pair<int,int> > & v, bool include_empty=false)const;

        /// @brief return total pixels within the radius about the direction
        /// @param dir center of cone
        /// @param radius in radians
        int total_pix(const astro::SkyDir&dir, double radius)const;

        /// @brief fill a vector indeces of 7 or 8 neighbors of given index
        /// Note: requires nested indexing!
        void findNeighbors(int index, std::vector<int> &neighbors)const;

        /// @brief return the photon density at the given direction
        double density(const astro::SkyDir&dir, bool smooth=false, int mincount = 0, int kernel = 0, double smooth_radius = 3)const;


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
        double sigma2()const{return m_sigma2;}
        double gamma2()const{return m_gamma2;}
        double frac2() const { return m_frac2;}

        static int cache_pix() {return cache_pix_counts;}
        static bool enable_cache(bool b);
  
        void setSigma(double s){m_sigma=s;}
        void setGamma(double g){m_gamma=g;}
        void setSigma2(double s){m_sigma2=s;}
        void setGamma2(double g){m_gamma2=g;}
        void setFrac2(double f){m_frac2=f;}

        typedef std::map<int,int> PixelMap;
        typedef PixelMap::iterator iterator;
        typedef PixelMap::const_iterator const_iterator;
        iterator begin(){return m_pixels.begin();}
        iterator end(){return m_pixels.begin();}
        const_iterator begin()const{return m_pixels.begin();}
        const_iterator end()const{return m_pixels.end();}
        size_t size()const{return m_pixels.size();}

        /// @access photon source ids
        void add_source(int source_id);
        std::vector<std::pair<int,int> > source()const{return m_source;}
        
        size_t density_cache_size() {return m_density_cache.size();}
        
    private:
        PixelMap m_pixels; 
        int m_nside;
        int m_event_class;
        double m_emin, m_emax;
        double m_sigma, m_gamma;
        double m_sigma2, m_gamma2, m_frac2;
        const healpix::Healpix* m_healpix;
        static int cache_pix_counts;
        static bool m_enable_cache;
        mutable std::map<int,double> m_density_cache;

        // @photon source ids
        std::vector<std::pair<int,int> > m_source;
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
        @param sdir direction for cone
        @param radius cone half-angle (radians)
        @param include_empty [false] include empty pixels
        */
        WeightedSkyDirList( const Band& band, const astro::SkyDir& sdir, double radius, bool include_empty=false);

        /// @brief implement SkyFunction: return weight, or zero
        double operator()(const astro::SkyDir& sdir)const;

        /// @brief total pixels in list, including empties
        int total_pix()const {return m_pix;}
        int counts()const {return m_counts;}

    private:
        const skymaps::Band& m_band; ///< the band used to make
        int m_pix;
        int m_counts;
    };


}

#endif

