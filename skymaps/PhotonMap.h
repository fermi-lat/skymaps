/** @file PhotonMap.h
@brief declare class PhotonMap

$Header$

*/
#ifndef skymaps_PhotonMap_h
#define skymaps_PhotonMap_h

#include "skymaps/SkySpectrum.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"
#include "skymaps/Gti.h"
#include <map>

namespace astro {class Photon;}

#include <string>

namespace skymaps {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class PhotonMap
    @brief a SkySpectrum that adapts the class map_tools::PhotonMap. 

    */

    class PhotonMap : public skymaps::SkySpectrum
        ,  public std::map<healpix::HealPixel, unsigned int> {
    public:
        /** @brief ctor defines energy binning
        @param emin [100] Minimum energy
        @param eratio [2.35] ratio between bins
        @param nlevels [8] number of levels
        @param minlevel [6] HealPixel level of first bin

        */
        PhotonMap(double emin=100, double eratio=2.35, int nlevels=8, int minlevel=6);

        /// Create PhotonMap object from saved fits file
        PhotonMap(const std::string & inputFile, const std::string & tablename);

        ///@brief data value for bin with given energy 
        ///@param e energy in MeV
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief integral for the energy limits, in the given direction
        /// Assume that request is for an energy bin.
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        virtual std::string name()const{return m_name;}

        ///@ allow different name
        void setName(const std::string name){m_name = name;}

        /// add a photon to the map with the given energy and direction
        void addPhoton(const astro::Photon& gamma);

        /// add a healpixel to the map with an associated count
        void addPixel(const healpix::HealPixel & px, int count);

        /// @return density for a given direction, in photons/area of the base pixel.
        double density (const astro::SkyDir & sd) const;

        //! Count the photons within a given pixel.
        double photonCount(const healpix::HealPixel & px, bool includeChildren=false,
            bool weighted=false) const;

        //! Count the photons within a given pixel, weighted with children.  Also return weighted direction.
        double photonCount(const healpix::HealPixel & px, astro::SkyDir & NewDir) const;

        ///  implement the SkyFunction class by returning density
        double operator()(const astro::SkyDir & sd) const{ return density(sd);}
#ifndef SWIG // don't understand why this causes problems
        ///   combines photonmaps and their gti's
        void operator+=(const skymaps::PhotonMap& other);
#endif

        /// the binning function: return a HealPixel corresponding to the 
        /// direction and energy
        healpix::HealPixel pixel(const astro::Photon& gamma);

        /** @brief extract a subset around a given direction. include selected pixels and all children.
        @param radius The maximum radius (deg). Set to >=180 for all
        @param vec the vector to fill with (healpixel, count ) pairs
        @param summary_level [-1] selection level: default is the minimum level 
        @param select_level [-1] set to return only pixels at this level
        @return the total number of photons (sum of count)
        */
        int extract(const astro::SkyDir& dir, double radius,
            std::vector<std::pair<healpix::HealPixel, int> >& vec,
            int summary_level = -1, int select_level = -1) const;

        /** @brief extract a subset around a given direction.  single level only.
        @param radius The maximum radius (deg). Set to >=180 for all
        @param vec the vector to fill with (healpixel, count ) pairs
        @param select_level [-1] return only pixels at this level. default is the minimum level
        @param include_all [false] True: return all possible select_level pixels within radius.  False: return only pixels found in current PhotonMap.
        @return the total number of photons (sum of count)
        */

        int extract_level(const astro::SkyDir& dir, double radius,
            std::vector<std::pair<healpix::HealPixel, int> >& vec,
            int select_level = -1, bool include_all = false) const;

        int photonCount()const { return m_photons;} ///< current number of photons
        int pixelCount()const { return m_pixels; } ///< current nubmer of pixesl

        int minLevel()const { return m_minlevel;} ///< minimum Healpixel level
        int levels()const {return m_levels;};  ///< number of energy bins

        /// @return a vector of the left edges of the energy bins
        std::vector<double> energyBins()const;


        /**@brief Write a PhotonMap object to a fits file
        @param outputFile Fully qualified fits output file name
        @param tablename Fits secondary extension name
        @param clobber Whether to delete an existing file first 
        */
        void write(const std::string & outputFile,
            const std::string & tablename="PHOTONMAP",
            bool clobber = true) const;

        /**@brief add GTI info to the current gti
        */
        void addgti(const skymaps::Gti& other);

        /**@brief Write a PhotonMap gti info to a fits file
        @param outputFile Fully qualified fits output file name
        */
        void writegti(const std::string & outputFile) const;
        /// @return a modifyable reference to gti info
        skymaps::Gti & gti() {return m_gti;};
#ifndef SWIG
        /// @return a const reference to the gti
        const skymaps::Gti & gti()const {return m_gti;};
#endif


    private:
        std::string m_name; ///< use name of file as descriptive name
        double m_emin;     ///< minimum energy for first bin
        double m_logeratio;   ///< log10(ratio between energy bins)
        int    m_levels;   ///< number of levels to create
        int    m_minlevel; ///< level number for first pixel 
        int    m_photons;  ///< total number of photons added
        int    m_pixels;   ///< keep track of total number of pixels
        skymaps::Gti m_gti;   ///< gti information associated with this photon map
        };


}

#endif
