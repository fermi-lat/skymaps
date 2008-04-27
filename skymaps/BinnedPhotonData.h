/** @file BinnedPhotonData.h
@brief declare class BinnedPhotonData

$Header$

*/
#ifndef skymaps_BinnedPhotonData_h
#define skymaps_BinnedPhotonData_h

#include "skymaps/SkySpectrum.h"
#include "skymaps/PhotonBinner.h"
#include "skymaps/Band.h"
#include "skymaps/Gti.h"
#include <list>
#include <iostream>

namespace astro {class Photon;}

#include <string>

namespace skymaps {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class BinnedPhotonData
    @brief a SkySpectrum containing the photon data used by pointlike

    */

    class BinnedPhotonData : public skymaps::SkySpectrum , public std::list<Band> {

    public:
        ///@brief default ctor: will use default binner
        BinnedPhotonData();
        
        ///@brief ctor created from external binner
        BinnedPhotonData(const skymaps::PhotonBinner& binner);

        /// Create  object from a saved fits file
        BinnedPhotonData(const std::string & inputFile, const std::string & tablename);

        ///@brief data value for bin with given energy 
        ///@param e energy in MeV
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief integral for the energy limits, in the given direction
        /// Assume that request is for an energy bin.
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        /// add a photon to the map with the given energy and direction
        void addPhoton(const astro::Photon& gamma);

        /// @return density for a given direction, in photons/area of the base pixel.
        double density (const astro::SkyDir & sd) const;

        /// @brief print out a summary of the contents
        void info(std::ostream& out = std::cout)const;
        
        /**@brief Write  to a fits file
        @param outputFile Fully qualified fits output file name
        @param tablename Fits secondary extension name
        @param clobber Whether to delete an existing file first 
        */
        void write(const std::string & outputFile,
            const std::string & tablename="BINNEDPHOTONS",
            bool clobber = true) const;

        
        /**@brief add GTI info to the current gti
        */
        void addgti(const skymaps::Gti& other);

        /**@brief Write a BinnedPhotonData gti info to a fits file
        @param outputFile Fully qualified fits output file name
        */
        void writegti(const std::string & outputFile) const;
        /// @return a modifyable reference to gti info
        skymaps::Gti & gti() {return m_gti;};

#ifndef SWIG  // don't understand why these cause SWIG problems

        ///   combines photonmaps and their gti's
        void operator+=(const skymaps::BinnedPhotonData& other);

        /// @return a const reference to the gti
        const skymaps::Gti & gti()const {return m_gti;};
#endif

        const PhotonBinner& bands()const {return m_binner;}

        int photonCount()const{return m_photons;}
        int pixelCount()const{ return 0;} ///TODO

    private:
        skymaps::PhotonBinner m_binner; ///< object that handles binning, has all the data
        skymaps::Gti m_gti;   ///< gti information associated with these data
        int m_photons;  ///< keep track of total number of photons in the database

    };

}


#endif
