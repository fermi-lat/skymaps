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
#include <map>
#include <iostream>

namespace astro {class Photon;}

#include <string>

namespace skymaps {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class BinnedPhotonData
    @brief a SkySpectrum that adapts the class map_tools::BinnedPhotonData. 

    */

    class BinnedPhotonData : public skymaps::SkySpectrum 
        , public std::map<int, std::map<unsigned int, unsigned int> > {

    public:
        ///@brief ctor created from binner
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


        /** @brief extract a subset around a given direction, corresponding to the band
        @param bin corresponds a photon with the desired band, used for center
        @param radius The maximum radius (deg). Set to >=180 for all
        @param vec the vector to fill with (ident, count ) pairs
         @return the total number of photons (sum of count)
        */
        int extract(const BinnedPhoton& bin, double radius,
            std::vector<std::pair<unsigned int, unsigned int> > & vec) const;

        void info(std::ostream& out = std::cout)const;
    private:
        skymaps::PhotonBinner m_binner; ///< object that handles binning
        skymaps::Gti m_gti;   ///< gti information associated with these data
        int m_photons;  ///< keep track of total number of photons in the database

    };

}


#endif
