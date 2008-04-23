/** @file PhotonBinner.h
@brief declare class PhotonBinner

$Header$

*/
#ifndef skymaps_PhotonBinner_h
#define skymaps_PhotonBinner_h

//#include "skymaps/BinnedPhoton.h"

namespace astro {class Photon;}
namespace skymaps {class BinnedPhoton;}

#include <string>
#include <vector>
#include <map>

#include "skymaps/Band.h"

namespace skymaps {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class PhotonBinner
    @brief 

    */
    class PhotonBinner { 
    public:
         /**@brief ctor  old style of combining front/back
        */
        PhotonBinner(bool combine=true);

        /** @brief ctor  takes arguments of the left edge bin energy in MeV
        */
        PhotonBinner(std::vector<double>& bins);

        /**@brief ctor  takes arguments for a power law binning 
        * @param emin  minimum energy in MeV
        * @param emax  ratio between bins
        * @param bins  number of bins
        */
        PhotonBinner(double emin, double ratio, int bins);

        skymaps::BinnedPhoton operator()(const astro::Photon& photon)const;

        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void setupbins();
        
        /**@brief level  returns the healpix level for an energy bin
        */
        int level(int band, int event_class) const;
        
        /**@brief ebins  returns the energy bins for output
        */
        std::vector<double> ebins() const {return m_bins;}

        int bands() const {return m_bins.size();}

        bool comb() {return m_comb;}

        double sigma(int level);

        double gamma(int level);

    private:

        ///! add a new band, if unique, return pointer to transient object
        const Band* addBand(const skymaps::Band& band)const;

        /// Maintain a list of bands: key is an int formed from the band info
        static std::map<int, skymaps::Band> s_bands;

        bool m_comb;                          //old style combine front/back events?
        std::vector<double> m_bins;           //the energy of each left bin edge
        std::vector<int> m_binlevelmap;       //the mapping between energy bins and healpix levels
        static std::vector<double> s_fenergy; //the mapping between energy and healpix levels for front events
        static std::vector<double> s_benergy; //the mapping between energy and healpix levels for back events
        static std::vector<double> s_sigma_level;
        static std::vector<double> s_gamma_level;
    };
}

#endif
