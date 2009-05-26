/** @file PhotonBinner.h
@brief declare class PhotonBinner

$Header$

*/
#ifndef skymaps_PhotonBinner_h
#define skymaps_PhotonBinner_h

//#include "skymaps/BinnedPhoton.h"

namespace astro {class Photon;}
//namespace skymaps {class BinnedPhoton;}

#include <string>
#include <vector>
#include <map>


#include "skymaps/Band.h"

namespace skymaps {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class PhotonBinner
    @brief manage binning of photons 
    */
    class PhotonBinner { 
    public:
         /**@brief ctor (default) 
        */
        PhotonBinner(double bins_per_decade=0);

        /** @brief ctor  takes arguments of the left edge bin energy in MeV
        */
        PhotonBinner(const std::vector<double>& bins);

        /**@brief ctor  takes arguments for a power law binning 
        * @param emin  minimum energy in MeV
        * @param ratio  ratio between bins
        * @param bins  number of bins
        */
        PhotonBinner(double emin, double ratio, int bins);
        virtual ~PhotonBinner(){};
        ///@brief bin a photon by returning an appropriate Band object
        virtual skymaps::Band operator()(const astro::Photon& photon)const;
           
        static void set_max_nside(int new_nside);
        static int  get_max_nside();

        static void set_sigma_scale(double sigscale);
        static double get_sigma_scale();

    private:
        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void setupbins();
        double m_bins_per_decade;
        
        std::vector<double> m_bins;           //the energy of each left bin edge
        std::vector<int> m_binlevelmap;       //the mapping between energy bins and healpix levels
        static std::vector<double> s_fenergy; //the mapping between energy and healpix levels for front events
        static std::vector<double> s_benergy; //the mapping between energy and healpix levels for back events
        static std::vector<double> s_sigma_level;
        static std::vector<double> s_gamma_level;
        static int max_nside;
        static double m_sigma_scale;
    };
}

#endif
