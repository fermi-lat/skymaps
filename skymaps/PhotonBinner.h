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

        /**@brief ctor user specified bins and nside to be used
        * @param edges left edges in MeV (note 1 TeV is always rightmost)
        * @param f_nside Healpix nside parameter for each energy band, front
        * @param b_nside Healpix nside parameter for each energy band, back
        */
        PhotonBinner(const std::vector<double>&edges, const std::vector<int>&f_nside, const std::vector<int>&b_nside);

        /**@brief ctor user specified bins and nside to be used. This version is for PSF event types
        * @param edges left edges in MeV (note 1 TeV is always rightmost)
        * @param psf0_nside Healpix nside parameter for each energy band, psf0
        * @param psf1_nside Healpix nside parameter for each energy band, psf1
        * @param psf2_nside Healpix nside parameter for each energy band, psf2
        * @param psf3_nside Healpix nside parameter for each energy band, psf3
        */
        PhotonBinner(const std::vector<double>&edges, 
            const std::vector<int>&psf0_nside, const std::vector<int>&psf1_nside,
            const std::vector<int>&psf2_nside, const std::vector<int>&psf3_nside);

        virtual ~PhotonBinner(){};
        ///@brief bin a photon by returning an appropriate Band object
        virtual skymaps::Band operator()(const astro::Photon& photon)const;
		
		///@brief get the key of the band a photon belongs in
		int get_band_key(const astro::Photon& photon);
        
        //@brief return an index for the event type: 0,1 for front, back, or 2,3,4,5 for PSF0..3 depending value 
        int event_type_index(int event_type)const;
        
        static void set_max_nside(unsigned int new_nside);
        static unsigned int  get_max_nside();

        static void set_min_nside(unsigned int new_nside);
        static unsigned int  get_min_nside();

        static void set_sigma_scale(double sigscale);
        static double get_sigma_scale();

    private:
        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void setupbins();
        double m_bins_per_decade;
        bool m_user_nside;      // should alwasy be set now
        bool m_psf_event_types; // switch between old front/back and new psf0..3
        
        std::vector<double> m_bins;           //the energy of each left bin edge
        std::vector<int> m_binlevelmap;       //the mapping between energy bins and healpix levels
        std::vector<int> m_fnside;           // nside parameters for front events
        std::vector<int> m_bnside;           // nside parameters for back events

        // New for PSF event types
        std::vector<int> m_psf0_nside;           // nside parameters for psf events
        std::vector<int> m_psf1_nside;           // nside parameters for psf events
        std::vector<int> m_psf2_nside;           // nside parameters for psf events
        std::vector<int> m_psf3_nside;           // nside parameters for psf events

        static std::vector<double> s_fenergy; //the mapping between energy and healpix levels for front events
        static std::vector<double> s_benergy; //the mapping between energy and healpix levels for back events
        static std::vector<double> s_sigma_level;
        static std::vector<double> s_gamma_level;
        static unsigned int max_nside;
        static unsigned int min_nside;
        static double m_sigma_scale;
        

        
    };
}

#endif
