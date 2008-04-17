/** @file EnergyBinner.h
    @brief declare class EnergyBinner

$Header$

*/

#include <vector>
#include <map>

namespace astro {
    class Photon;
}

namespace skymaps {

    /**
    @class EnergyBinner
    @brief generic energy binning scheme, defines connection to HEALPix representation

    */
    class EnergyBinner {
    
    public:
        
        /**@brief ctor  old style of combining front/back
        */
        EnergyBinner(bool combine);

        /** @brief ctor  takes arguments of the left edge bin energy in MeV
        */
        EnergyBinner(std::vector<double>& bins);

        /**@brief ctor  takes arguments for a power law binning 
        * @param emin  minimum energy in MeV
        * @param emax  max energy in MeV
        * @param bins  number of bins
        */
        EnergyBinner(double eminf, double emaxf, double eminb, double emaxb, int bins);

        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void setupbins();
        
        /**@brief level  returns the healpix level for an energy bin
        */
        int level(int band) {return band<0?-1:m_binlevelmap[band];}

        /**@brief level  returns the healpix level for a photon
        */
        int level(const astro::Photon& p);
        
        /**@brief bin  returns the energy bin for a photon
        */
        int band(const astro::Photon& p);

        /**@brief ebins  returns the energy bins for output
        */
        std::vector<double> ebins() {return m_bins;}

    private:
        bool m_comb;                          //combine front/back events?
        std::vector<double> m_bins;           //the energy of each left bin edge
        std::vector<int> m_binlevelmap;       //the mapping between energy bins and healpix levels
        static std::vector<double> s_fenergy; //the mapping between energy and healpix levels for front events
        static std::vector<double> s_benergy; //the mapping between energy and healpix levels for back events
        
    };


}