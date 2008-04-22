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
        
        static EnergyBinner* Instance(bool combine);

        static EnergyBinner* Instance(std::vector<double>& bins);

        static EnergyBinner* Instance(double emin, double ratio, int bins);

        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void setupbins();
        
        /**@brief level  returns the healpix level for an energy bin
        */
        int level(int band, int event_class) const{return (band<0&&abs(event_class)<2)?-1:m_binlevelmap[2*band+event_class];}

        /**@brief level  returns the healpix level for a photon
        */
        int level(const astro::Photon& p) const;
        
        /**@brief bin  returns the energy bin for a photon
        */
        int band(const astro::Photon& p) const;

        /**@brief ebins  returns the energy bins for output
        */
        std::vector<double> ebins() const {return m_bins;}

        int bands() const {return m_bins.size();}

        bool comb() {return m_comb;}

        double sigma(int level);

        double gamma(int level);

    protected:
        /**@brief ctor  old style of combining front/back
        */
        EnergyBinner(bool combine);

        /** @brief ctor  takes arguments of the left edge bin energy in MeV
        */
        EnergyBinner(std::vector<double>& bins);

        /**@brief ctor  takes arguments for a power law binning 
        * @param emin  minimum energy in MeV
        * @param emax  ratio between bins
        * @param bins  number of bins
        */
        EnergyBinner(double emin, double ratio, int bins);

    private:
        bool m_comb;                          //old style combine front/back events?
        std::vector<double> m_bins;           //the energy of each left bin edge
        std::vector<int> m_binlevelmap;       //the mapping between energy bins and healpix levels
        static std::vector<double> s_fenergy; //the mapping between energy and healpix levels for front events
        static std::vector<double> s_benergy; //the mapping between energy and healpix levels for back events
        static std::vector<double> s_sigma_level;
        static std::vector<double> s_gamma_level;
        static EnergyBinner* s_instance;       //Singleton
    };


}