/** @file SpectralFunction.h
@brief declaration of class SpectralFunction

$Header$
*/

#ifndef skymaps_SpectralFunction_h
#define skymaps_SpectralFunction_h

#include "astro/SkyDir.h"

#include <vector>

namespace skymaps{
class Band;
class SkySpectrum;

    /** @class SpectralFunction
        @brief evaluate a function of energy

        Note that the set_exposures method must be called before creating an object
    */
    class SpectralFunction {
    public:

        /** @brief ctor
        @param type The spectral function type
        @param pars specified parameters, appropriate for the type

        */
        SpectralFunction(const std::string& name, const std::vector<double>& pars);

        virtual ~SpectralFunction(){}

        /// @brief evaluate function itself
        double operator()(double energy)const{return value(energy);}

        /// @brief the function--note virtual to allow a polymorphic subclass to override

        virtual double value(double energy)const;

        /// @brief return expected counts in band
        ///
        /// @param dir the direction for exposure
        /// @param band the band object
        /// Note that this depends on exposure and effective area
        /// It will to the integral over the band
        double expected(const astro::SkyDir& dir,const skymaps::Band& band)const; 
    
        const std::string& name()const{return m_name;}
        
        std::vector<double>pars()const{return m_pars;}

        void set_parameters(const std::vector<double>& pars){m_pars=pars;}

        static void set_exposures(const skymaps::SkySpectrum* front, const skymaps::SkySpectrum* back);

        /// @brief access to exposure object
        static const skymaps::SkySpectrum* exposure(int n);
        static void set_simpson(int n);

    private:
        int m_index;
        std::string m_name; 
        std::vector<double> m_pars;
        void setup(); 
        
        static int s_n;///< Simpson's rule count
        static double s_e0; 
        static double s_flux_scale;


        /// list of exposure objects
        static std::vector<const skymaps::SkySpectrum*> s_exposures;

    };

}
#endif
