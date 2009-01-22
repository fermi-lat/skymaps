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
class Exposure;

    /** @class SpectralFunction
        @brief evaluate a function of energy


    */
    class SpectralFunction {
    public:
        typedef enum {PowerLaw, ExpCutoff, BrokenPowerLaw} Type;

        SpectralFunction(Type type, const std::vector<double>& pars);

        /// @brief evaluate function itself
        double operator()(double energy)const{return value(energy);}
        /// @brief the function

        double value(double energy)const;

        /// @brief return expected counts in band
        ///
        /// @param dir the direction for exposure
        /// @param band the band object
        /// Note that this depends on exposure and effective area
        /// It will to the integral over the band
        double expected(const astro::SkyDir& dir,const skymaps::Band& band)const; 
    

        Type type()const{return m_type;}
        
        std::vector<double>pars()const{return m_pars;}

        void set_parameters(const std::vector<double>& pars){m_pars=pars;}

        static void set_exposures(const skymaps::Exposure* front, const skymaps::Exposure* back);

        /// @brief access to exposure object
        static const skymaps::Exposure* exposure(int n);
        static void set_simpson(int n);

    private:
        Type m_type;
        std::vector<double> m_pars;
        
        static int s_n;///< Simpson's rule count
        static double s_e0; 
        static double s_flux_scale;


        /// list of exposure objects
        static std::vector<const skymaps::Exposure*> s_exposures;

    };

}
#endif
