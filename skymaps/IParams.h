/** @file IParams.h
    @brief declare class IParams

$Header$

*/
#ifndef skymaps_IParams_h
#define skymaps_IParams_h

#include <string>
#include <vector>

namespace skymaps {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IParams
    @brief a class that manages the parameters of the PSF from a database
*/

    class IParams {
public:
    
    /** @brief sigma returns PSF parameter sigma
        @param energy 
        @param event_class 0=front, 1=back
        */
    static double sigma(double energy, int event_class);

    /** @brief gamma returns PSF parameter gamma
        @param energy 
        @param event_class 0=front, 1=back
        */
    static double gamma(double energy, int event_class);

    //sets front and back parameterization for sigma
    // s = (a[0]^2+a[2]*(E/100)^(a[1])+a[3]*(E/100)^(2*a[1]))^0.5
    static void set_fp(std::vector<double> params);
    static void set_bp(std::vector<double> params);

    //sets gamma values for front and back set for 4/decade from 100MeV->100GeV
    //with a bin for <100 MeV and one for >100GeV (15 in all)
    static void set_fgam(std::vector<double> gams);
    static void set_bgam(std::vector<double> gams);

    //set energy bins
    static void set_elist(std::vector<double> elist);

    //run the initialization through the embed python module
    //the file is psf_defaults.py
    //the default reads from CALDB files
    static void init();

    //returns 
    static std::vector<double> params(int event_class);

private:

    /** @brief ctor not used
    */
    IParams(std::string& name);

    static std::vector<double> s_elist; //energy bin list 4/decade
    static std::vector<double> s_fparams; //front width parameterization
    static std::vector<double> s_bparams; //back width parameterization
    static std::vector<double> s_fgam; //front tail parameterization
    static std::vector<double> s_bgam; //back tail parameterization

    static bool s_init; //checks to see whether initialization has been performed
};

}
#endif

