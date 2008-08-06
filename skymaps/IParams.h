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
    IParams(std::string& name);
    
    static double sigma(double energy, int event_class);
    static double gamma(double energy, int event_class);

    static void set_fp(std::vector<double> params){s_fparams=params;}
    static void set_bp(std::vector<double> params){s_bparams=params;}


    static std::vector<double> params(int event_class);

private:

    static std::vector<double> s_elist;
    static std::vector<double> s_fparams;
    static std::vector<double> s_bparams;

};

}
#endif

