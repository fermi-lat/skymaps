/** @file IParams.cxx
    @brief implement class IParams

$Header$

*/

#include "skymaps/IParams.h"
#include <cmath>

using namespace skymaps;

namespace {
    
    //kluge - all gamma parameterization
    double f_p[] = {0.01165278, -1.45237162, 1.67211108, 2.28589082};
    double b_p[] = {0.02614202, -1.55092155, 5.18346847, 3.68094493};

}

std::vector<double> IParams::s_fparams(f_p,f_p+4);
std::vector<double> IParams::s_bparams(b_p,b_p+4);

IParams::IParams(std::string& name)
{
   //TODO: setup database, for now use defaults from allGamma v14r8 
}

double IParams::sigma(double energy, int event_class){
    double deg = 0;
    if(event_class) {
        //from parameterization of PSF function
        deg = sqrt(IParams::s_bparams[0]*IParams::s_bparams[0]+IParams::s_bparams[2]*pow(energy/100,IParams::s_bparams[1])+IParams::s_bparams[3]*pow(energy/100,2*IParams::s_bparams[1]));
    } else {
        deg = sqrt(IParams::s_fparams[0]*IParams::s_fparams[0]+IParams::s_fparams[2]*pow(energy/100,IParams::s_fparams[1])+IParams::s_fparams[3]*pow(energy/100,2*IParams::s_fparams[1]));
    }
    return deg*M_PI/180;
}

double IParams::gamma(double energy, int event_class) {
    //return average for now
    return 2.25;
}