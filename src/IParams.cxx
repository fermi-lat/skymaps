/** @file IParams.cxx
@brief implement class IParams

$Header$

*/

#include "skymaps/IParams.h"
#include "skymaps/PhotonBinner.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>

using namespace skymaps;

namespace {

    //kluge - all gamma parameterization
    //double f_p[] = {0.01165278, -1.45237162, 1.67211108, 2.28589082};
    //double b_p[] = {0.02614202, -1.55092155, 5.18346847, 3.68094493};

    double elist[] = {0,100,215,464,1000,2154,4642,10000,21544,46416,100000,1000000};

    double fgam[]= {2.25,2.500,2.086,2.132,2.035,1.995,1.860,2.108,2.108,2.141,2.25,2.25};
    double bgam[]= {2.25,2.243,2.097,2.073,2.192,2.100,2.049,1.983,2.250,2.250,2.25,2.25};
    //first light -> 07/01/2008-07-27/2008
    double f_p[] = {0.0, -1.359, 1.638, 0.189};
    double b_p[] = {0.043, -1.575, 6.476, 2.602};
}


std::vector<double> IParams::s_elist(elist,elist+12);
std::vector<double> IParams::s_fparams(f_p,f_p+4);
std::vector<double> IParams::s_bparams(b_p,b_p+4);
void IParams::set_fp(std::vector<double> params){s_fparams=params;}
void IParams::set_bp(std::vector<double> params){s_bparams=params;}

IParams::IParams(std::string& name)
{
    //TODO: setup database, for now use defaults from allGamma v14r8 
}

double IParams::sigma(double energy, int event_class){
    double deg = 0;
    if(event_class) {
        //from parameterization of PSF function
        deg = sqrt(s_bparams[0]*s_bparams[0]+s_bparams[2]*pow(energy/100,IParams::s_bparams[1])+s_bparams[3]*pow(energy/100,2*s_bparams[1]));
    } else {
        deg = sqrt(s_fparams[0]*s_fparams[0]+s_fparams[2]*pow(energy/100,IParams::s_fparams[1])+s_fparams[3]*pow(energy/100,2*s_fparams[1]));
    }
    return deg*M_PI/180;
}

double IParams::gamma(double energy, int event_class) {
    const std::vector<double>& elisti (s_elist);

    std::vector<double>::const_iterator itold=
            std::lower_bound(elisti.begin(), elisti.end(), energy, std::less<double>());
    int level = itold-elisti.begin()-1;

    //return average for now
    return event_class==0?fgam[level]:bgam[level];
    //return 2.25;
}

std::vector<double> IParams::params(int event_class) 
{
    return event_class?s_bparams:s_fparams;
}
