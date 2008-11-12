/** @file IParams.cxx
@brief implement class IParams

$Header$

*/

#include "skymaps/IParams.h"
#include "skymaps/PhotonBinner.h"
#include "embed_python/Module.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>

using namespace skymaps;

namespace {

    //default - all gamma parameterization v15r0 
    double f_p[] = {0.01165278, -1.45237162, 1.67211108, 2.28589082};
    double b_p[] = {0.02614202, -1.55092155, 5.18346847, 3.68094493};

    double elist[] = {0   , 100, 178, 316, 562,1000,1778,3162,5623,10000,17783,31623,56234,100000,1000000};

    double fgam[]=   {2.65,2.65,2.13,1.97,2.08,2.20,2.43,2.51,2.50, 2.17, 1.94, 1.73, 1.66,  1.62};
    double bgam[]=   {2.59,2.59,2.09,1.90,1.90,1.64,1.63,1.72,1.76, 1.78, 1.68, 1.71, 1.73,  1.78};

}


std::vector<double> IParams::s_elist(elist,elist+15);
std::vector<double> IParams::s_fparams(f_p,f_p+4);
std::vector<double> IParams::s_bparams(b_p,b_p+4);
std::vector<double> IParams::s_fgam(fgam,fgam+15);
std::vector<double> IParams::s_bgam(bgam,bgam+15);
bool IParams::s_init(false);

void IParams::set_elist(std::vector<double> elist) {s_elist=elist;}
void IParams::set_fp(std::vector<double> params){s_fparams=params;}
void IParams::set_bp(std::vector<double> params){s_bparams=params;}
void IParams::set_fgam(std::vector<double> gams) {s_fgam=gams;}
void IParams::set_bgam(std::vector<double> gams) {s_bgam=gams;}

IParams::IParams(std::string& name)
{
    //TODO: setup database, for now use defaults from allGamma v14r8 
}

double IParams::sigma(double energy, int event_class){
    //check to see if initialized, if not do it
    if(!s_init)  init();
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
    //check to see if initialized, if not do it
    if(!s_init)  init();

    const std::vector<double>& elisti (s_elist);

    //find nearest energy bin
    std::vector<double>::const_iterator itold=
            std::lower_bound(elisti.begin(), elisti.end(), energy, std::less<double>());
    int level = itold-elisti.begin()-1;

    //TODO: implement weighting function
    return event_class==0?s_fgam[level]:s_bgam[level];
}

std::vector<double> IParams::params(int event_class) 
{
    return event_class?s_bparams:s_fparams;
}

void IParams::init() {
    s_init=true;
    int argc=0;
    char** argv=((char**)0); // avoid warning message
    std::string python_path("../python");
    try {
    embed_python::Module setup(python_path , "psf_defaults",  argc, argv);
    std::vector<double> ps;
    //get energy list from CALDB file
    setup.getList("energy",ps);
    set_elist(ps);
    ps.clear();
    //get front sigma paramterization
    setup.getList("fparams",ps);
    set_fp(ps);
    ps.clear();
    //get back sigma parameterization
    setup.getList("bparams",ps);
    set_bp(ps);
    ps.clear();
    //now, just get the gamma values
    setup.getList("fgam",ps);
    set_fgam(ps);
    ps.clear();
    setup.getList("bgam",ps);
    set_bgam(ps);
    ps.clear();
    } catch(const std::exception& e){
        std::cout << "Caught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        std::cout << "Using default PSF values" << std::endl;
    }
}