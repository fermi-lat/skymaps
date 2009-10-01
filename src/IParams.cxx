/** @file IParams.cxx
@brief implement class IParams

$Header$

*/

#include "skymaps/IParams.h"
#include "skymaps/ExposureWeighter.h"
#include "skymaps/PhotonBinner.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

using namespace skymaps;

namespace {

    //default - all gamma parameterization v15r0 
    double f_p[] = {0.01165278, -1.45237162, 1.67211108, 2.28589082};
    double b_p[] = {0.02614202, -1.55092155, 5.18346847, 3.68094493};

    double elist[] = {1   , 100, 178, 316, 562,1000,1778,3162,5623,10000,17783,31623,56234,100000,1000000};

    double fsig[]=   {75.088,1.423,0.868,0.532,0.329,0.208,0.139,0.103,0.085,0.077,0.074,0.073,0.073,0.073,0.073};
    double bsig[]=   {617.153,3.896,2.07,1.107,0.599,0.338,0.214,0.162,0.146,0.139,0.137,0.137,0.137,0.137,0.137};
    double fgam[]=   {2.65,2.65,2.13,1.97,2.08,2.20,2.43,2.51,2.50, 2.17, 1.94, 1.73, 1.66,  1.62, 1.62};
    double bgam[]=   {2.59,2.59,2.09,1.90,1.90,1.64,1.63,1.72,1.76, 1.78, 1.68, 1.71, 1.73,  1.78, 1.78};

    double mf(1.5);
    double mb(1.2);

    bool verbose(false);
}


std::vector<double> IParams::s_elist(elist,elist+15);
std::vector<double> IParams::s_fparams(f_p,f_p+4);
std::vector<double> IParams::s_bparams(b_p,b_p+4);
std::vector<double> IParams::s_fsig(fsig,fsig+15);
std::vector<double> IParams::s_bsig(bsig,bsig+15);
std::vector<double> IParams::s_fgam(fgam,fgam+15);
std::vector<double> IParams::s_bgam(bgam,bgam+15);
std::vector<double> IParams::s_fwght(0);
std::vector<double> IParams::s_bwght(0);
std::string         IParams::s_CALDB("");
std::string         IParams::s_livetimefile("");
bool IParams::s_init(false);
astro::SkyDir s_dir;

void IParams::set_elist(std::vector<double> elist) {s_elist=elist;}
void IParams::set_fp(std::vector<double> params){s_fparams=params;}
void IParams::set_bp(std::vector<double> params){s_bparams=params;}
void IParams::set_fsig(std::vector<double> sigs) {s_fsig=sigs;}
void IParams::set_bsig(std::vector<double> sigs) {s_bsig=sigs;}
void IParams::set_fgam(std::vector<double> gams) {s_fgam=gams;}
void IParams::set_bgam(std::vector<double> gams) {s_bgam=gams;}
void IParams::set_fwght(std::vector<double> wght) {s_fwght=wght;}
void IParams::set_bwght(std::vector<double> wght) {s_bwght=wght;}
void IParams::set_CALDB(const std::string& dir) {s_CALDB=dir;}
void IParams::set_livetimefile(const std::string& ltfile) {s_livetimefile = ltfile;}
void IParams::set_skydir(const astro::SkyDir& dir) {s_dir = dir;}







IParams::IParams()
{
    //TODO: setup database, for now use defaults from allGamma v14r8 
}


double IParams::sigma(double energy, int event_class){

    //check to see if initialized, if not do it
    //if(!s_init)  init();

    const std::vector<double>& elisti (s_elist);

    //find nearest energy bin
    std::vector<double>::const_iterator itold=
        std::lower_bound(elisti.begin(), elisti.end(), energy, std::less<double>());
    int level = itold-elisti.begin()-1;
    if (level==elisti.size()-1) level--;
    if (energy <= s_elist[0]) { level = 0; } // Kerr - does power law extrapolation at low energies
    double emin(s_elist[level]),emax(s_elist[level+1]);
    emin>0?emin:emin=1;
    double smin(event_class==0?s_fsig[level]:s_bsig[level]),smax(event_class==0?s_fsig[level+1]:s_bsig[level+1]);
    
    //weighting function
    double m = (log(smax)-log(smin))/(log(emax)-log(emin));
    double sbar = exp(m*(log(energy)-log(emax))+log(smax));
    return sbar*M_PI/180;
}

double IParams::gamma(double energy, int event_class) {
    //check to see if initialized, if not do it
    //if(!s_init)  init();

    if (energy <= s_elist[0]) { // Kerr -- handle "underflow" by returning parameters for low energy
        return event_class==0 ? s_fgam[0] : s_bgam[0];
    }

    const std::vector<double>& elisti (s_elist);

    //find nearest energy bin
    std::vector<double>::const_iterator itold=
        std::lower_bound(elisti.begin(), elisti.end(), energy, std::less<double>());
    int level = itold-elisti.begin()-1;
    if (level==elisti.size()-1) level--;
    double emin(s_elist[level]),emax(s_elist[level+1]);
    emin>0?emin:emin=1;
    double gmin(event_class==0?s_fgam[level]:s_bgam[level]),gmax(event_class==0?s_fgam[level+1]:s_bgam[level+1]);
    //weighting function
    double m = (gmax-gmin)/(log(emax)-log(emin));
    double gbar = m*(log(energy)-log(emax))+gmax;
    return gbar;
}

double IParams::effaWeight(double energy, int event_class){

    //check to see if initialized, if not do it
    //if(!s_init)  init();

    const std::vector<double>& elisti (s_elist);

    //find nearest energy bin
    std::vector<double>::const_iterator itold=
        std::lower_bound(elisti.begin(), elisti.end(), energy, std::less<double>());
    int level = itold-elisti.begin()-1;
    if (level==elisti.size()-1) level--;
    if (energy <= s_elist[0]) { level = 0; } // Kerr - does power law extrapolation at low energies
    double emin(s_elist[level]),emax(s_elist[level+1]);
    emin>0?emin:emin=1;
    double wmin(event_class==0?s_fwght[level]:s_bwght[level]),wmax(event_class==0?s_fwght[level+1]:s_bwght[level+1]);
    
    //weighting function
    double m = (log(wmax)-log(wmin))/(log(emax)-log(emin));
    double wbar = exp(m*(log(energy)-log(emax))+log(wmax));
    return wbar*M_PI/180;
}

std::vector<double> IParams::params(int event_class) 
{
    return event_class?s_bparams:s_fparams;
}

double IParams::scale(double energy, int event_class) {
    double p0(0),p1(0);
    if(event_class) {
        p0=0.096;
        p1=0.00130;
    }else {
        p0=0.058;
        p1=0.000377;
    }
    return sqrt(p0*p0*pow(energy/100.,-1.6)+p1*p1)*180/M_PI;
}

void IParams::init(const std::string& name, const std::string& clevel, const std::string& file) {
    s_init=true;

    //Seems to crash python, so the routine to read fits has been
    //changed to C++ implementation

    const tip::Table * ptablef(0);
    const tip::Table * ptableb(0);
    std::string psf_table("RPSF");

    if( file.empty() ){
        // no overriding filename
        
        if( s_CALDB.empty()){
            const char* c(::getenv("CALDB") );
            if( c==0){
                std::cerr << "IParams:: CALDB is not set" << std::endl;
                throw std::invalid_argument("IParams:: CALDB is not set");
            }
            s_CALDB = std::string(c);
        }
    }

    std::string fcaldb(s_CALDB+"/bcf/psf/psf_"+name+"_"+clevel+"_front.fits");
    std::string bcaldb(s_CALDB+"/bcf/psf/psf_"+name+"_"+clevel+"_back.fits");

    try{
        ptablef = tip::IFileSvc::instance().readTable(fcaldb,psf_table);
    }catch(const std::exception&){
        std::cerr << "IParams: init() "<< psf_table << " not found in file " << fcaldb << std::endl;
        throw;
    }
    if( verbose) std::cout << "Using front PSFs: " << fcaldb << std::endl;
    try{
        ptableb = tip::IFileSvc::instance().readTable(bcaldb,psf_table);
    }catch(const std::exception&){
        std::cerr << "IParams: init() "<< psf_table << " not found in file " << bcaldb << std::endl;
        throw;
    }
    if( verbose) std::cout << "Using back PSFs: " << bcaldb << std::endl;
    const tip::Table& ftable(*ptablef);  // reference for convenience
    const tip::Table& btable(*ptableb);
    tip::Table::ConstIterator itor = ftable.begin();
    std::vector<double> energy_lo,energy_hi,fsigmas,fgammas,bsigmas,bgammas,cost_lo,cost_hi;
    std::vector<double> en,fs,fg,fw,bs,bg,bw;

    (*itor)["ENERG_LO"].get(energy_lo);
    (*itor)["ENERG_HI"].get(energy_hi);
    (*itor)["SIGMA"].get(fsigmas);
    (*itor)["GCORE"].get(fgammas);
    (*itor)["CTHETA_LO"].get(cost_lo);
    (*itor)["CTHETA_HI"].get(cost_hi);

    itor = btable.begin();
    (*itor)["SIGMA"].get(bsigmas);
    (*itor)["GCORE"].get(bgammas);
    
    delete ptablef;
    delete ptableb;

    // Build something to weight the PSF in average; if static variable giving a
    // livetime cube file is set, uses full exposure -- n.b. set the static SkyDir
    // as well or else get the default location!
    std::string faeff_str(s_CALDB+"/bcf/ea/aeff_"+name+"_"+clevel+"_front.fits");
    std::string baeff_str(s_CALDB+"/bcf/ea/aeff_"+name+"_"+clevel+"_back.fits");
    ExposureWeighter ew(faeff_str,baeff_str,s_livetimefile);

    //iterate through energies
    for(unsigned int i(0);i<energy_lo.size();++i) {
        double w0(0),w1(0),s0(0),s1(0),g0(0),g1(0);
        double enr( sqrt(energy_lo[i]*energy_hi[i]) );

        //iterate through cos-th bins
        for(unsigned int j(0);j<cost_lo.size();++j) {
            
            double wf ( ew(cost_lo[j],cost_hi[j],energy_lo[i],energy_hi[i],0,s_dir) );
            double wb ( ew(cost_lo[j],cost_hi[j],energy_lo[i],energy_hi[i],1,s_dir) );

            wf = wf > 0 ? wf : 1e-20; //guard against 0 effective area
            wb = wb > 0 ? wb : 1e-20;

            w0+=wf;
            w1+=wb;
            s0+=wf*fsigmas[energy_lo.size()*j+i];
            s1+=wb*bsigmas[energy_lo.size()*j+i];
            g0+=wf*fgammas[energy_lo.size()*j+i];
            g1+=wb*bgammas[energy_lo.size()*j+i];
        }

        en.push_back(enr);
        s0/=w0;
        s1/=w1;
        g0/=w0;
        g1/=w1;
        s0*=scale(enr,0);
        s1*=scale(enr,1);
        fs.push_back(s0);
        bs.push_back(s1);
        fg.push_back(g0);
        bg.push_back(g1);
	fw.push_back(w0);
	bw.push_back(w1);
    }
    //setup parameters
    set_fsig(fs);
    set_bsig(bs);
    set_fgam(fg);
    set_bgam(bg);
    set_fwght(fw);
    set_bwght(bw);
    set_elist(en);
}



void IParams::print_parameters() {
    std::cout << "Energy     F_sig    B_sig    F_gam    B_gam" << std::endl;
    for (int j(0); j < s_elist.size(); ++j) {
        std::cout << s_elist[j] << "    " << s_fsig[j] << "    " << s_bsig[j] << "    " << s_fgam[j] << "    " << s_bgam[j] << std::endl;
    }
}
