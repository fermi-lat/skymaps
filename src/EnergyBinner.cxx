/** @file EnergyBinner.cxx
@brief implement EnergyBinner
$Header$
*/

#include "skymaps/EnergyBinner.h"
#include "astro/Photon.h"

using namespace skymaps;

namespace {
    //default energy to pixel level conversion
    double front_list[] ={0,0,0,0,0,
        42.555, 100, 235, 552.26, 1297.79 ,3050, 7167.2, 16842.7, 39580.1};

    double back_list[] ={0,0,0,0,0,
        79.14, 186, 437.1, 1027.19, 2413.9, 5672.8, 13330.8};

    int base_level = 6;

    double fminE = 42.553;
    double bminE = 186;

    double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}

    double gamma_list[] ={0,0,0,0,0,
        2.25,  2.27,  2.22,  2.31,  2.30,  2.31,  2.16,  2.19,  2.07};
    double sigma_list[] ={0,0,0,0,0,
        0.343, 0.4199,0.4249 ,0.4202 ,0.4028 ,0.4223 ,0.4438 ,0.5113 ,0.5596 };
}

//static list of energy to pixel level conversions
std::vector<double> EnergyBinner::s_fenergy(front_list,front_list+14);
std::vector<double> EnergyBinner::s_benergy(back_list,back_list+12);
std::vector<double> EnergyBinner::s_gamma_level(gamma_list,gamma_list+14);
std::vector<double> EnergyBinner::s_sigma_level(sigma_list,sigma_list+14);
EnergyBinner* EnergyBinner::s_instance=0; 

EnergyBinner* EnergyBinner::Instance(bool combine)
{
    if(s_instance==0){
        s_instance = new EnergyBinner(combine);
    }
    return s_instance;
}

EnergyBinner* EnergyBinner::Instance(std::vector<double> &bins)
{
    if(s_instance==0){
        s_instance = new EnergyBinner(bins);
    }
    return s_instance;
}

EnergyBinner* EnergyBinner::Instance(double emin, double ratio, int bins)
{
    if(s_instance==0){
        s_instance = new EnergyBinner(emin, ratio, bins);
    }
    return s_instance;
}

EnergyBinner::EnergyBinner(bool combine):m_comb(combine) 
{
    double fact = 2.35;
    double curfact = 1;
    for(int i(0);i<8;++i) {
        m_bins.push_back(fminE*curfact);
        curfact*=fact;
    }
    setupbins();
}

EnergyBinner::EnergyBinner(std::vector<double>& bins):m_bins(bins),m_comb(false)
{
    setupbins();
}

EnergyBinner::EnergyBinner(double emin, double ratio, int bins):m_comb(false)
{
    double curfact = 1.;
    for(int i(0);i<bins;++i) {
        m_bins.push_back(emin*curfact);
        curfact*=ratio;
    }
    setupbins();
}

void EnergyBinner::setupbins() {
    //clear current mapping
    m_binlevelmap.clear();
    for(int i(0);i<m_bins.size();++i) {
        int flag(0);
        //front
        int j(0);
        for(;(j<s_fenergy.size())&&(flag==0);++j) {
            if(s_fenergy[j]>=m_bins[i]) {
                m_binlevelmap.push_back(j);
                flag = 1;
            }
        }
        if(flag==0) m_binlevelmap.push_back(j-1);
        //odd i->back
        j=0;
        flag=0;
        for(;(j<s_benergy.size())&&(flag==0);++j) {
            if(s_benergy[j]>=m_bins[i]) {
                m_binlevelmap.push_back(j);
                flag = 1;
            }
        }
        if(flag==0) m_binlevelmap.push_back(j-1);
    }
}

int EnergyBinner::band(const astro::Photon& p) const{
    int offset = p.eventClass();
    double energy = p.energy();
    //return -1 if below minimum energy
    if(energy<m_bins[0]) return -1;
    //only traverse through front or back bins
    int i(0);
    for(;i<m_bins.size();++i) {
        if(energy<m_bins[i]) return i-1;
    }
    return i-1;
}

int EnergyBinner::level(const astro::Photon& p) const
{
    return level(band(p),p.eventClass());
}

double EnergyBinner::sigma(int level)
{
    return scale_factor(level)*sigma_list[level];
}

double EnergyBinner::gamma(int level)
{
    return gamma_list[level];
}