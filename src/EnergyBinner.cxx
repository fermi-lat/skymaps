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
        42.5, 100, 235, 552.26, 1297.79 ,3050, 7167.2, 16842.7, 39580.1};

    double back_list[] ={0,0,0,0,0,
        79.14, 186, 437.1, 1027.19, 2413.9, 5672.8, 13330.8};
}

//static list of energy to pixel level conversions
std::vector<double> EnergyBinner::s_fenergy(front_list,front_list+14);
std::vector<double> EnergyBinner::s_benergy(back_list,back_list+12);


EnergyBinner::EnergyBinner(bool combine):m_comb(combine) 
{
    double eminf = front_list[5];
    double eminb = back_list[5];
    double fact = 2.35;
    double curfact = 1;
    for(int i(0);i<9;++i) {
        m_bins.push_back(eminf*curfact);
        m_bins.push_back(eminb*curfact);
        curfact*=fact;
    }
    setupbins();
}
 
EnergyBinner::EnergyBinner(std::vector<double>& bins):m_bins(bins),m_comb(false)
{
    setupbins();
}

EnergyBinner::EnergyBinner(double eminf, double emaxf, double eminb, double emaxb, int bins):m_comb(false)
{
    double flogfact = pow(emaxf/eminf,1./bins);
    double fcurfact = 1.;
    double blogfact = pow(emaxb/eminb,1./bins);
    double bcurfact = 1.;
    for(int i(0);i<bins;++i) {
        m_bins.push_back(eminf*fcurfact);
        m_bins.push_back(eminb*bcurfact);
        fcurfact*=flogfact;
        bcurfact*=blogfact;
    }
    setupbins();
}

void EnergyBinner::setupbins() {
    //clear current mapping
    m_binlevelmap.clear();
    for(int i(0);i<m_bins.size();++i) {
        int flag(0);
        //even i->front, odd i->back
        if(i%2==0) {
            int j(0);
            for(;(j<s_fenergy.size())&&(flag==0);++j) {
                if(s_fenergy[j]>=m_bins[i]) {
                    m_binlevelmap.push_back(j);
                    flag = 1;
                }
            }
            if(flag==0) m_binlevelmap.push_back(j-1);
        } else {
            int j(0);
            for(;(j<s_benergy.size())&&(flag==0);++j) {
                if(s_benergy[j]>=m_bins[i]) {
                    m_binlevelmap.push_back(j);
                    flag = 1;
                }
            }
            if(flag==0) m_binlevelmap.push_back(j-1);
        }
    }
}

int EnergyBinner::band(const astro::Photon& p) {
    int offset = p.eventClass();
    double energy = p.energy();
    //return -1 if below minimum energy
    if(energy<m_bins[offset]) return -1;
    //only traverse through front or back bins
    int i(0);
    for(;i<m_bins.size()/2;++i) {
        if(energy<m_bins[2*i+offset]) return m_comb?i-1:(2*(i-1)+offset);
    }
    return m_comb?i-1:2*(i-1)+offset;
}

int EnergyBinner::level(const astro::Photon& p) 
{
    return level(m_comb?2*band(p)+p.eventClass():band(p));
}