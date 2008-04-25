/** @file PhotonBinner.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/PhotonBinner.h"
#include "skymaps/Band.h"
#include "astro/Photon.h"
#include <algorithm>
#include <functional>
#include <map>
#include <stdexcept>
#include <iostream>
#include <iomanip>
using namespace skymaps;

namespace {
    //default energy to pixel level conversion
    double front_list[] ={0,0,0,0,0,
        42.555, 100, 235,   552.26, 1297.79 , 3050,  7167.2, 16842.7, 39580., 1e6};
    // note that back is only to level 11, wired-in values above 6500 MeV, ratio of 1.85 below
    double back_list[] ={0,0,0,0,0,
        79.14,  186, 437.1, 1027.19, 2413.9,  6500., 21000, 1e6};

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
std::vector<double> PhotonBinner::s_fenergy(front_list,front_list+sizeof(front_list)/sizeof(double));
std::vector<double> PhotonBinner::s_benergy(back_list,  back_list+sizeof(back_list)/sizeof(double));
std::vector<double> PhotonBinner::s_gamma_level(gamma_list,gamma_list+14);
std::vector<double> PhotonBinner::s_sigma_level(sigma_list,sigma_list+14);


PhotonBinner::PhotonBinner(bool combine):m_comb(combine) 
{
    double fact = 2.35;
    double curfact = 1;
    for(int i(0);i<8;++i) {
        m_bins.push_back(fminE*curfact);
        curfact*=fact;
    }
    setupbins();
}

PhotonBinner::PhotonBinner(std::vector<double>& bins):m_bins(bins),m_comb(false)
{
    setupbins();
}

PhotonBinner::PhotonBinner(double emin, double ratio, int bins):m_comb(false)
{
    double curfact = 1.;
    for(int i(0);i<bins;++i) {
        m_bins.push_back(emin*curfact);
        curfact*=ratio;
    }
    setupbins();
}

void PhotonBinner::add(const astro::Photon& p)
{
    double energy ( p.energy() );
    int event_class (  p.eventClass() );
    if( event_class<0) event_class=0; // should not happen?

    const std::vector<double>& elist = event_class<=0? s_fenergy : s_benergy;
    std::vector<double>::const_iterator it=
        std::lower_bound(elist.begin(), elist.end(), energy, std::less<double>());
    int level ( it-elist.begin()-1);

    unsigned int nside ( 1<<level );
    event_class = 0; // combine front, back

    Band* b = addBand(
        Band(nside, event_class, elist[level], elist[level+1], 
            s_sigma_level[level]*scale_factor(level), s_gamma_level[level])
        ); 

    b->add(p.dir());

}

Band* PhotonBinner::addBand(const Band& band) 
{
    int key(band);
    std::pair<iterator, bool> q = insert(std::make_pair(key,band));
    iterator it = q.first;
    Band & b= it->second;
    return & b;
}



void PhotonBinner::setupbins() {
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

const Band& PhotonBinner::band(int index)const
{
    const_iterator it( find(index) );
    if( it== end()){
        throw std::invalid_argument("PhotonBinner::band -- band id not found");
    }
    return it->second;

}

int PhotonBinner::level(int band, int event_class) const
{
    return (band<0&&abs(event_class)<2)?-1:m_binlevelmap[2*band+event_class];
}



double PhotonBinner::sigma(int level)
{
    return scale_factor(level)*sigma_list[level];
}

double PhotonBinner::gamma(int level)
{
    return gamma_list[level];
}
void PhotonBinner::info(std::ostream& out)const
{
    int total_pixels(0), total_photons(0);
    out << " nside type   emin    emax    sigma   pixels   photons\n";

    for( const_iterator it=begin();  it!=end(); ++it)
    {
        const Band& band = it->second;
        int pixels(band.size()), photons(band.photons());
        out 
            <<std::setw(6) << band.nside()
            <<std::setw(4) << band.event_class()
            <<std::setw(8) << int(band.emin()+0.5)
            <<std::setw(8) << int(band.emax()+0.5)
            <<std::setw(8) << int(band.sigma()*180/M_PI*3600+0.5)
            <<std::setw(10)<< pixels
            <<std::setw(10)<< photons 
            <<std::endl;
        total_photons += photons; total_pixels+=pixels;
    }
    out << " total"
        <<std::setw(38)<<total_pixels
        <<std::setw(10)<<total_photons << std::endl;
}