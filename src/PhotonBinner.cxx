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
#include <iterator>
using namespace skymaps;

namespace {
    //default energy to pixel level conversion
    double front_list[] ={0,0,0,0,10,
        42.555, 100, 235,   552.26, 1297.79 , 3050,  7167.2, 16842.7, 39580., 1e6};
    // note that back is only to level 11, wired-in values above 6500 MeV, ratio of 1.85 below
    double back_list[] ={0,0,0,0,10,
        79.14,  186, 437.1, 1027.19, 2413.9,  6500., 21000, 1e6};

    int base_level = 6;

    double fminE = 42.553;
    double bminE = 186;

    double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}

    double gamma_list[] ={0,0,0,0, 2.2,
        2.25,  2.27,  2.22,  2.31,  2.30,  2.31,  2.16,  2.19,  2.07};
    double sigma_list[] ={0,0,0,0,0.5,
         0.343, 0.335, 0.319, 0.431,  0.449,  0.499,  0.566,  0.698,  0.818
        };
    double infinite(1e6); // largest energy
    double sqr(double x){return x*x;}

}

//static list of energy to pixel level conversions
std::vector<double> PhotonBinner::s_fenergy(front_list,front_list+sizeof(front_list)/sizeof(double));
std::vector<double> PhotonBinner::s_benergy(back_list,  back_list+sizeof(back_list)/sizeof(double));
std::vector<double> PhotonBinner::s_gamma_level(gamma_list,gamma_list+14);
std::vector<double> PhotonBinner::s_sigma_level(sigma_list,sigma_list+14);


PhotonBinner::PhotonBinner(double bins_per_decade)
: m_bins_per_decade(bins_per_decade)
{

    double eratio( 2.35);
    if( bins_per_decade<=0){

    }else{
        eratio = pow(10., 1./bins_per_decade);
        double emin(10.), emax(1e5); // energy range
        m_bins.push_back(0);
        for( double e(emin); e< 1.01*emax; e*=eratio){
            m_bins.push_back(e);
        }
        m_bins.push_back(infinite);
    }
}

PhotonBinner::PhotonBinner(const std::vector<double>& bins)
{
    // surround the user list with 30 and 1e6 to catch under and overflows
    m_bins.push_back(30);
    std::copy(bins.begin(), bins.end(), std::back_insert_iterator<std::vector<double> >(m_bins));
    m_bins.push_back(1e6); 
    setupbins();
}

PhotonBinner::PhotonBinner(double emin, double ratio, int bins)
{
    double curfact = 1.;
    for(int i(0);i<bins;++i) {
        m_bins.push_back(emin*curfact);
        curfact*=ratio;
    }
    setupbins();
}

skymaps::Band PhotonBinner::operator()(const astro::Photon& p)const
{
    double energy ( p.energy() );
    if( energy>=infinite) energy=0.9999 * infinite;

    int event_class (  p.eventClass() );
    if( event_class<0 || event_class>15) event_class=0; // should not happen?

    // setup for old-style levels, with sigma and gamma for the given energy
    const std::vector<double>& elist ( event_class<=0? s_fenergy : s_benergy );

    std::vector<double>::const_iterator it=
        std::lower_bound(elist.begin(), elist.end(), energy, std::less<double>());
    int level ( it-elist.begin()-1);
    double sigma( s_sigma_level[level] * scale_factor(level) )
        ,  gamma( s_gamma_level[level] )
        ,  elow( elist[level] )
        ,  ehigh( elist[level+1] );
    unsigned int nside ( 1<<level );

    if( m_bins_per_decade>0){
        // no, new binning
        bool high(ehigh==1e6); 
        double oldebar(sqrt(elow*ehigh));
        it= std::lower_bound(m_bins.begin(), m_bins.end(), energy, std::less<double>());
        elow = *(it-1); ehigh = *(it);
        double ebar(sqrt(elow*ehigh));

        // get old parameters for this energy
        std::vector<double>::const_iterator itold=
            std::lower_bound(elist.begin(), elist.end(), ebar, std::less<double>());
        level = itold-elist.begin()-1;
        sigma = s_sigma_level[level] * scale_factor(level);
        gamma = s_gamma_level[level] ;
        nside = 1<<level;
        nside = 3*nside/2;  // double number of pixels

        //  kluge, from Marshall's fits
        if( event_class==0){
            sigma = sqrt( sqr(1.62*pow(ebar/100, -0.8)) + sqr(0.011) )*M_PI/180;
        }else if(event_class==1){
            sigma = sqrt( sqr(2.34*pow(ebar/100, -0.8)) + sqr(0.024) )*M_PI/180;
        }
    }
    return  Band(nside, event_class, elow, ehigh, sigma, gamma);
}


void PhotonBinner::setupbins() {
#if 1

#else
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
#endif
}


