/** @file PhotonBinner.cxx
@brief implement class BinnedPhotonData 

$Header$
*/

#include "skymaps/PhotonBinner.h"
#include "skymaps/Band.h"
#include "skymaps/IParams.h"
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

//minimum nside (default 1 is base Healpix pixelization)
unsigned int PhotonBinner::min_nside(1);
void PhotonBinner::set_min_nside(unsigned int new_nside){min_nside = new_nside;}
unsigned int PhotonBinner::get_min_nside(){return min_nside;}

//maximum nside (default 8192 set by 32 bit limit)
unsigned int PhotonBinner::max_nside(8192);
void PhotonBinner::set_max_nside(unsigned int new_nside){max_nside = new_nside;}
unsigned int PhotonBinner::get_max_nside(){return max_nside;}

//sigma scale - constant factor by which to multiply 1/sigma to get nside
//default value gives a pixel side of about sigma/2
double PhotonBinner::m_sigma_scale(2*M_PI/3);
void PhotonBinner::set_sigma_scale(double sigscale){m_sigma_scale = sigscale;}
double PhotonBinner::get_sigma_scale(){return m_sigma_scale;}


PhotonBinner::PhotonBinner(double bins_per_decade)
: m_bins_per_decade(bins_per_decade)
, m_user_nside(false)
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
: m_bins_per_decade(1)
, m_user_nside(false)
{
    // surround the user list with 30 and 1e6 to catch under and overflows
    m_bins.push_back(30);
    std::copy(bins.begin(), bins.end(), std::back_insert_iterator<std::vector<double> >(m_bins));
    m_bins.push_back(infinite);
    setupbins();
}

PhotonBinner::PhotonBinner(const std::vector<double>&edges, const std::vector<int>&f_nside, const std::vector<int>&b_nside)
: m_user_nside(true)
{
    std::copy(edges.begin(), edges.end(), std::back_insert_iterator<std::vector<double> >(m_bins));
    std::copy(f_nside.begin(), f_nside.end(), std::back_insert_iterator<std::vector<int> >(m_fnside));
    std::copy(b_nside.begin(), b_nside.end(), std::back_insert_iterator<std::vector<int> >(m_bnside));
    if (m_bins.front() > 0) {
        m_bins.insert(m_bins.begin(),0);
        m_fnside.insert(m_fnside.begin(),min_nside);
        m_bnside.insert(m_bnside.begin(),min_nside);
    }
    if (m_bins.back() < infinite) {
        m_bins.push_back(infinite);
        m_fnside.push_back(max_nside);
        m_bnside.push_back(max_nside);
    }
}

PhotonBinner::PhotonBinner(double emin, double ratio, int bins)
: m_user_nside(false)
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

    if (m_user_nside) {
        std::vector<double>::const_iterator it = std::lower_bound(m_bins.begin(), m_bins.end(), energy, std::less<double>());
        int nside(event_class==0?m_fnside[it - m_bins.begin() -1]:m_bnside[it - m_bins.begin() -1]);
        return Band(nside,event_class,*(it-1),*(it),0,0);
    }
    // setup for old-style levels, with sigma and gamma for the given energy

	double elow,ehigh,sigma,gamma,nside;
    if( m_bins_per_decade>0){
        // new binning
		// eew: not sure what below line was for...
        // bool high(ehigh==1e6); 
		std::vector<double>::const_iterator it= 
			std::lower_bound(m_bins.begin(), m_bins.end(), energy, std::less<double>());
        elow = *(it-1); ehigh = *(it);
        double ebar = sqrt(elow*ehigh);

        //  set sigma/gamma from CALDB; nside inversion prop to sigma
        sigma = IParams::sigma(ebar,event_class);
        gamma = IParams::gamma(ebar,event_class);
        nside = m_sigma_scale/sigma;
        nside=nside>max_nside?max_nside:nside;
        nside=nside<min_nside?min_nside:nside;
    }else{
        const std::vector<double>& elist ( event_class<=0? s_fenergy : s_benergy );

        std::vector<double>::const_iterator it=
            std::lower_bound(elist.begin(), elist.end(), energy, std::less<double>());
        int level ( it-elist.begin()-1);
        sigma = s_sigma_level[level] * scale_factor(level); 
        gamma = s_gamma_level[level];
        elow  = elist[level];
        ehigh = elist[level+1];
        nside = 1<<level;
	}
    return  Band(nside, event_class, elow, ehigh, sigma, gamma);
}

int PhotonBinner::get_band_key(const astro::Photon& p)
{
	double energy(p.energy());
	if(energy>=infinite) energy=0.9999 * infinite;
	int event_class(p.eventClass());
	if( event_class<0 || event_class >15) event_class = 0;

    if (m_user_nside) {
		//Expect this to be a rare use case, just make the Band to get the key
        std::vector<double>::const_iterator it = std::lower_bound(m_bins.begin(), m_bins.end(), energy, std::less<double>());
        int nside(event_class==0?m_fnside[it - m_bins.begin() -1]:m_bnside[it - m_bins.begin() -1]);
        return int(Band(nside,event_class,*(it-1),*(it),0,0));
    }	
	double elow;
    if( m_bins_per_decade>0){
		std::vector<double>::const_iterator it =
			std::lower_bound(m_bins.begin(), m_bins.end(), energy, std::less<double>());
        elow = *(it-1);
    }else{
        const std::vector<double>& elist ( event_class<=0? s_fenergy : s_benergy );

        std::vector<double>::const_iterator it=
            std::lower_bound(elist.begin(), elist.end(), energy, std::less<double>());
        int level ( it-elist.begin()-1);
		elow = elist[level];
	}
	int key (event_class + 10*static_cast<int>(elow+.5));
	return key;
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


