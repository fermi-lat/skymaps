/** @file Convolution.cxx
@brief Convolves healpix maps with a sky function

@author M. Roth 

$Header$
*/
#include "skymaps/Convolution.h"
#include "healpix/AlmOp.h"
#include "healpix/HealPixel.h"
#include "skymaps/PsfFunction.h"
#include "skymaps/IParams.h"

// EAC, update include path
#include "healpix/alm_filter_tools.h"

#include <map>

#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>

using namespace healpix;

namespace{
    void ShowPercent(int sofar, int total, int found)
    {
        static int toskip(50), skipped(0);
        if(++skipped<toskip) return; skipped=0;
        static int lastpercent(-1);
        int percent( static_cast<int>(100 * sofar / total +0.5) );
        if( percent==lastpercent) return;
        lastpercent=percent;
        char   s[50];
        sprintf(s, "%d%%,  %d found.", percent, found);
        std::cout << s;
        if (sofar < total)
        {
            for (size_t j = 0; j < strlen(s); ++j)
                std::cout << "\b";
        }
        else
            std::cout << std::endl;
    }

    //klugy scale factor to visualize pmap
    double klug = 1;
    int minE = 555;
}

using namespace skymaps;

Convolution::Convolution (const skymaps::SkySpectrum &sf, const skymaps::SkySpectrum &ker, double energy, int level) 
: m_level(level)
, m_sf(&sf)
, m_ker(&ker)
{
    createConv(sf,ker,energy);
}

Convolution::Convolution(const skymaps::SkySpectrum &sf, double energy, int level) 
: m_level(level)
, m_sf(&sf)
, m_ker(0)
{
    createConv(sf,energy);
}

void Convolution::createConv(const skymaps::SkySpectrum& sf, double energy) {
    iterator it = this->find((int)energy);
    if(it==end()) {
        std::cout << "*********************************************************" << std::endl;
        std::cout << "     Convolution: level " << m_level << " and energy " << energy << std::endl << std::endl;
        skymaps::PsfFunction psf(2.25);
        std::cout << "          Allocating map and harmonic storage...";
        Map<double> sfm(m_level);
        Map<double> psfm(m_level);
        double sigmasq = IParams::sigma(energy,0)*IParams::sigma(energy,0)+IParams::sigma(energy,1)*IParams::sigma(energy,1);
        sigmasq *= (M_PI*M_PI)/(180*180);
        //point spread function parameter for energy band
        //Nyquist frequency for spherical harmonics is half the number of discrete latitudes
        AlmOp<xcomplex<double> > sfh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
        AlmOp<xcomplex<double> > psfh(2*psfm.map()->Nside(),2*psfm.map()->Nside());
        std::cout << "done!" << std::endl;
        std::cout << "          Populating sky and psf maps...";
        //put psf at north pole (azimuthal symmetry)
        //iterate over each latitudnal ring
        for(int i = 1;i<4*sfm.map()->Nside();i++) {
            int startpix,ringpix;
            double theta;
            bool shift;
            sfm.map()->get_ring_info2(i,startpix,ringpix,theta,shift);
            double u = 0.5*(4*tan(theta/2)*tan(theta/2))/(sigmasq);
            u = psf(u);
            for(int j=0;j<ringpix;++j) {
                ShowPercent(j+startpix,sfm.map()->Npix(),j+startpix);
                HealPixel hp(sfm.map()->ring2nest(j+startpix),m_level);
                (*psfm.map())[startpix+j] = u;
                (*sfm.map())[startpix+j] = sf.value(hp,energy)*klug;
            }
        }
        std::cout << "done!" << std::endl;
        std::cout << "          Calculating sky harmonics...";
        map2alm_iter(*sfm.map(),*sfh.Alms(),0);
        std::cout << "done!" << std::endl;
        std::cout << "          Calculating psf harmonics...";
        map2alm_iter(*psfm.map(),*psfh.Alms(),0);
        std::cout << "done!" << std::endl;
        mf_constantnoise(*sfh.Alms(),*psfh.Alms());
        std::cout << "          Performing convolution...";
        alm2map(*sfh.Alms(),*sfm.map());
        std::cout << "done!" << std::endl;
        this->insert(std::map<int,healpix::Map<double> >::value_type((int)energy,sfm));
    }
}

void Convolution::createConv(const skymaps::SkySpectrum& sf, const skymaps::SkySpectrum& ker, double energy) {
    iterator it = this->find((int)energy);
    if(it==end()) {
        std::cout << "*********************************************************" << std::endl;
        std::cout << "     Convolution: level " << m_level << " and energy " << energy << std::endl << std::endl;
        std::cout << "          Allocating map and harmonic storage...";
        Map<double> sfm(m_level);
        Map<double> kerm(m_level);

        AlmOp<xcomplex<double> > sfh(2*sfm.map()->Nside(),2*sfm.map()->Nside());
        AlmOp<xcomplex<double> > kerh(2*kerm.map()->Nside(),2*kerm.map()->Nside());
        std::cout << "done!" << std::endl;
        std::cout << "          Populating kernel and sky function maps...";
        int idx(0);
        int npix(sfm.map()->Npix());
        for(int i(0);i<npix;++i) {
            ShowPercent(i,npix,i);
            HealPixel hp(i,m_level);
            idx = (sfm.map()->nest2ring(i));

            (*sfm.map())[idx]  = sf(hp,energy);
            (*kerm.map())[idx] = ker(hp,energy);
        }

        // map2alm_iter(*sfm.map(),*sfh.Alms(),0); // Kerr -- is this redundant?
        // map2alm_iter(*kerm.map(),*kerh.Alms(),0);
        std::cout << "done!" << std::endl;
        // Get normalizations
        double sf_norm = sfm.map()->average()*sfm.map()->Npix();
        double ker_norm = kerm.map()->average()*kerm.map()->Npix();
        std::cout << "             map integral: " << sf_norm << std::endl;
        std::cout << "             ker integral: " << ker_norm << std::endl;
        std::cout << "          Calculating sky harmonics...";
        map2alm_iter(*sfm.map(),*sfh.Alms(),0);
        std::cout << "done!" << std::endl;
        std::cout << "          Calculating psf harmonics...";
        map2alm_iter(*kerm.map(),*kerh.Alms(),0);
        std::cout << "done!" << std::endl;
        mf_constantnoise(*sfh.Alms(),*kerh.Alms());
        std::cout << "          Performing convolution...";
        alm2map(*sfh.Alms(),*sfm.map());
        std::cout << "done!" << std::endl;
        // Check normalization
        double sf_norm_pc = sfm.map()->average()*sfm.map()->Npix();
        double scale(sf_norm/sf_norm_pc);
        std::cout << "             map integral: " << sf_norm_pc << " percent difference: " << (sf_norm_pc-sf_norm)/sf_norm*100 << std::endl;
        std::cout << "          Renormalizing...";
        for (int j = 0; j < sfm.map()->Npix(); ++j) {
            double& pval = sfm.map()->operator[](j);
            pval *= scale;
        }
        double sf_norm_renorm = sfm.map()->average()*sfm.map()->Npix();
        std:: cout << "done!" << std::endl;
        std::cout << "             map integral: " << sf_norm_renorm << " percent difference: " << (sf_norm_renorm-sf_norm)/sf_norm*100 << std::endl;
        this->insert(std::map<int,healpix::Map<double> >::value_type((int)energy,sfm));
    }
}

double Convolution::value(const astro::SkyDir& dir,double e) const{
    if((this)->size()<2) return value(dir);
    int e1(layer(e,true));//, e2(layer(e,false));
    int pixel_index = find(e1)->second.cmap()->nest2ring(healpix::HealPixel(dir,m_level).index());
    /*double f1((*find(e1)->second.cmap())[pixel_index]),f2((*find(e2)->second.cmap())[pixel_index]);
    if(f1<=0&&f2<=0) return 0;
    f1<=0?f1=1e-40:0;
    f2<=0?f2=1e-40:0;
    double alpha ( log(f1/f2)/log((double)(e2)/e1) );
    return f1*pow(e/e1,alpha);*/
    return (*find(e1)->second.cmap())[pixel_index];
}

double Convolution::integral(const astro::SkyDir& dir, double a, double b) const{
    double fa(value(dir, a)), fb(value(dir,b));
    double q ( 1. - log(fa/fb)/log(b/a) );
    return fa* a * (pow(b/a, q)-1)/q;
}

int Convolution::layer(double e, bool front) const{
    const_iterator it = begin();
    const_iterator stop = end();
    int laste = it->first;
    for(; it!=stop; ++it){
        if( it->first>e ) break;
    }
    if(it==stop) {
        --it;
        --it;
    }
    return front?it->first:(++it)->first;
}

double Convolution::value(const astro::SkyDir &dir) const {
    
    const_iterator it=begin();
    int pix_index=it->second.cmap()->nest2ring(healpix::HealPixel(dir,m_level).index());
    double sum((*it->second.cmap())[pix_index]);
    ++it;
    while(it!=end())
    {
        sum+=(*it->second.cmap())[pix_index];
        ++it;
    }
    return sum;
}

std::string Convolution::name()const
{
    std::stringstream label;
    label << "Convolution of " << m_sf->name() << " and " <<( m_ker==0? "PSF" : m_ker->name());
    
    return label.str();
}


