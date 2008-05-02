/** @file likelihood_main.cxx

$Header$

*/
//#include "skymaps/PhotonMap.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/Exposure.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/PhotonBinner.h"
#include "astro/Photon.h"
#include "astro/SkyDir.h"


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

using astro::Photon;
using astro::SkyDir;

using namespace skymaps;


int main(int , char** )
{

    int rc(0);
    try{

#if 0 // test code for Exposure--need a cube for built-in test

        std::string testexp("d:\\users\\kerrm\\Comparison\\expCube.fits");
        Exposure exp(testexp);
        double t = exp(SkyDir());
        std::cout << "exposure check: " << t << std::endl;
#endif

#if 1  // test code for the DiffuseFunction
        std::string path( ::getenv("EXTFILESSYS"));
//        DiffuseFunction df( path + "/galdiffuse/GP_gamma_v0r0p1.fits");
        DiffuseFunction df( path + "/galdiffuse/GP_gamma.fits");
        SkyDir gc(0,0, SkyDir::GALACTIC);
        double e1(1000), e2(2000);
        double t1 ( df(gc, e1) );  // differential at 1000
        double t2 ( df(gc, e2) ); // check integral
        double power ( log(t1/t2)/log(e2/e1) );
        double t12 = df(gc, e1, e2); // integral?
        double check ( t1 * e1 *(pow(e2/e1, 1-power)-1)/(1-power) );
        double diff( t12/check-1 );
        if( fabs(diff)> 1e-3 ) {
            std::cout << "Fail diffuse function check:" << check << std::endl;
            rc =1 ;
        }
        
#endif
        std::cout << "\ntesting BinnedPhotonData to create all standard Band objects" << std::endl;
    PhotonBinner binner;
    BinnedPhotonData* bpd= new BinnedPhotonData(binner);
    bpd->addPhoton(Photon(SkyDir(0,0),50, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),150, 0, 0));
    bpd->addPhoton(Photon(SkyDir(3,0),150, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),300, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),600, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),1000, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),2000, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),4000, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),10000, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),20000, 0, 0));
    bpd->addPhoton(Photon(SkyDir(0,0),100000, 0, 0));
    // a few back guys
    bpd->addPhoton(Photon(SkyDir(0,0),150, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),300, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),1000, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),100000, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),2000, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),3000, 0, 1));
    bpd->addPhoton(Photon(SkyDir(0,0),10000, 0, 1));
    bpd->info();
    // now check that we can find somthing
    typedef  std::vector<std::pair< SkyDir,  int> > PixelVec;
    PixelVec vec;

    BinnedPhotonData::const_iterator it=bpd->begin();
    ++it; ++it;  // to 3rd band
    const Band& b (*it );
    int n = b.query_disk( SkyDir(),  5.*M_PI/180 ,vec);
    std::cout << "found " << vec.size() << " pixels in band "<<int(b) << std::endl;
    for( PixelVec::const_iterator it = vec.begin(); it !=vec.end(); ++it){
        SkyDir r (it->first);
        std::cout << "\t( "<< r.ra() << ", " << r.dec() << " ): " << it->second << std::endl;
    }

    std::cout << "density : " << bpd->density(SkyDir()) << std::endl;
    std::cout << "value at 110 MeV: " << bpd->value(SkyDir(), 110) << std::endl;

    std::cout << "Writing a FITS file: " << std::endl;

    bpd->write("binned.fits");

    std::cout << "reading it back" << std::endl;
    BinnedPhotonData back(std::string("binned.fits"));
    back.info();

    }catch(const std::exception& e){
        std::cerr << "Caught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        rc=1;
    }
     return rc;
}

