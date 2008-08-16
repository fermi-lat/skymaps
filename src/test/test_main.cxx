/** @file likelihood_main.cxx

$Header$

*/
//#include "skymaps/PhotonMap.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/Exposure.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/PhotonBinner.h"

#include "skymaps/ExposureCube.h"
#include "healpix/HealpixArrayIO.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"


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

#if 1
        {
        // try to sample a different image
        SkyImage in("G:\\kamae\\high_resolution_av.fits");
        in.reimage(SkyDir(180,0,SkyDir::GALACTIC), "D:\\common\\temp\\anticenter05a.fits", 0.05, 20, "CAR", true);
        }

#endif

#if 0
        // test code for ExposureCube -- make an exposure cube
        bool zenith(false);
        ExposureCube ex;  // default parameters
        std::string infile("F:\\glast\\data\\first_light\\dubois\\FT2_236511638-239811535_merged.fits");
        std::string outfile("D:\\common\\first_light\\earth\\sky_cube.fits");
        std::string outtable("Exposure");
        tip::Table * scData = tip::IFileSvc::instance().editTable(infile, "SC_DATA");
        if(zenith) ex.useZenith();
        std::cout << "\nCreating "<< (zenith?"zenith":"sky") <<"frame exposure cube:\nat:\t" << outfile <<"\nfrom:\t" << infile << std::endl;
        ex.load(scData);
        healpix::HealpixArrayIO::instance().write(ex.data(), outfile, outtable);
#endif

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

        int bins_per_decade(0);
        std::cout << "\ntesting BinnedPhotonData to create Band objects, bins/decade = "<<bins_per_decade << std::endl;
        PhotonBinner binner(bins_per_decade); // bins per decade
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

        std::cout << "Testing query_disk:\n  band        ra     dec  count \n" ;
        std::cout << std::fixed << std::right;
        using std::setw;
        using std::setprecision;
        BinnedPhotonData::const_iterator it=bpd->begin();
        for(int n(0) ; it!=bpd->end(); ++it,++n){
            PixelVec vec;
            const Band& b (*it );
            int m = b.query_disk( SkyDir(),  5.*M_PI/180 ,vec);
            std::vector<int>neighbors;
            b.findNeighbors(0, neighbors);
            std::cout << setw(4) << n ;
            for( PixelVec::const_iterator it = vec.begin(); it !=vec.end(); ++it){
                SkyDir r (it->first);
                std::cout << "\t"<< setw(8)<<setprecision(3)<< r.ra() << setw(8)<< r.dec() << setw(4) << it->second << std::endl;
            }
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

