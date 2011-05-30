#include "skymaps/PythonUtilities.h"
#include "skymaps/WeightedSkyDir.h"

#include "CLHEP/Vector/ThreeVector.h"
#include <stdexcept>
#include <assert.h>
#include <cstdlib>
#include <string.h>
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/BaryTimeComputer.h"
extern "C" {
double *glastscorbit (char *, double, int *) ;       /* GLAST-specific */
}

namespace skymaps {

void PythonUtilities::val_grid (double* rvals, int rvals_size,
                                const std::vector<double>& lons, const std::vector<double>&lats, 
                                const astro::SkyDir& roi_center, const astro::SkyFunction& myfunc)
{
    astro::SkyDir rot_axis(astro::SkyDir(roi_center.l()-90,0));
    CLHEP::Hep3Vector rot_axisv(rot_axis.dir());
    double rot_extent(roi_center.b()*(M_PI/180));
    for (std::vector<double>::const_iterator it_lon = lons.begin(); it_lon != lons.end(); ++it_lon){
        for (std::vector<double>::const_iterator it_lat = lats.begin(); it_lat != lats.end(); ++it_lat){
            astro::SkyDir sd(astro::SkyDir(*it_lon,*it_lat));
            sd().rotate(rot_extent,rot_axisv);
            *rvals++ = myfunc(astro::SkyDir(sd.ra(),sd.dec(),astro::SkyDir::GALACTIC));
        }
    }
}


void PythonUtilities::rot_grid(double *rvals, int rvals_size,
                          const std::vector<astro::SkyDir> sdirs,
                          const astro::SkyDir& roi_center)
{
    assert(rvals_size == sdirs.size()*2); 
    astro::SkyDir rot_axis(astro::SkyDir(roi_center.l()+90,0));
    CLHEP::Hep3Vector rot_axisv(rot_axis.dir());
    double rot_extent = roi_center.b()*(M_PI/180);
    for (std::vector<astro::SkyDir>::const_iterator it(sdirs.begin()); it != sdirs.end(); ++it) {
        astro::SkyDir sd(astro::SkyDir(it->l(),it->b()));
        sd().rotate(rot_extent,rot_axisv);
        *rvals++ = sd.ra();
        *rvals++ = sd.dec();
    }
}

void PythonUtilities::arclength(double *rvals, int rvals_size,
                          const std::vector<astro::SkyDir> sdirs,
                          const astro::SkyDir& roi_center)
{
    if (rvals_size != sdirs.size()) throw 20;
    for (std::vector<astro::SkyDir>::const_iterator it(sdirs.begin()); it != sdirs.end(); ++it) {
        *rvals++ = roi_center.difference(*it);
    }
}

void PythonUtilities::get_float_col(double *rvals, int rvals_size,const std::string file_name,
                                    const std::string table_name, const std::string col_name)
{
    const tip::Table * table = tip::IFileSvc::instance().readTable(file_name, table_name, "");
    tip::Table::ConstIterator itend(table->end());
    for(tip::Table::ConstIterator it(table->begin()); it != itend; ++ it) {
        (*it)[col_name].get(*rvals);
        rvals++;
    }
    delete table;
}


void PythonUtilities::tdb2met(double *rvals, int rvals_size,
                              double ra, double dec,
                              const std::string sc_file,
                              double tol) {

        const timeSystem::BaryTimeComputer & tbc = timeSystem::BaryTimeComputer::getComputer("JPL DE405");
        char *ft2_name;
        ft2_name = new char[sc_file.size()+1];
        strcpy(ft2_name,sc_file.c_str());
        int error(0);
        double met_offset(51910+64.184/86400);
        double guess_frac(600./86400);
        for (int i=0; i<rvals_size; i++) {
            double mjd_tdb(*(rvals+i));
            timeSystem::AbsoluteTime original("TDB",int(mjd_tdb),(mjd_tdb-int(mjd_tdb))*86400);
            // do a bisection search to the required tolerance
            double hi(mjd_tdb+guess_frac),lo(mjd_tdb-guess_frac);
            double* sc_position;
            bool success(false);
            double pos_tol(double(1000*tol)/86400);
            for (int j=0; j<100; j++) {
                double new_time((hi+lo)/2);
                timeSystem::AbsoluteTime guess("TT",int(new_time),(new_time-int(new_time))*86400);
                if( (j==0) || ((hi-lo)>pos_tol)) {
                    sc_position = glastscorbit(ft2_name,(new_time-met_offset)*86400,&error);
                //if (error==0) throw 20; // don't know error code
                }
                std::vector<double> spos(sc_position,sc_position+3);
                tbc.computeBaryTime(ra,dec,spos,guess);
                double ddiff(original.computeElapsedTime("TDB",guess).getDuration("Sec"));
                if (abs(ddiff) < tol) {
                    rvals[i] = (new_time-met_offset)*86400;
                    success = true;
                    break;
                }
                if (ddiff < 0) hi = new_time;
                else if (ddiff > 0) lo = new_time;
                else throw 20;
            }
            if (!success) throw 20;
        }
    delete[] ft2_name;
}

void PythonUtilities::met2tdb(double *rvals, int rvals_size,
                        double ra, double dec,
                        const std::string sc_file) {
    const timeSystem::BaryTimeComputer & tbc = timeSystem::BaryTimeComputer::getComputer("JPL DE405");
    char *ft2_name;
    ft2_name = new char[sc_file.size()+1];
    strcpy(ft2_name,sc_file.c_str());
    int error(0);
    double met_offset(51910+64.184/86400);
    timeSystem::AbsoluteTime glast_tdb_origin("TDB",51910,64.184);
    for (int i=0; i<rvals_size; i++) {
        double met_mjd(rvals[i]/86400.+met_offset);
        timeSystem::AbsoluteTime at_met("TT",int(met_mjd),(met_mjd-int(met_mjd))*86400);
        double* sc_position = glastscorbit(ft2_name,rvals[i],&error);
        std::vector<double> spos(sc_position,sc_position+3);
        tbc.computeBaryTime(ra,dec,spos,at_met);
        rvals[i] = at_met.computeElapsedTime("TDB",glast_tdb_origin).getDuration("Sec");
    }
    delete[] ft2_name;
}

void PythonUtilities::met2geo(double *rvals, int rvals_size,
                        double ra, double dec,
                        const std::string sc_file) {
    const timeSystem::BaryTimeComputer & tbc = timeSystem::BaryTimeComputer::getComputer("JPL DE405");
    char *ft2_name;
    ft2_name = new char[sc_file.size()+1];
    strcpy(ft2_name,sc_file.c_str());
    int error(0);
    double met_offset(51910+64.184/86400);
    timeSystem::AbsoluteTime glast_geo_origin("TT",51910,64.184);
    for (int i=0; i<rvals_size; i++) {
        double met_mjd(rvals[i]/86400.+met_offset);
        timeSystem::AbsoluteTime at_met("TT",int(met_mjd),(met_mjd-int(met_mjd))*86400);
        double* sc_position = glastscorbit(ft2_name,rvals[i],&error);
        std::vector<double> spos(sc_position,sc_position+3);
        tbc.computeGeoTime(ra,dec,spos,at_met);
        rvals[i] = at_met.computeElapsedTime("TT",glast_geo_origin).getDuration("Sec");
    }
    delete[] ft2_name;
}

void PythonUtilities::get_wsdl_weights(double *rvals, int rvals_size,
                                       const BaseWeightedSkyDirList& wsdl) {
  assert(rvals_size=wsdl.size());
  for (BaseWeightedSkyDirList::const_iterator it = wsdl.begin(); it != wsdl.end(); ++it) {
        *rvals++ = it->weight();
  }
}

void PythonUtilities::set_wsdl_weights(const std::vector<double>& weights,
                                       BaseWeightedSkyDirList& wsdl) {
  assert(weights.size() == wsdl.size()); 
  std::vector<double>::const_iterator wt = weights.begin();
  for (BaseWeightedSkyDirList::iterator wsd = wsdl.begin(); wsd != wsdl.end(); ++wsd) {
    wsd->set_weight(*wt++);
  }
}

}
