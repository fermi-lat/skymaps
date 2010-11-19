#include "skymaps/PythonUtilities.h"

#include "CLHEP/Vector/ThreeVector.h"
#include <stdexcept>
#include <assert.h>
//#include "timeSystem/GlastTimeHandler.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace skymaps {

void PythonUtilities::val_grid (double* rvals, int rvals_size,
                                const std::vector<double>& lons, const std::vector<double>&lats, 
                                const astro::SkyDir& roi_center, const astro::SkyFunction& myfunc)
{
    astro::SkyDir rot_axis(astro::SkyDir(roi_center.l()-90,0));
    Hep3Vector rot_axisv(rot_axis.dir());
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
    Hep3Vector rot_axisv(rot_axis.dir());
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

}
