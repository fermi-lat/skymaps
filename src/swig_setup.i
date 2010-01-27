%module(docstring="Interface to skymaps") skymaps
%{
#include <stdexcept>
#include <vector>
#include <utility>

#include "astro/EarthCoordinate.h"
#include "astro/Photon.h"
#include "astro/IGRField.h"
#include "astro/JulianDate.h"
#include "astro/SolarSystem.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"
#include "healpix/HealpixMap.h"

#include "skymaps/SkyImage.h"

#include "skymaps/BinnedPhotonData.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/CompositeSkyFunction.h"
#include "skymaps/PhotonMap.h"
#include "skymaps/Exposure.h"
#include "skymaps/ExposureWeighter.h"
#include "skymaps/Convolution.h"
#include "skymaps/PsfFunction.h"
#include "skymaps/LivetimeCube.h"
#include "skymaps/PySkyFunction.h"
#include "skymaps/Gti.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/Background.h"
#include "skymaps/IsotropicPowerLaw.h"
#include "skymaps/PhotonBinner.h"
#include "skymaps/IParams.h"
#include "skymaps/HealpixDiffuseFunc.h"
#include "skymaps/PsfSkyFunction.h"
#include "skymaps/IsotropicSpectrum.h"
#include "skymaps/BandBackground.h"
#include "skymaps/SpectralFunction.h"
#include "skymaps/SkyIntegrator.h"
#include "skymaps/ComplexSkySpectrum.h"
#include "skymaps/SmoothedSkySpectrum.h"
#include "skymaps/PsfSkySpectrum.h"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/EulerAngles.h"

#include "CLHEP/Vector/ThreeVector.h"

%}
%include stl.i
%exception {
   try {
      $action
   } catch (const std::exception & eObj) {
      if( strcmp(eObj.what(),"StopIteration")==0 ){
          PyErr_SetString(PyExc_StopIteration, const_cast<char *>(eObj.what()));
      } else if(strcmp(eObj.what(),"IndexError")==0 ){
          PyErr_SetString(PyExc_IndexError, const_cast<char *>(eObj.what()));
      } else {
          PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(eObj.what()));
      }
      return NULL;
   }
}


%include std_vector.i
%include std_pair.i
%template(DoublePair) std::pair<double, double>;
%template(StringVector) std::vector<std::string>;
%template(DoubleVector) std::vector<double>;
%template(FloatVector) std::vector<float>;
%template(LongVector) std::vector<long>;
%template(IntVector) std::vector<int>;

%include CLHEP/Vector/ThreeVector.h
namespace CLHEP {
 class HepRotation {
public:
};
}
%extend CLHEP::Hep3Vector{
// for convenience: make it behave like array of 3 elements
   double __getitem__(size_t i) {
      switch (i){
      case 0: return self->x();
      case 1: return self->y();
      case 2: return self->z();
      case 3: throw std::range_error("StopIteration"); //must be exactly this string
      default: throw std::range_error("IndexError");
      }
   }
   size_t __len__() {      return 3;       }
}
%extend astro::SkyDir{
// for convenience: make it behave like array of 3 elements
   double __getitem__(size_t i) {
      switch (i){
      case 0: return self->dir().x();
      case 1: return self->dir().y();
      case 2: return self->dir().z();
      case 3: throw std::range_error("StopIteration"); //must be exactly this string
      default: throw std::range_error("IndexError");
      }
   }
   size_t __len__() { return 3; }
// make it pickle   
   %insert("python") %{
def __getstate__(self):  return (self.ra(),self.dec())
def __setstate__(self,tup):
    self.this = _skymaps.new_SkyDir(tup[0],tup[1])
    self.thisown=1
def __str__(self):  return 'SkyDir(%.3f,%.3f)' % (self.ra(),self.dec())
def __repr__(self): return self.__str__()

%}
}

%inline{
// Create a double * object, needed for SkyProj ctor
// modifed from http://www.swig.org/Doc1.1/HTML/Typemaps.html#n7
double * double2(double x, double y) {
	double *f = (double *) malloc(2*sizeof(double));
	f[0] = x;
	f[1] = y;
	return f;
}
// Destroy a double* object (is it called?)
void delete_double2(double *d) {
	free(d);
}
}

%extend CLHEP::HepRotation{
   double __getitem__(size_t i){
   switch(i){
      case 0: return self->xx(); case 1: return self->xy(); case 2: return self->xz();
      case 3: return self->yx(); case 4: return self->yy(); case 5: return self->yz();
      case 6: return self->zx(); case 7: return self->zy(); case 8: return self->zz();
      case 9: throw std::range_error("StopIteration"); //must be exactly this string
      default: throw std::range_error("IndexError");
      }
   }
   size_t __len__() {return 9;}
}

%extend healpix::HealpixMap{
   float __getitem__(size_t i)            {return (*self)[i];}
   double __call__(const astro::SkyDir& d){return (*self)(d);}
   size_t __len__() {return self->size();}
   void __setitem__(int i, double v){ (*self)[i] = v;}
}  
%include astro/SkyProj.h
%include astro/SkyDir.h
%include astro/Photon.h
%include astro/IGRField.h

%include astro/JulianDate.h
%include astro/SolarSystem.h

%include healpix/Healpix.h
%include healpix/HealpixMap.h

// fails now
//%include healpix/HealPixel.h

    
%include skymaps/BinnedPhotonData.h
%extend skymaps::BinnedPhotonData{
// provide access to the Band objects as an array
   skymaps::Band * __getitem__(size_t i){ 
      if( i == (*self).size() ) throw std::range_error("StopIteration");
      if( i<0 || i > self->size() ) throw std::range_error("IndexError");
      skymaps::BinnedPhotonData::iterator it= self->begin();
      for(int j(0); j!=i; ++j, ++it);
      return &(*it); // note address of
   }
   size_t __len__() {      return self->size();       }
}


%include skymaps/HealpixDiffuseFunc.h

%include skymaps/PhotonMap.h
%include skymaps/SkyImage.h

%include skymaps/Convolution.h
%include skymaps/PsfFunction.h
%include skymaps/PsfSkyFunction.h

%include skymaps/PySkyFunction.h
%include skymaps/Background.h
%include skymaps/Band.h
%include skymaps/SmoothedSkySpectrum.h
%include skymaps/PsfSkySpectrum.h

%extend skymaps::WeightedSkyDirList{
// provide access to the WeightedSKyDir objects as an array
   skymaps::WeightedSkyDir * __getitem__(size_t i){ 
      if( i == (*self).size() ) throw std::range_error("StopIteration");
      if( i<0 || i > self->size() ) throw std::range_error("IndexError");
      return &(*self)[i];
      //skymaps::WeightedSkyDirList::iterator it= self->begin();
      //for(int j(0); j!=i; ++j, ++it);
      //return &(*it); // note address of
   }
   size_t __len__() {      return self->size();       }
}

%include skymaps/PhotonBinner.h

%include skymaps/Gti.h
%extend skymaps::Gti{
   double computeOntime(){return (*self).GtiBase::computeOntime(); }
   int getNumIntervals() const {return (*self).GtiBase::getNumIntervals();}
   size_t __len__() {return self->getNumIntervals();}
   void insertInterval(double tstart, double tstop){ (*self).GtiBase::insertInterval(tstart, tstop); }
   void intersection(const skymaps::Gti& other){
     self->operator&=(other);
     }
     // note that the natural "union" is a C++ reserved word
   void combine(const skymaps::Gti& other){ self->operator|=(other);  }
   skymaps::Gti * __getitem__(int i){
      if( i == self->getNumIntervals() ) throw std::range_error("StopIteration");
      if(i<0) i+=self->getNumIntervals();
      if( i<0 || i > self->getNumIntervals() ) throw std::range_error("IndexError");
      skymaps::GtiBase::ConstIterator it(self->begin());
      for(int j(0); j!=i;++j,++it);
      skymaps::Gti * new_gti = new skymaps::Gti(self->applyTimeRangeCut(it->first,it->second));
      return new_gti;

      }

}

%include skymaps/CompositeSkySpectrum.h

%extend skymaps::CompositeSkySpectrum{

    double average(const astro::SkyDir& dir, double angle, double tolerance)
    {
        return self->average( dir, angle, tolerance);
    }
}

%extend skymaps::Background{

    void setEnergy(double e)
    {
        self->setEnergy(e);
    }
}

%extend skymaps::SmoothedSkySpectrum{
    void setEnergy(double e)
    {
        self->setEnergy(e);
    }
}

%include skymaps/CompositeSkyFunction.h
%extend skymaps::CompositeSkyFunction{

    double average(const astro::SkyDir& dir, double angle, double tolerance)
    {
        return self->average( dir, angle, tolerance);
    }
}

%include skymaps/SpectralFunction.h


%feature("kwargs");  // using keywords for these
%include skymaps/IsotropicSpectrum.h
%include skymaps/DiffuseFunction.h
%include skymaps/LivetimeCube.h
%include skymaps/Exposure.h
%include skymaps/EffectiveArea.h
%include skymaps/IsotropicPowerLaw.h
%include skymaps/IParams.h
%include skymaps/BandBackground.h
%include skymaps/SkyIntegrator.h
%include skymaps/ComplexSkySpectrum.h
%include skymaps/SmoothedSkySpectrum.h
%include skymaps/ExposureWeighter.h







