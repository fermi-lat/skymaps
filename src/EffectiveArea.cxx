/** @file EffectiveArea.cxx
@brief implement EffectiveArea

$Header$
*/

#include "skymaps/EffectiveArea.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace{
// this code is copied from the latResponse package, mostly authored by J. Chiang 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    size_t binIndex(double x, const std::vector<float> & xx) {
        std::vector<float>::const_iterator ix = 
            std::upper_bound(xx.begin(), xx.end(), x);
        return ix - xx.begin();
    }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    class Array {
    public:
        Array(const std::vector<float> & values, size_t nx) 
            : m_values(values), m_nx(nx) {}
            float operator()(size_t iy, size_t ix) const {
                return m_values.at(iy*m_nx + ix);
            }
    private:
        const std::vector<float> & m_values;
        size_t m_nx;
    };
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**  
 * @class Bilinear
 *
 * @brief Bilinear interpolator.  Hold references to the vectors so
 * that introspection and greater transparency is available during
 * debugging.
 *
 */

class Bilinear {

public:

   Bilinear(const std::vector<float> & x, const std::vector<float> & y, 
            const std::vector<float> & values);

   Bilinear(const std::vector<float> & x, const std::vector<float> & y, 
            const std::vector<float> & values, 
            float xlo, float xhi, float ylo, float yhi);

   double operator()(float x, float y) const;

private:

   std::vector<float> m_x;
   std::vector<float> m_y;
   std::vector<float> m_values;
   
};


    Bilinear::Bilinear(const std::vector<float> & x, 
                   const std::vector<float> & y,
                   const std::vector<float> & values) 
   : m_x(x), m_y(y), m_values(values) {}

Bilinear::Bilinear(const std::vector<float> & x, 
                   const std::vector<float> & y,
                   const std::vector<float> & values,
                   float xlo, float xhi, float ylo, float yhi) {
   m_x.resize(x.size() + 2);
   std::copy(x.begin(), x.end(), m_x.begin() + 1);
   m_x.front() = xlo;
   m_x.back() = xhi;

   m_y.resize(y.size() + 2);
   std::copy(y.begin(), y.end(), m_y.begin() + 1);
   m_y.front() = ylo;
   m_y.back() = yhi;

   Array array(values, x.size());
   m_values.push_back(array(0, 0));
   for (size_t i(0); i < x.size(); i++) {
      m_values.push_back(array(0, i));
   }
   m_values.push_back(array(0, x.size()-1));
   for (size_t j(0); j < y.size(); j++) {
      m_values.push_back(array(j, 0));
      for (size_t i(0); i < x.size(); i++) {
         m_values.push_back(array(j, i));
      }
      m_values.push_back(array(j, x.size()-1));
   }
   m_values.push_back(array(y.size()-1, 0));
   for (size_t i(0); i < x.size(); i++) {
      m_values.push_back(array(y.size()-1, i));
   }
   m_values.push_back(array(y.size()-1, x.size()-1));
}

double Bilinear::operator()(float x, float y) const {
   typedef std::vector<float>::const_iterator const_iterator_t;

   const_iterator_t ix(std::upper_bound(m_x.begin(), m_x.end(), x));
   if (ix == m_x.end() && x != m_x.back()) {
      throw std::invalid_argument("Bilinear::operator: x out of range");
   }
   if (x == m_x.back()) {
      ix = m_x.end() - 1;
   } else if (x <= m_x.front()) {
      ix = m_x.begin() + 1;
   }
   int i(ix - m_x.begin());
    
   const_iterator_t iy(std::upper_bound(m_y.begin(), m_y.end(), y));
   if (iy == m_y.end() && y != m_y.back()) {
      throw std::invalid_argument("Bilinear::operator: y out of range");
   }
   if (y == m_y.back()) {
      iy = m_y.end() - 1;
   } else if (y <= m_y.front()) {
      iy = m_y.begin() + 1;
   }
   int j(iy - m_y.begin());

   double tt((x - m_x.at(i-1))/(m_x.at(i) - m_x.at(i-1)));
   double uu((y - m_y.at(j-1))/(m_y.at(j) - m_y.at(j-1)));

   size_t xsize(m_x.size());

   double y1(m_values.at(xsize*(j-1) + (i-1)));
   double y2(m_values.at(xsize*(j-1) + (i)));
   double y3(m_values.at(xsize*(j) + (i)));
   double y4(m_values.at(xsize*(j) + (i-1)));

   double value = ( (1. - tt)*(1. - uu)*y1 
                    + tt*(1. - uu)*y2  
                    + tt*uu*y3 
                    + (1. - tt)*uu*y4 ); 
   return value;
}
}//anon namespace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using namespace skymaps;
   
class EffectiveArea::FitsTable {

public:

   FitsTable(const std::string & filename,
             const std::string & extname,
             const std::string & tablename,
             size_t nrow=0);

   FitsTable();

   FitsTable(const FitsTable & other);

   ~FitsTable();
      
   /// @brief lookup a value from the table
   /// @param logenergy log10(energy)
   /// @param costh cos(theta)
   /// @param interpolate [true] if true, make linear
   /// interpolation. Otherwise, return value for given cell.
   double value(double logenergy, double costh, bool interpolate=true) const;
    
   double maximum() const {
      return m_maxValue;
   }
   
   double minCosTheta() const {
      return m_minCosTheta;
   }

   static void getVectorData(const tip::Table * table,
                             const std::string & fieldName,
                             std::vector<float> & values,
                             size_t nrow=0);

protected:

   /// Disable copy assignment operator.
   FitsTable & operator=(const FitsTable &) {
      return *this;
   }

private:

   Bilinear * m_interpolator;

   std::vector<float> m_logEnergies; 
   std::vector<float> m_mus; 
   std::vector<float> m_values;

   std::vector<float> m_ebounds;
   std::vector<float> m_tbounds;
   
   float m_minCosTheta;

   float m_maxValue;

};

EffectiveArea::FitsTable::FitsTable(const std::string & filename,
                     const std::string & extname,
                     const std::string & tablename,
                     size_t nrow) : m_interpolator(0) {

   const tip::Table * table(tip::IFileSvc::instance().readTable(filename, 
                                                                extname));

   std::vector<float> elo, ehi;
   getVectorData(table, "ENERG_LO", elo, nrow);
   getVectorData(table, "ENERG_HI", ehi, nrow);
   for (size_t k(0); k < elo.size(); k++) {
      m_ebounds.push_back(std::log10(elo.at(k)));
      m_logEnergies.push_back(std::log10(std::sqrt(elo.at(k)*ehi.at(k))));
   }
   m_ebounds.push_back(std::log10(ehi.back()));

   std::vector<float> mulo, muhi;
   getVectorData(table, "CTHETA_LO", mulo, nrow);
   getVectorData(table, "CTHETA_HI", muhi, nrow);
   for (size_t i(0); i < muhi.size(); i++) {
      m_tbounds.push_back(mulo.at(i));
      m_mus.push_back((m_tbounds.at(i) + muhi.at(i))/2.);
   }
   m_tbounds.push_back(muhi.back());

   m_minCosTheta = mulo.front();

   getVectorData(table, tablename, m_values, nrow);
   m_maxValue = m_values.front();
   for (size_t i(1); i < m_values.size(); i++) {
      if (m_values.at(i) > m_maxValue) {
         m_maxValue = m_values.at(i);
      }
   }

// Replicate nasty THF2 and RootEval::Table behavior from handoff_response,
// by passing xlo, xhi, ylo, yhi values
   float xlo, xhi, ylo, yhi;
   m_interpolator = new Bilinear(m_logEnergies, m_mus, m_values,
                                 xlo=0., xhi=10., ylo=-1., yhi=1.);

   delete table;
}

EffectiveArea::FitsTable::FitsTable() : m_interpolator(0) {}

EffectiveArea::FitsTable::FitsTable(const FitsTable & rhs) 
   : m_interpolator(0), m_logEnergies(rhs.m_logEnergies), m_mus(rhs.m_mus),
     m_values(rhs.m_values), m_ebounds(rhs.m_ebounds),
     m_tbounds(rhs.m_tbounds), m_minCosTheta(rhs.m_minCosTheta), 
     m_maxValue(rhs.m_maxValue) {
   m_interpolator = new Bilinear(m_logEnergies, m_mus, m_values,
                                 0, 10, -1, 1);
}

EffectiveArea::FitsTable::~FitsTable() { 
   delete m_interpolator;
}

double EffectiveArea::FitsTable::
value(double logenergy, double costh, bool interpolate) const {
   if (interpolate) {
      if (costh > m_mus.back()) {
         costh = m_mus.back();
      }
      return (*m_interpolator)(logenergy, costh);
   }

   double maxloge(*(m_logEnergies.end() - 2)); 
   if (logenergy >= maxloge) { // use last bin
      logenergy = maxloge;
   }
   if (logenergy <= m_logEnergies.at(1)) { // use first bin
      logenergy = m_logEnergies.at(1);
   }

   size_t ix = binIndex(logenergy, m_ebounds);
   size_t iy = binIndex(costh, m_tbounds);
   if (iy == 0) {
      iy = 1;
   } 
   size_t indx = (iy - 1)*m_logEnergies.size() + ix - 1;

   return m_values.at(indx);
}

void EffectiveArea::FitsTable::getVectorData(const tip::Table * table,
                              const std::string & fieldName,
                              std::vector<float> & values,
                              size_t nrow) {
   values.clear();

   tip::Table::ConstIterator it(table->begin());
   tip::ConstTableRecord & row(*it);

   for (size_t i(0); i < nrow; i++) {
      ++it;
   }
   row[fieldName].get(values);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




std::string EffectiveArea::s_CALDB;
void EffectiveArea::set_CALDB(std::string CALDB){ s_CALDB=CALDB;}

EffectiveArea::EffectiveArea(std::string irf_name)
: m_simple(false)
, m_aeffTable(0)
{
    if( irf_name=="simple"){
        m_simple= true;
        return;
    }
    if( s_CALDB.empty()){
        const char* c(::getenv("CALDB") );
        if( c==0){
            throw std::invalid_argument("EffectiveArea:: CALDB is not set");
        }
        s_CALDB = std::string(c);
    }
    std::string infile(s_CALDB+"/bcf/ea/aeff_"+irf_name+".fits");
    static std::string table_name("EFFECTIVE AREA");
    try{
        //const tip::Table * table = tip::IFileSvc::instance().readTable(infile, table_name, "");
        m_aeffTable = new FitsTable(infile, table_name, "EFFAREA" );
    }catch(const std::exception& e){
        std::cerr << "EffectiveArea: could not open " << infile<< "["<<table_name <<"]" << std::endl;
        throw;
    }

}

EffectiveArea::~EffectiveArea()
{
    delete m_aeffTable;
}

double EffectiveArea::value(double energy, double costheta)const
{
    if(m_simple) {
        static double ctmin(0.2);
        return costheta>ctmin? 4000.*(costheta-ctmin)/(1.-ctmin) : 0 ;
    }
    bool interpolate;
    return m_aeffTable->value(std::log10(energy), costheta, interpolate=true)*1e4;
}


