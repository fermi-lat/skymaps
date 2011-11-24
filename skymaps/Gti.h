/**
 * @file Gti.h
 * @brief Handle GTIs. Derives from GtiBase, adding a constructor and 
 * an accept() method.
 * @author J. Chiang
 *
 * Copied by T. Burnett
 * $Header$
 */

#ifndef skymaps_Gti_h
#define skymaps_Gti_h

#include "skymaps/GtiBase.h"

namespace tip {
   class Table;
}

namespace skymaps {
/**
 * @class Gti
 * @brief A more useful and complete implementation of GtiBase
 * @author J. Chiang
 *
 */

    class Gti : public skymaps::GtiBase {

public:

   Gti() : skymaps::GtiBase() {}

   Gti(const std::string & filename, const std::string & extension="GTI") 
      : skymaps::GtiBase(filename, extension) {}

   Gti(const tip::Table & gtiTable);

   Gti(const skymaps::GtiBase & gti) : skymaps::GtiBase(gti) {}

   Gti(std::vector<double> starts,std::vector<double> stops) : skymaps::GtiBase(starts,stops) {}

   bool accept(double time) const;

   void writeExtension(const std::string & filename) const;

   /// @brief Given a start and stop time this method recomputes the
   /// the GTIs.
   /// @return A new Gti object with the new good-time intervals.
   /// @param start The start time of the time range cut (MET seconds)
   /// @param stop The stop time of the time range cut (MET seconds)
   Gti applyTimeRangeCut(double start, double stop) const;

   /// @return The minimum lower bound of the GTIs (MET seconds)
   double minValue() const;

   /// @return The maximum upper bound of the GTIs (MET seconds)
   double maxValue() const;

};
} // namespace skymaps

#endif // astro_Gti_h
