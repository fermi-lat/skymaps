/** \file GtiBase.cxx
    \brief Implementation of encapsulation of the concept of a GTI. May be constructed from a GTI extension.
    \author James Peachey, HEASARC/GSSC
*/
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "skymaps/GtiBase.h"

#include "st_facilities/FileSys.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

//namespace skymaps {
using namespace skymaps;
GtiBase::GtiBase(): m_intervals() {}

GtiBase::GtiBase(std::vector<double> & starts, std::vector<double> & stops): m_intervals() {
    assert(starts.size()==stops.size());
    std::vector<double>::const_iterator it1(starts.begin());
    std::vector<double>::const_iterator it2(stops.begin());
    for (; it1 != starts.end(); ++it1, ++it2) {
        m_intervals.insert(Interval_t(*it1,*it2));
    }
    consolidate();
}

  GtiBase::GtiBase(const std::string & file_name, const std::string & ext_name): m_intervals() {
    using namespace st_facilities;

    // Get container of file names from the supplied input file.
    FileSys::FileNameCont file_cont = FileSys::expandFileList(file_name);

    // Iterate over input files.
    for (FileSys::FileNameCont::iterator itor = file_cont.begin(); itor != file_cont.end(); ++itor) {
      // Open GTI extension.
      std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(*itor, ext_name));

      // Fill container with intervals from the extension.
      for (tip::Table::ConstIterator itor = gti_table->begin(); itor != gti_table->end(); ++itor) {
        double start = (*itor)["START"].get();
        double stop = (*itor)["STOP"].get();
        if (start > stop) {
          std::ostringstream os;
          os << "GtiBase: In file " << file_name << ", record " << itor->getIndex() << " is invalid: " <<
            "start time " << start << " > stop time " << stop;
          throw std::runtime_error(os.str());
        }
        insertInterval(Interval_t(start, stop));
      }
    }
    consolidate();
  }

  double GtiBase::getFraction(double tstart, double tstop, ConstIterator & gti_pos) const {
    double fraction = 0.;

    for (; gti_pos != m_intervals.end(); ++gti_pos) {
      // Check if this interval ends before GTI starts and return 0. fraction if it does.
      if (tstop <= gti_pos->first) break;

      // Check if this interval is completely contained in the GTI and return 1 if it does.
      if (tstart >= gti_pos->first && tstop <= gti_pos->second) {
        if (tstop == gti_pos->second) ++gti_pos;
        fraction = 1.;
        break;
      }

      // Check if there is some overlap and add that overlap.
      if (tstart < gti_pos->second) {
        double start = tstart > gti_pos->first ? tstart : gti_pos->first;
        double stop = tstop < gti_pos->second ? tstop : gti_pos->second;
        fraction += (stop - start) / (tstop - tstart);
      }

      // Check if this GTI still has some part which might overlap some future interval.
      // If it does, break to avoid incrementing the GTI iterator.
      if (tstop < gti_pos->second) break;
    }

    return fraction;
  }

  GtiBase GtiBase::operator &(const GtiBase & old_gti) const {
    GtiBase new_gti;

    ConstIterator it1 = m_intervals.begin();
    ConstIterator it2 = old_gti.m_intervals.begin();

    // Iterate until either set of intervals is finished.
    while(it1 != m_intervals.end() && it2 != old_gti.m_intervals.end()) {
      // See if interval 1 comes before interval 2.
      if (it1->second <= it2->first) ++it1;
      // See if interval 2 comes before interval 1.
      else if (it2->second <= it1->first) ++it2;
      else {
        // They overlap, so find latest start time.
        double start = it1->first > it2->first ? it1->first : it2->first;

        // And earliest stop time.
        double stop = it1->second < it2->second ? it1->second: it2->second;

        new_gti.insertInterval(Interval_t(start, stop));

        // Skip to the next interval in the series for whichever interval ends earliest.
        if (it1->second < it2->second) ++it1; else ++it2;
      }
    }
    new_gti.consolidate();

    return new_gti;
  }

  GtiBase GtiBase::operator |(const GtiBase & gti) const {
    GtiBase new_gti = *this;
    new_gti |= gti;
    return new_gti;
  }

  GtiBase & GtiBase::operator &=(const GtiBase & gti) {
    *this = *this & gti;
    return *this;
  }

  GtiBase & GtiBase::operator |=(const GtiBase & gti) {
    for (ConstIterator itor = gti.begin(); itor != gti.end(); ++itor) {
      insertInterval(*itor);
    }
    consolidate();
    return *this;
  }

  bool GtiBase::operator ==(const GtiBase & gti) const { return m_intervals == gti.m_intervals; }

  bool GtiBase::operator !=(const GtiBase & gti) const { return m_intervals != gti.m_intervals; }

  GtiBase::Iterator GtiBase::begin() { return m_intervals.begin(); }

  GtiBase::Iterator GtiBase::end() { return m_intervals.end(); }

  GtiBase::ConstIterator GtiBase::begin() const { return m_intervals.begin(); }

  GtiBase::ConstIterator GtiBase::end() const { return m_intervals.end(); }

  void GtiBase::insertInterval(double tstart, double tstop) {
    m_intervals.insert(Interval_t(tstart, tstop));
    consolidate();
  }

  int GtiBase::getNumIntervals() const { return m_intervals.size(); }

  void GtiBase::setNumIntervals(int) {}

  double GtiBase::computeOntime() const {
    double on_time = 0.;
    for (IntervalCont_t::const_iterator itor = m_intervals.begin(); itor != m_intervals.end(); ++itor)
      on_time += itor->second - itor->first;
    return on_time;
  }

  void GtiBase::write(std::ostream & os) const {
    IntervalCont_t::const_iterator itor = m_intervals.begin();
    if (m_intervals.end() != itor) {
      os << "[" << itor->first << ", " << itor->second << "]";
      ++itor;
    }
    for (; itor != m_intervals.end(); ++itor) {
      os << std::endl;
      os << "[" << itor->first << ", " << itor->second << "]";
    }
  }

  void GtiBase::consolidate() {
    // Get iterator pointing to the first interval.
    Iterator current = m_intervals.begin();
    if (m_intervals.end() != current) {
      // Get iterator pointing to the next interval.
      Iterator next = current;

      // Consider each pair of intervals in the container.
      for (++next; next != m_intervals.end(); ++next) {
        if (current->first < next->first && next->first <= current->second) {
          // Next interval begins in the middle of the current interval.
          // If next interval continues past end of current interval, stretch
          // the current interval to cover the combined range.
          if (current->second < next->second) current->second = next->second;
          // Next interval is no longer needed in any case.
          m_intervals.erase(next);
          next = current;
        } else if (next->first < current->first && current->first <= next->first) {
          // Current interval begins in the middle of the next interval.
          // If current interval continues past end of next interval, stretch
          // the next interval to cover the combined range.
          if (next->second < current->second) next->second = current->second;
          // Current interval is no longer needed in any case.
          m_intervals.erase(current);
        }
        current = next;
      }
    }
  }

  void GtiBase::insertInterval(Interval_t interval) {
    // See if an interval already exists which has this start time.
    Iterator found = m_intervals.find(interval.first);
    if (m_intervals.end() == found) {
      m_intervals.insert(interval);
    } else {
      if (interval.second > found->second) found->second = interval.second;
    }
  }

  std::ostream & operator <<(std::ostream & os, const GtiBase & gti) {
    gti.write(os);
    return os;
  }

//}
