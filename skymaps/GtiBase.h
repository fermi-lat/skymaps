/** \file GtiBase.h
    \brief Encapsulation of the concept of a GTI. May be constructed from a GTI extension.

    \author James Peachey, HEASARC/GSSC
    Copied from evtbin::Gti by T.Burnett
    $Header$
*/
#ifndef skymaps_GtiBase_h
#define skymaps_GtiBase_h

#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace skymaps {

  /** \class GtiBase
      \brief Encapsulation of the concept of a GTI. May be constructed from a GTI extension.
  */
  class GtiBase {
    public:
      typedef std::pair<double, double> Interval_t;
      typedef std::map<double, double> IntervalCont_t;
      typedef IntervalCont_t::iterator Iterator;
      typedef IntervalCont_t::const_iterator ConstIterator;

      /// \brief Construct a GTI object which contains no intervals.
      GtiBase();

      /** \brief Construct a GTI object which contains intervals specified by vectors.
          \param starts A series of start times.
          \param stops A series of stop times.
      */
      GtiBase(std::vector<double> & starts, std::vector<double> & stops);

      /** \brief Construct a GTI object, reading its intervals from the given file.
          \param file_name The input file.
          \param ext_name The extension in the file containing the GTI information.
      */
      GtiBase(const std::string & file_name, const std::string & ext_name = std::string("GTI"));

      /** \brief Compute fraction of time interval bounded by tstart and tstop which falls
          within one or more of the Good time intervals.
          \param tstart Interval start time.
          \param tstop Interval start time.
          \param gti_pos Iterator pointing to the first good time interval which might contain
          part of the given interval (a hint on which GTI to use for increased efficiency).
      */
      double getFraction(double tstart, double tstop, ConstIterator & gti_pos) const;

      /** \brief From this GtiBase and a second, produce a third GtiBase which is the intersection of
          the two inputs.
          \param gti The other gti.
      */
      GtiBase operator &(const GtiBase & gti) const;

      /** \brief From this GtiBase and a second, produce a third GtiBase which is the union of
          the two inputs.
          \param gti The other gti.
      */
      GtiBase operator |(const GtiBase & gti) const;

      /** \brief Replace this GtiBase with the intersection of this GtiBase and the one supplied as an argument.
          \param gti The other gti.
      */
      GtiBase & operator &=(const GtiBase & gti);

      /** \brief Replace this GtiBase with the union of this GtiBase and the one supplied as an argument.
          \param gti The other gti.
      */
      GtiBase & operator |=(const GtiBase & gti);

      /** \brief Compare two sets of time intervals. If they are identical, this returns true.
          \param gti The GtiBase object with which this one will be compared.
      */
      bool operator ==(const GtiBase & gti) const;

      /** \brief Compare two sets of time intervals. If they differ in any way, this returns true.
          \param gti The GtiBase object with which this one will be compared.
      */
      bool operator !=(const GtiBase & gti) const;

      /// Return iterator pointing to first time interval in GTI.
      Iterator begin();

      /// Return iterator pointing to one past last time interval in GTI.
      Iterator end();

      /// Return const iterator pointing to first time interval in GTI.
      ConstIterator begin() const;

      /// Return const iterator pointing to one past last time interval in GTI.
      ConstIterator end() const;

      /** \brief Add an interval to this GTI object.
          \param tstart The start time of the interval.
          \param tstop The start time of the interval.
      */
      void insertInterval(double tstart, double tstop);

      /// \brief Return the number of intervals in this GTI container.
      int getNumIntervals() const;

      /** \brief Set the number of intervals in this GTI container.
          \param num_intv The new size of the GTI container.
      */
      void setNumIntervals(int num_intv);

      /// \brief Compute and return the ontime (total of all intervals).
      double computeOntime() const;

      /// \brief Write itervals to a stream.
      void write(std::ostream & os) const;

    protected:
      /// \brief Merge overlapping intervals.
      void consolidate();

      /** \brief For internal use only. Add an interval to this GTI object. Does not call consolidate to ensure 
            intervals are merged to form minimal set.
          \param tstart The start time of the interval.
          \param tstop The start time of the interval.
      */
      void insertInterval(Interval_t interval);

    private:
      IntervalCont_t m_intervals;
  };

  std::ostream & operator <<(std::ostream & os, const GtiBase & gti);
}

#endif
