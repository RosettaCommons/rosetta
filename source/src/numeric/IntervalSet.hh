// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/IntervalSet.hh
/// @brief  definition and implementation of IntervalSet class
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_numeric_IntervalSet_hh
#define INCLUDED_numeric_IntervalSet_hh

// Numeric headers
#include <numeric/IntervalSet.fwd.hh>
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/vector0.hh>

// C++ Headers
#include <ostream>

namespace numeric {

template <typename T>
class IntervalSet {

private:

	utility::vector0<T> endpoints_;

public:

	inline
	IntervalSet()
	{}

	inline
	IntervalSet(
		T start,
		T end
	):
		endpoints_(2)
	{
		assert(start <= end);

		endpoints_[0] = start;
		endpoints_[1] = end;
	}

	inline
	IntervalSet(
		T start1,
		T end1,
		T start2,
		T end2
	):
		endpoints_(4)
	{
		assert(start1 <= end1 && end1 <= start2 && start2 <= end2);

		endpoints_[0] = start1;
		endpoints_[1] = end1;
		endpoints_[2] = start2;
		endpoints_[3] = end2;
	}

	/// @brief vector of interval set end points
	utility::vector0<T> const &
	endpoints() const
	{
		return endpoints_;
	}

	/// @brief clear the contents
	inline
	void
	clear()
	{
		endpoints_.clear();
	}

	/// @brief set the inverval set to contain a single interval
	inline
	void
	set(
		T start,
		T end
	)
	{
		endpoints_.resize(2);

		assert(start <= end);

		endpoints_[0] = start;
		endpoints_[1] = end;
	}

	/// @brief set the interval set to contain two intervals
	inline
	void
	set(
		T start1,
		T end1,
		T start2,
		T end2
	)
	{
		endpoints_.resize(4);

		assert(start1 <= end1 && end1 <= start2 && start2 <= end2);

		endpoints_[0] = start1;
		endpoints_[1] = end1;
		endpoints_[2] = start2;
		endpoints_[3] = end2;
	}

	/// @brief add an interval to the end of the set
	inline
	void
	push_back(
		T start,
		T end
	)
	{
		int const size = endpoints_.size();
		endpoints_.resize(size+2);

		assert(start <= end && (size == 0 || start >= endpoints_[size-1]));

		endpoints_[size] = start;
		endpoints_[size+1] = end;
	}

	/// @brief calculate the total length of all the intervals
	inline
	T
	length()
	{
		T len = 0;
		int const size = endpoints_.size();
		for ( int i = 0; i < size; i += 2 ) {
			len += endpoints_[i+1] - endpoints_[i];
		}

		return len;
	}

	/// @brief determine if a point is within one of the intervals
	inline
	bool
	is_inside(
		T point
	)
	{
		int const size = endpoints_.size();
		for ( int i = 0; i < size; i += 2 ) {
			if ( point >= endpoints_[i] && point <= endpoints_[i+1] ) {
				return true;
			}
		}

		return false;
	}

	/// @brief calculate the intersection of two IntervalSets
	/// @param[in] right - the second Interval set
	/// @return a new IntervalSet
	IntervalSet
	operator&(
		IntervalSet const & right
	)
	{
		IntervalSet<T> newintervals;
		utility::vector0<T> const & endpoints_right = right.endpoints();
		int const size1 = endpoints_.size();
		int const size2 = endpoints_right.size();
		int i1 = 0;
		int i2 = 0;

		while ( i1 < size1 && i2 < size2 ) {

			// Determine which interval starts first
			if ( endpoints_[i1] < endpoints_right[i2] ) {

				if ( endpoints_[i1+1] < endpoints_right[i2] ) {
					// If there's no overlap, advance to the next endpoints_ interval
					i1 += 2;
				} else if ( endpoints_[i1+1] < endpoints_right[i2+1] ) {
					// If there's some overlap, push the intersection
					newintervals.push_back(endpoints_right[i2], endpoints_[i1+1]);
					// Advance to the next endpoints_ interval
					i1 += 2;
				} else {
					// If there's total overlap, push the whole endpoints_right interval
					newintervals.push_back(endpoints_right[i2], endpoints_right[i2+1]);
					// Advance to the next endpoints_right interval
					i2 += 2;
				}
			} else {

				if ( endpoints_right[i2+1] < endpoints_[i1] ) {
					// If there's no overlap, advance to the next endpoints_right interval
					i2 += 2;
				} else if ( endpoints_right[i2+1] < endpoints_[i1+1] ) {
					// If there's some overlap, push the intersection
					newintervals.push_back(endpoints_[i1], endpoints_right[i2+1]);
					// Advance to the next endpoints_right interval
					i2 += 2;
				} else {
					// If there's total overlap, push the whole endpoints_ interval
					newintervals.push_back(endpoints_[i1], endpoints_[i1+1]);
					// Advance to the next endpoints_ interval
					i1 += 2;
				}
			}
		}

		return newintervals;
	}

	/// @brief pick a random number uniformly from all the intervals
	/// @param[in] RG - random number generator to use
	/// @return random number
	T
	random_point(random::RandomGenerator & RG)
	{
		T const len = length();
		int const size = endpoints_.size();
		int i = 0;

		T rand_0_1(RG.uniform());

		assert(rand_0_1 >= 0 && rand_0_1 <= 1);

		T randnum = endpoints_[i] + rand_0_1*len;

		while ( i+3 < size && randnum > endpoints_[i+1] ) {
			randnum += endpoints_[i+2] - endpoints_[i+1];
			i += 2;
		}

		if ( randnum == endpoints_[i+1] && i+3 < size ) {
			int j = 1;
			while ( i+2*j+3 < size && endpoints_[i+2*j+2] == endpoints_[i+2*j+3] ) j++;
			int randpoint = RG.random_range(0, j);
			if ( randpoint > 0 ) randnum = endpoints_[i+2*randpoint];
		}

		return randnum;
	}
};

template <typename T>
inline
std::ostream &
operator<<(
	std::ostream & output,
	const IntervalSet<T> & interval
)
{
	utility::vector0<T> const & endpoints = interval.endpoints();

	for ( int i = 0; i < (signed)endpoints.size(); i += 2 ) {
		if ( i != 0 ) output << " ";
		output << "[" << endpoints[i] << ", " << endpoints[i+1] << "]";
	}

	return output;
}

// PyRosetta WorkAround
class IntervalSet_Double : public IntervalSet<double>
{
};


} // namespace numeric

#endif // INCLUDED_numeric_IntervalSet_HH
