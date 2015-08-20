// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/FormFactor.hh
/// @brief Collects histogram of distances between atoms of given SAXS-types, e.g. distances between  ARG-CEN and TRP-CEN
/// @details The histogram is used in fast SAXS spectrum evaluation (slightly approximate)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_saxs_DistanceHistogram_hh
#define INCLUDED_core_scoring_saxs_DistanceHistogram_hh

#include <core/scoring/saxs/DistanceHistogram.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>


#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {

/// @brief
class DistanceHistogram: public utility::pointer::ReferenceCount {
public:

	/// @brief  Create a histogram that collects distances between two predefined atom types
	DistanceHistogram(void) :
		factor_(10.0),
		max_dist_(300),
		h_(factor_*max_dist_, 0) // apl -- I really think you want to initialize the counts to 0
	{}

	/// @brief  Returns the number of counts for a given distance
	Size operator()(Real distance) const { return h_[(Size)(distance*factor_)]; }

	/// @brief  Returns the number of counts for a given distance
	Size get(Real distance) const  { return h_[(Size)(distance*factor_)]; }

	/// @brief  Returns the number of counts for a given bin
	Size get(Size bin) const  { return h_[bin]; }

	/// @brief tells waht distance falls into a certain bin
	Real distance(Size bin) const { return bin/factor_; }

	/// @brief Adds a distance observation to the histogram
	void inline insert(Real distance) {
	  Size b = (Size)(distance*factor_);
	  if( b>last_nonempty_bin_ ) last_nonempty_bin_ = b;
	  if ( b>0) {
	    ++h_[b];
		}
	}

	/// @brief Returns the size of this histogram
	Size size() { return h_.size(); }

	/// @brief Clears this histogram by filling each cell with 0.0
	void zeros() { for(Size i=1;i<=h_.size();i++) { h_[i] = 0; } last_nonempty_bin_ = 0; }

	/// @brief Returns the total number of counts in this histogram
	Size total() const {
	  Size cnt = 0;
	  for(Size i=1;i<=h_.size();i++)
	    cnt+=h_[i];
	  return cnt;
	}

	Size last_nonempty_bin() const { return last_nonempty_bin_; }

private:
	Real factor_;
	Real max_dist_;
	Size last_nonempty_bin_;
	utility::vector1<Size> h_;
};

} // core
} // scoring
} // saxs

#endif /* INCLUDED_core_scoring_saxs_DistanceHistogram_HH */
