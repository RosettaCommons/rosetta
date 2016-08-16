// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaPool.hh
/// @brief provides a single pool used by quota
/// @details a QuotaSelector is constructed with a few pools. At fragment selection procedure
///  it tries to fit each fragment candidate into the pools by a round robin procedure
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_QuotaPool_hh
#define INCLUDED_protocols_frag_picker_quota_QuotaPool_hh

#include <protocols/frag_picker/quota/QuotaPool.fwd.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace frag_picker {
namespace quota {

/// @brief represents a single pool used by quota selector
class QuotaPool: public CandidatesCollector {
public:
	/// @brief Creates a pool of a given size and name
	/// @param name - name assigned to this pool. This in general may be any string that
	/// later allows one control pool's behavior from a flag file
	QuotaPool(std::string const & pool_name, Real quota_fraction):
		pool_name_(pool_name),
		quota_fraction_(quota_fraction)
	{}

	virtual ~QuotaPool() {};

	virtual bool could_be_accepted(ScoredCandidate) const = 0;

	/// @brief Says how many fragments (in total) may fit into this pool
	virtual Size total_size() const = 0;

	/// @brief Says how many fragments are currently in this pool
	virtual Size current_size() const = 0;

	/// @brief Says how many fragments can still be inserted into this pool
	virtual Size size_left() const = 0;

	/// @brief Makes the pool empty by removing all candidates
	virtual void clear() = 0;

	/// @brief Push a fragment candidate into the pool container
	virtual void push(ScoredCandidate) = 0;

	/// @brief  Check how many candidates have been already collected for a given position
	/// @details This is a very special case - collector will be used only for a given position.
	/// Thus it returns the total number of inserted candidates, as count_candidates() does
	virtual Size count_candidates() const = 0;

	/// @brief returns the name assigned to this quota pool
	inline
	std::string const &
	get_pool_name() const  {
		return pool_name_;
	}

	/// @brief prints information on which fragments can be accepted by this pool and how many of them
	/// @details base class' impementation says the capacity that has left
	virtual void show_availability(std::ostream & where) const {
		where << pool_name_<<" : "<< size_left()<<std::endl;
	}

	/// @brief returns the fraction of this quota pool in the entire population of fragments
	Real get_fraction() const { return quota_fraction_; }

	/// @brief Sets the fraction of this quota pool in the entire population of fragments
	virtual void set_fraction(Real new_fraction) { quota_fraction_ = new_fraction; }

	/// @brief provides the score for a candidate that was used to sort a quota pool
	/// @details This base class returns the most recent total score for a fragment
	inline virtual Real quota_score(ScoredCandidate candidate) const {
		return candidate.second->get_most_recent_total_score();
	}
private:
	std::string pool_name_;
	Real quota_fraction_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_QuotaPool_HH */
