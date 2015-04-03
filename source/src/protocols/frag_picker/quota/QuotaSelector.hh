// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/QuotaSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_QuotaSelector_hh
#define INCLUDED_protocols_frag_picker_quota_QuotaSelector_hh

#include <protocols/frag_picker/quota/QuotaSelector.fwd.hh>
#include <protocols/frag_picker/quota/QuotaCollector.hh>
#include <protocols/frag_picker/quota/QuotaPool.hh>

#include <protocols/frag_picker/FragmentSelectingRule.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <core/types.hh>

// C++ headers
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace quota {

/// @brief selects a given number of fragments using a quota scheme
class QuotaSelector: public FragmentSelectingRule {
public:

	/// @brief  Constructor sets the desired number of fragments.
	QuotaSelector(Size,Size,QuotaCollectorOP);

	virtual ~QuotaSelector() {
	}

	/// @brief  Selects desired number of fragments from a given candidates
	virtual void select_fragments(
				ScoredCandidatesVector1 const& in,
				ScoredCandidatesVector1 & out)
	{
		select_fragments_25_200(in,out);
	}

protected:
	QuotaCollectorOP collector_;
	Size q_pos_;
	inline Size round(Real x) { return Size(x > 0.0 ? x + 0.5 : x - 0.5); }

	Size next_from_pool(
		 ScoredCandidatesVector1 const&,
		 Size recently_taken,
		 std::set<Size> & in_use
	);

	void push_the_limits(utility::vector1<Size> & q_limits,Size target_total);

	void push_the_limits_to_the_winner(utility::vector1<Size> & q_limits,Size target_total);

	void push_the_limits(utility::vector1<Size> &,Size,utility::vector1<Real> &);

	virtual void select_fragments_200(
				ScoredCandidatesVector1 const& in,
				ScoredCandidatesVector1 & out);

	virtual void select_fragments_25_200(
				ScoredCandidatesVector1 const& in,
				ScoredCandidatesVector1 & out);
private:
    utility::vector1<std::string> tags_;
    std::map<std::string,Size> tag_map_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_QuotaSelector_HH */
