// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CompositeFragmentSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CompositeFragmentSelector_hh
#define INCLUDED_protocols_frag_picker_CompositeFragmentSelector_hh

//#include <protocols/frag_picker/CompositeFragmentSelector.fwd.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/vector1.hh>

//#include <vector>


namespace protocols {
namespace frag_picker {

/// @brief selects fragments by running several selectors
class CompositeFragmentSelector: public FragmentSelectingRule {
public:

	/// @brief  Constructor sets the desired number of fragments.
	CompositeFragmentSelector(core::Size frags_per_pos) : FragmentSelectingRule(frags_per_pos) {
	}

	/// @brief  Selects desired number of fragments from a given set of candidates
	void select_fragments( ScoredCandidatesVector1 const&, ScoredCandidatesVector1& ) override;

	~CompositeFragmentSelector() override = default;

	void add_selector(FragmentSelectingRuleOP new_selector) { selectors_.push_back( new_selector ); }
private:
	utility::vector1<FragmentSelectingRuleOP> selectors_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_CompositeFragmentSelector_HH */
