// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/BestTotalScoreSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_BestTotalScoreSelector_hh
#define INCLUDED_protocols_frag_picker_BestTotalScoreSelector_hh

#include <protocols/frag_picker/BestTotalScoreSelector.fwd.hh>

#include <protocols/frag_picker/FragmentSelectingRule.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CommonFragmentComparators.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

// utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

/// @brief selects a given number of fragments using a quota scheme
class BestTotalScoreSelector: public FragmentSelectingRule {
public:

	/// @brief  Constructor sets the desired number of fragments.
	BestTotalScoreSelector(Size frags_per_pos,
		scores::FragmentScoreManagerOP scoring_scheme) :
		FragmentSelectingRule(frags_per_pos), comparator_(scoring_scheme) {

	}

	virtual ~BestTotalScoreSelector() {
	}

	/// @brief  Selects desired number of fragments from a given candidates
	virtual void select_fragments( ScoredCandidatesVector1 const&, ScoredCandidatesVector1& );

private:
	CompareTotalScore comparator_;
};


} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_BestTotalScoreSelector_HH */
