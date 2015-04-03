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

#include <protocols/frag_picker/BestTotalScoreSelector.hh>

#include <protocols/frag_picker/FragmentSelectingRule.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static thread_local basic::Tracer trBestTotalScoreSelector(
		"protocols.frag_picker.BestTotalScoreSelector");

void BestTotalScoreSelector::select_fragments(
   ScoredCandidatesVector1 const& input_candidates,
	 ScoredCandidatesVector1& output_selection )
{

	Size n = frags_per_pos();
	trBestTotalScoreSelector.Debug << "Selecting " << n << "fragments from "
			<< input_candidates.size() << " candidates" << std::endl;

	output_selection = input_candidates;
	if ( n > output_selection.size() ) {
		return;
	}
	std::sort( output_selection.begin(), output_selection.end(), comparator_ );
	output_selection.erase( output_selection.begin()+n,output_selection.end() );
}

} // frag_picker
} // protocols

