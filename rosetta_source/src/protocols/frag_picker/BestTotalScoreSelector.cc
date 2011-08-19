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
// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <utility>

namespace protocols {
namespace frag_picker {

static basic::Tracer trBestTotalScoreSelector(
		"protocols.frag_picker.BestTotalScoreSelector");

void BestTotalScoreSelector::select_fragments(utility::vector1<std::pair<
		FragmentCandidateOP, scores::FragmentScoreMapOP> >& input_canditates,
		utility::vector1<std::pair<FragmentCandidateOP,
				scores::FragmentScoreMapOP> >& output_selection) {

	Size n = frags_per_pos();
	if (n > input_canditates.size())
		n = input_canditates.size();
	trBestTotalScoreSelector.Debug << "Selecting " << n << "fragments from "
			<< input_canditates.size() << " candidates" << std::endl;

	std::sort(input_canditates.begin(), input_canditates.end(), comparator_);
	for (Size i = 1; i <= n; i++) {
		output_selection.push_back(input_canditates[i]);
	}
}

} // frag_picker
} // protocols

