// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/PValuedFragmentScoreManager.hh
/// @brief  Score manager that stores statistics of all scores. This allows estimate p-value for a score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_PValuedFragmentScoreManager_hh
#define INCLUDED_protocols_frag_picker_scores_PValuedFragmentScoreManager_hh

// package headers
#include <protocols/frag_picker/scores/AdaptiveScoreHistogram.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief holds particular score components, weights and calculates the total score for a fragment candidate
/// @details a fragment picker object needs exactly one fragment manager to pick fragments. Adding new scoring methods
/// is supposed to be done FragmentPicker, which calls proper method from this class.
class PValuedFragmentScoreManager: public FragmentScoreManager {
public:

	/// @brief creates an empty manager
	PValuedFragmentScoreManager() {};

	/// @brief calculates all the small scores for a given fragment
	/// @brief results are properly stored inside a FragmentScoreMap object
	bool score_fragment_from_cache(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief calculates all the small scores for a given fragment
	/// @brief results are properly stored inside a FragmentScoreMap object
	bool score_fragment(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief prints a flat table with all scores for all the fragments in a given vector
	/// @details If the manager allows for annotations, they will be printed as well
	void describe_fragments(utility::vector1<std::pair<FragmentCandidateOP,
		scores::FragmentScoreMapOP> > const&, std::ostream&);

private:
	utility::vector1<utility::vector1<AdaptiveScoreHistogramOP> > statistics_;
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_PValuedFragmentScoreManager_HH */

