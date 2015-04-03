// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ScoreEValuator.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ScoreEValuator_hh
#define INCLUDED_protocols_frag_picker_scores_ScoreEValuator_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/ProfileScore.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  Computes E-Value for a Profile object
class ScoreEValuator: public ProfileScore {

public:
	ScoreEValuator(Size priority, Real lowest_acceptable_value,
			bool use_lowest, sequence::SequenceProfileOP query_profile,
			sequence::ScoringSchemeOP profile_scoring, Size largest_chunk_size) :
		ProfileScore(priority, lowest_acceptable_value, use_lowest, query_profile,
				profile_scoring, largest_chunk_size) {
		max_rand_ = 12;
	}

	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	//bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map, std::ostream& out);
private:
	Size max_rand_;
};

class MakeScoreEValuator: public MakeFragmentScoringMethod {
public:

	MakeScoreEValuator() :
		MakeFragmentScoringMethod("ScoreEValuator") {
	}

	using MakeFragmentScoringMethod::make;

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP);
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_ScoreEValuator_HH */
