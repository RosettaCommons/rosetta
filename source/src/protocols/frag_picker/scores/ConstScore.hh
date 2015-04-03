// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ConstScore.hh
/// @brief  score that adds a constant (so basically it does nothing!)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ConstScore_hh
#define INCLUDED_protocols_frag_picker_scores_ConstScore_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  ConstScore adds a constant to the total score for each position
/// @details The total ConstScore for a fragment = n_frag_res * score_const
class ConstScore: public FragmentScoringMethod {
public:
	/// @brief the value used to score each position
	/// @details You don't have to change the value
	/// here as you can always rescale the score by a weight provided in a score configuration file
	static const int CONST_SCORE = 1;

	ConstScore(Size priority, Real lowest_acceptable_value, bool use_lowest) :
		FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest,"ConstScore") {
	}

	/// @brief Computes the score i.e. returns a constant times the number of residues in a fragment
	/// @details the method returns ALWAYS TRUE which is inconsistent with other methods derived from
	/// FragmentScoringMethod base class. The other option (i.e. returning false when a score is too high)
	/// doesn't make any sense in this case.
	virtual bool score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	    empty_map->set_score_component((Real) f->get_length() * CONST_SCORE, id_);
	    return true;
	}
};

/// @brief  Maker class that produces a new ConstScore object
class MakeConstScore: public MakeFragmentScoringMethod {
public:

	MakeConstScore() :
		MakeFragmentScoringMethod("ConstScore") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value,
			bool use_lowest, FragmentPickerOP, std::string /* params */) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ConstScore(priority,
				lowest_acceptable_value, use_lowest) );
	}
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ConstScore_HH */
