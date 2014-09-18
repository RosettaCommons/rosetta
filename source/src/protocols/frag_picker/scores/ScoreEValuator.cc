// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite && is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ScoreEValuator.cc
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ScoreEValuator.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// option key includes
// AUTO-REMOVED #include <core/init/init.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>
#include <basic/Tracer.hh>

#include <algorithm>

#include <utility/vector1.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

namespace protocols {
namespace frag_picker {
namespace scores {


using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trProfScore(
		"protocols.frag_picker.scores.ScoreEValuator");

void ScoreEValuator::do_caching(VallChunkOP chunk) {

	ProfileScore::do_caching(chunk);
}

void ScoreEValuator::clean_up() {

	ProfileScore::clean_up();
}

bool ScoreEValuator::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	Real totalScore = 0;
	Real mean = 0;
	Real stdev = 0;

	utility::vector1<Size> columnsQ;
	utility::vector1<Size> columnsV;
	for (Size i = 1; i <= f->get_length(); i++) {
		assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		assert(f->get_first_index_in_vall()
				+ i - 1<= scores_[1].size());
		totalScore
				+= scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall()
						+ i - 1];
		columnsQ.push_back(f->get_first_index_in_query() + i - 1);
		columnsV.push_back(f->get_first_index_in_vall() + i - 1);
	}

	for (Size i_rand = 1; i_rand <= max_rand_; ++i_rand) {
		numeric::random::random_permutation(columnsQ.begin(), columnsQ.end(), numeric::random::rg());
		numeric::random::random_permutation(columnsV.begin(), columnsV.end(), numeric::random::rg());
		Real s = 0;
		for (Size i = 1; i <= columnsQ.size(); i++) {
			s = scores_[columnsQ[i]][columnsV[i]];
			mean += s;
			stdev += s * s;
		}
	}
	mean /= ((Real) max_rand_);
	stdev /= ((Real) max_rand_);
	stdev = sqrt(stdev - mean * mean);

	totalScore = (totalScore - mean) / stdev;

	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

FragmentScoringMethodOP MakeScoreEValuator::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker) {

	if (option[frags::scoring::profile_score].user()) {
		core::sequence::ScoringSchemeFactory ssf;
		core::sequence::ScoringSchemeOP ss(ssf.get_scoring_scheme(
				option[frags::scoring::profile_score]()));
		Size len = picker->get_vall()->get_largest_chunk_size();
		trProfScore << "Profile scoring method is: "
				<< option[frags::scoring::profile_score]() << std::endl;
		return (FragmentScoringMethodOP) new ScoreEValuator(priority,
				lowest_acceptable_value, use_lowest, picker->get_query_seq(), ss, len);
	}
	utility_exit_with_message(
			"[ERROR] Undefined profile scoring method. Provide it with frags::scoring_scheme flag");

	return NULL;
}

} //scores
} // frag_picker
} // protocols
