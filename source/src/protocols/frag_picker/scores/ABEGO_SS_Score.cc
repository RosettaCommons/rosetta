// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite && is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ABEGO_SS_Score.cc
/// @brief  scores a fragment by secondary structure similarity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/ABEGO_SS_Score.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Map.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <core/fragment/SecondaryStructure.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// project headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

static thread_local basic::Tracer trABEGO_SS_Score(
		"protocols.frag_picker.scores.ABEGO_SS_Score");

bool ABEGO_SS_Score::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	Real totalScore = 0.0;
	for (Size i = 1; i <= f->get_length(); i++) {
		VallChunkOP chunk = f->get_chunk();
		VallResidueOP r = chunk->at(f->get_first_index_in_vall() + i - 1);
		char s(r->ss());
		Size bin = quota::torsion2big_bin_id(r->phi(),r->psi(),r->omega());
		for(Size j=1;j<=maps_.size();j++) {
		    if(maps_[j]->check_status(s,bin)) {
			totalScore += 1.0 - ratios_[i][j];
			break;
		    }
		}
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}


void ABEGO_SS_Score::do_caching(VallChunkOP chunk) {

	std::string tmp = chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;

	trABEGO_SS_Score.Debug << "caching ABEGO-SS score for " << chunk->get_pdb_id()
			<< " of size " << chunk->size() << std::endl;

	utility::vector1<Size> chunk_bins_;
	for (Size j = 1; j <= chunk->size(); ++j) {
	    VallResidueOP r = chunk->at(j);
	    char s(r->ss());
	    Size bin = quota::torsion2big_bin_id(r->phi(),r->psi(),r->omega());
	    bool if_found = false;
	    for(Size i=1;i<=maps_.size();i++) {
		if(maps_[i]->check_status(s,bin)) {
		    chunk_bins_.push_back(i);
		    if_found = true;
		    break;
		}
		continue;
	    }
	    if(!if_found) {
		trABEGO_SS_Score.Warning << "Can't find a feature label for the following combination of ss,abego: "
		    <<s<<","<<maps_[1]->abego_char(bin)<<std::endl;
		assert(false);
	    }
	}
	for (Size i = 1; i <= query_len_; ++i) {
		for (Size j = 1; j <= chunk->size(); ++j) {
			scores_[i][j] = 1.0 - ratios_[i][ chunk_bins_[j] ];
		}
	}
	trABEGO_SS_Score.Debug << "precomputed matrix of scores " << scores_.size()
			<< "x" << chunk->size() << std::endl;
}

bool ABEGO_SS_Score::cached_score(FragmentCandidateOP f,
		FragmentScoreMapOP empty_map) {

	std::string tmp = f->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(f->get_chunk());

	Real totalScore = 0;
	for (Size i = 1; i <= f->get_length(); ++i) {
		assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		assert(f->get_first_index_in_vall()	+ i - 1<= scores_[1].size());
		totalScore += scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1];
	}

	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

FragmentScoringMethodOP MakeABEGO_SS_Score::make(Size priority, Real lowest_acceptable_value, bool use_lowest,
			FragmentPickerOP picker, std::string /*prediction_file*/) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Size vall_max_len = picker->get_vall()->get_largest_chunk_size();

	if ( !option[ in::file::torsion_bin_probs ].user() ) {
		utility_exit_with_message("Error: no input file specified for ABEGO_SS score; use in::file::torsion_bin_probs option");
	}
	return (FragmentScoringMethodOP) new ABEGO_SS_Score(priority,
		lowest_acceptable_value, use_lowest,option[ in::file::torsion_bin_probs ](),
		vall_max_len);
}

} // scores
} // frag_picker
} // protocols


