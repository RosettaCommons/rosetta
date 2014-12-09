// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.cc
/// @brief  scores a fragment by weighting L1 profile distances by residue type
/// @author Robert Vernon


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScoreDistWeight.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// project headers
#include <basic/Tracer.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>

// utils
#include <basic/prof.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trProfScoreDistWeight(
		"protocols.frag_picker.scores.ProfileScoreDistWeight");

void ProfileScoreDistWeight::do_caching(VallChunkOP chunk) {

	std::string tmp = chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;
	//Size size_q = query_profile_->length();

	trProfScoreDistWeight.Debug << "caching profile score for " << chunk->get_pdb_id()
			<< " of size " << chunk->size() << std::endl;
	PROF_START( basic::FRAGMENTPICKING_PROFILE_CAHING );
	//for (Size i = 1; i <= size_q; ++i) {
	//std::cout << "A " << query_sequence_ << " " << query_sequence_.length() << std::endl;
	for (Size i = 0; i < query_sequence_.length(); ++i) {
		//std::cout << "B" << chunk->size() << std::endl;
		//utility::vector1<Real> query_prof_row = query_profile_->prof_row(i);
		Size seqpos_res_id(aa_order_map_.find(query_sequence_[i])->second);


		for (Size j = 1; j <= chunk->size(); ++j) {
			//std::cout << "C" << chunk->size() << std::endl;

			//utility::vector1<Real> tmplt_prof_row = chunk->at(j)->profile();
			Real score(0.0);

			Size tmplt_res_id (aa_order_map_.find(chunk->at(j)->aa())->second);

			for (Size s = 1; s <= 3; s++) {
				Real ss_weight(0.0);
				if (s == 1)
					ss_weight = query_ss_->helix_fraction(i+1);
				if (s == 2)
					ss_weight = query_ss_->strand_fraction(i+1);
				if (s == 3)
					ss_weight = query_ss_->loop_fraction(i+1);

				//for (Size v = 1; v <= 20; v++) {
				//	Real distance_weight(0.0);
				//
				//	distance_weight = distance_weights_[s][seqpos_res_id][v];
				//	score += tmplt_prof_row[v]*ss_weight*distance_weight;
				//}

				Real distance_weight(0.0);

				distance_weight = distance_weights_[s][seqpos_res_id][tmplt_res_id];
				score += ss_weight*distance_weight;

					//score += std::abs(tmplt_prof_row[v] - query_prof_row[q])
					//	*ss_weight
					//	*distance_weight;
					//	}

				//std::cout << "SCOREIJ_FINAL " << score << std::endl;
			}

			//std::cout << "SCOREIJ_FINAL " << score << " " <<  i+1 << " " << j << " " << query_sequence_[i] << " " << chunk->at(j)->aa() << " " << query_ss_->helix_fraction(i+1) << " " << query_ss_->strand_fraction(i+1) << " " << query_ss_->loop_fraction(i+1) << " " << distance_weights_[1][seqpos_res_id][tmplt_res_id] << " " << distance_weights_[2][seqpos_res_id][tmplt_res_id] << " " << distance_weights_[3][seqpos_res_id][tmplt_res_id] << std::endl;

			scores_[i+1][j] = score;
		}
	}

	PROF_STOP( basic::FRAGMENTPICKING_PROFILE_CAHING );
	trProfScoreDistWeight.Debug << "precomputed matrix of scores " << scores_.size()
			<< "x" << chunk->size() << std::endl;
}

bool ProfileScoreDistWeight::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	std::string tmp = f->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(f->get_chunk());

	Real totalScore = 0.0;//f->get_length() * 20.0;
	for (Size i = 1; i <= f->get_length(); i++) {
		assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		assert(f->get_first_index_in_vall()
				+ i - 1<= scores_[1].size());
		totalScore += scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1];
		//std::cout << "TOTALSCORE " << totalScore << " " << scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1] << " " <<  f->get_first_index_in_query() + i - 1 << " " << f->get_first_index_in_vall() + i - 1 << std::endl;
	}
	totalScore *= (Real) 100.0;
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore < lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

FragmentScoringMethodOP MakeProfileScoreDistWeight::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string prediction_id) {

	//query_sequence_ = picker->get_query_seq_string();

	//std::cout << "QUERY_SEQUENCE " << query_sequence_ << std::endl;

	Size len = picker->get_vall()->get_largest_chunk_size();

	//std::istringstream line_stream(config_line);
	//std::string score_name;
	//Size p;
	//Real weight;
	//Real lowest;
	//std::string prediction_id;
	//line_stream >> score_name >> p >> weight >> lowest >> prediction_id;

	trProfScoreDistWeight << "Profile scoring method is: DistWeight" << std::endl;
	core::fragment::SecondaryStructureOP query_prediction( picker->get_query_ss(prediction_id) );
	if( ! query_prediction ) {
		utility_exit_with_message("Unable to find secondary structure prediction for " + prediction_id );
	}
	return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ProfileScoreDistWeight(priority,
																															lowest_acceptable_value, use_lowest, picker->get_query_seq(), query_prediction, picker->get_query_seq_string(),len) );
}

} //scores
} // frag_picker
} // protocols
