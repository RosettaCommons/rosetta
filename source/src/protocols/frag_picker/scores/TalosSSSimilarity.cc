// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/TalosSSSimilarity.cc
/// @brief  scores a fragment by secondary structure similarity
/// @author rvernon@u.washington.edu

#include <protocols/frag_picker/scores/TalosSSSimilarity.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/fragment_scoring_utilities.hh>

#include <core/fragment/SecondaryStructure.hh>


// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// project headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trTalosSSSimilarity(
		"protocols.frag_picker.scores.TalosSSSimilarity");

bool TalosSSSimilarity::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	utility::vector1< Real > values;
	values.resize(f->get_length());

	Real totalScore = 0.0;
	for (Size i = 1; i <= f->get_length(); i++) {
		//mjo commenting out 'ss_weight' because it is unused and causes a warning
		//Real ss_weight(0.0);
		VallChunkOP chunk = f->get_chunk();

		char s(chunk->at(f->get_first_index_in_vall() + i - 1)->ss());
		Size ss_id(3);
		if (s == 'H') ss_id = 1;
		if (s == 'E') ss_id = 2;
		if (s == 'L') ss_id = 3;

		values[i] = raw_probs_[f->get_first_index_in_query()+i-1][ss_id];

		//totalScore += ss_weight;
	}

	//H
 	std::sort( values.begin(), values.end() );
 	totalScore = 0.0;
 	for (Size i = 1; i <= f->get_length(); i++) {
 		//~1 at 1, ~0.05 at f->get_length, 0.5 at 0.7*f->get_length()
 		Real sigmoid_weight( 1 / ( 1 + exp( (10*( (Real) i ) / f->get_length()) - 7 ) ) );
 		totalScore += sigmoid_weight*values[i];
 	}//H


	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}


void TalosSSSimilarity::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;

	do_caching_simple(chunk);
	for(Size fl=1;fl<=cache_.size();fl++) {
	    if(cache_[fl].size() != 0) {
		trTalosSSSimilarity.Trace << "caching secondary score for " << chunk->get_pdb_id()
			<< " of size " << chunk->size() << " for fragment size "<<fl<<std::endl;
		rolling_score(scores_,fl,cache_[fl]);
	    }
	}
}

void TalosSSSimilarity::do_caching_simple(VallChunkOP chunk) {

	assert(query_ss_);
	utility::vector1<Size> chunk_ss_id( chunk->size() );
	for (Size j = 1; j <= chunk->size(); ++j) {
	    char s(chunk->at(j)->ss());
	    if (s == 'H') chunk_ss_id[j] = 1;
	    if (s == 'E') chunk_ss_id[j] = 2;
	    if (s == 'L') chunk_ss_id[j] = 3;
	}

	for (Size i = 1; i <= query_len_; ++i) {
		for (Size j = 1; j <= chunk->size(); ++j) {
			scores_[i][j] = raw_probs_[i][chunk_ss_id[j]];
		}
	}
	trTalosSSSimilarity.Debug << "precomputed matrix of scores " << scores_.size()
			<< "x" << chunk->size() << std::endl;
}


bool TalosSSSimilarity::cached_score(FragmentCandidateOP f,
		FragmentScoreMapOP empty_map) {

	//return score(f, empty_map);

	/*
	std::string & tmp = f->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(f->get_chunk());
	*/


	Real totalScore = cache_[f->get_length()][f->get_first_index_in_query()][f->get_first_index_in_vall()];

	totalScore /= (Real) f->get_length();

	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
	   return false;
	return true;
}

TalosSSSimilarity::TalosSSSimilarity(Size priority, Real lowest_acceptable_value, bool use_lowest,
			core::fragment::SecondaryStructureOP query_prediction, std::string prediction_name,
			Size sequence_length, utility::vector1<Size> & frag_sizes, Size longest_vall_chunk) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
				"TalosSSSimilarity") , prediction_name_(prediction_name) {
		query_len_ = sequence_length;
		query_ss_ = query_prediction;

		H_mult_ = option[frags::seqsim_H](); // Default is 1.0
		E_mult_ = option[frags::seqsim_E](); // Default is 1.0
		L_mult_ = option[frags::seqsim_L](); // Default is 1.0

		for (Size i = 1; i <= query_len_; ++i) {
			utility::vector1<Real> row(longest_vall_chunk);
			scores_.push_back(row);
			utility::vector1<Real> prow(3);
			//prow[1] = 1 - H_mult_*query_prediction->helix_fraction(i);
			//prow[2] = 1 - E_mult_*query_prediction->strand_fraction(i);
			//prow[3] = 1 - L_mult_*query_prediction->loop_fraction(i);

			Real hf;
			Real sf;
			Real lf;
			Real confidence;
			if ( i<=query_prediction->total_residue() ) {
				hf = query_prediction->helix_fraction(i);
				sf = query_prediction->strand_fraction(i);
				lf = query_prediction->loop_fraction(i);
				confidence = query_prediction->confidence(i);
			} else {
				hf = sf = lf = 1.0/3.0;
				confidence = 0.0;
			}

			Real average( (hf+sf+lf)/ 3.0 );
			Real sdev( sqrt( (pow(hf-average,2) + pow(sf-average,2) + pow(lf-average,2))/3.0 ) );

			// Real highest_ss_pred(query_prediction->helix_fraction(i));

// 			if (highest_ss_pred < query_prediction->strand_fraction(i)) {
// 				highest_ss_pred = query_prediction->strand_fraction(i);
// 			}
// 			if (highest_ss_pred < query_prediction->loop_fraction(i)) {
// 				highest_ss_pred = query_prediction->loop_fraction(i);
// 			}

// 			if ( highest_ss_pred > 0 ) {
// 				hf = hf / highest_ss_pred;
// 				sf = sf / highest_ss_pred;
// 			}

			if ( ( sdev <= 0.1 ) || (confidence <= 0.1) ) {
				prow[1] = 0.0;//-3*( 1 / ( 1 + exp(-5*hf + 5) ) );
				prow[2] = 0.0;//-3*( 1 / ( 1 + exp(-5*sf + 5) ) );
				prow[3] = 0.0;// - query_prediction->loop_fraction(i);
			} else {
				//prow[1] = -3*( 1 / ( 1 + exp(-5*hf + 5) ) )*confidence;
				//prow[2] = -3*( 1 / ( 1 + exp(-5*sf + 5) ) )*confidence;
				//prow[3] = 0.0;// - query_prediction->loop_fraction(i);

				prow[1] = (-2*( 1 / ( 1 + exp(-7*hf + 5) ) ))*sqrt(confidence);
				prow[2] = (-2*( 1 / ( 1 + exp(-7*sf + 5) ) ))*sqrt(confidence);
				prow[3] = (-2*( 1 / ( 1 + exp(-7*lf + 5) ) ))*sqrt(confidence);


				//prow[1] =
				//prow[2] =
				//prow[3] = 0.0;// - query_prediction->loop_fraction(i);
			}

			//sep_2fx20.6_px_m2
			//+2 when prediction is 0.0 (SS Penalty)
			//-2 when prediction is 1.0
			// 0 when prediction is 0.33
			//prow[1] = -(2 * ( (1 / ( 1 + exp((-20*hf + 6)))) + hf) - 2);
			//prow[2] = -(2 * ( (1 / ( 1 + exp((-20*sf + 6)))) + sf) - 2);
			//prow[3] = 0.0;// - query_prediction->loop_fraction(i);


			//BEST SO FAR!
			//prow[1] = -3*( 1 / ( 1 + exp(-5*hf + 5) ) );
			//	prow[2] = -3*( 1 / ( 1 + exp(-5*sf + 5) ) );
			//prow[3] = 0.0;// - query_prediction->loop_fraction(i);


			//sep_2fx20.6_px_m2
			//+2 when prediction is 0.0 (SS Penalty)
			//-2 when prediction is 1.0
			// 0 when prediction is 0.33
			//prow[1] = -1*( (1 / ( 1 + exp((-25*hf + 5)))) + (1 / ( 1 + exp((-12.5*hf + 9)))) - 1);
			//prow[2] = -1*( (1 / ( 1 + exp((-25*sf + 5)))) + (1 / ( 1 + exp((-12.5*sf + 9)))) - 1);
			//prow[3] = 0.0;// - query_prediction->loop_fraction(i);


			//if (prow[1] < 0.0) {
			//	prow[1] = 0.0;
			//}
			//if (prow[2] < 0.0) {
			//	prow[2] = 0.0;
			//}
			raw_probs_.push_back(prow);
		}

		create_cache(frag_sizes,query_len_,longest_vall_chunk,cache_);
}

} // scores
} // frag_picker
} // protocols


