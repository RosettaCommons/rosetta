// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/PartialSecondarySimilarity.cc
/// @brief  scores a fragment by secondary structure similarity
///         but throws out the worst 20% of the residue matches
/// @author Robert Vernon (rvernon@u.washington.edu)

#include <protocols/frag_picker/scores/PartialSecondarySimilarity.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <core/fragment/SecondaryStructure.hh>


// project headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

static THREAD_LOCAL basic::Tracer trPartialSecondarySimilarity(
	"protocols.frag_picker.scores.PartialSecondarySimilarity");

void PartialSecondarySimilarity::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;

	//assert(query_ss_);
	utility::vector1<Size> chunk_ss_id( chunk->size() );
	for ( Size j = 1; j <= chunk->size(); ++j ) {
		char s(chunk->at(j)->ss());
		if ( s == 'H' ) chunk_ss_id[j] = 1;
		if ( s == 'E' ) chunk_ss_id[j] = 2;
		if ( s == 'L' ) chunk_ss_id[j] = 3;
	}

	for ( Size i = 1; i <= query_len_; ++i ) {
		for ( Size j = 1; j <= chunk->size(); ++j ) {
			scores_[i][j] = raw_probs_[i][chunk_ss_id[j]];
		}
	}
	trPartialSecondarySimilarity.Debug << "precomputed matrix of scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}

bool PartialSecondarySimilarity::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	utility::vector1< Real > values;
	values.resize(f->get_length());

	for ( Size i = 1; i <= f->get_length(); i++ ) {
		values[i] = 0.0;
		VallChunkOP chunk = f->get_chunk();

		char s(chunk->at(f->get_first_index_in_vall() + i - 1)->ss());
		if ( s == 'H' ) {
			values[i] =  raw_probs_[f->get_first_index_in_query() + i - 1][1];
		} else {
			if ( s == 'E' ) {
				values[i] =  raw_probs_[f->get_first_index_in_query() + i - 1][2];
			} else {
				if ( s == 'L' ) {
					values[i] =  raw_probs_[f->get_first_index_in_query() + i - 1][3];
				}
			}
		}
	}

	std::sort( values.begin(), values.end() );
	Real totalScore = 0.0;

	for ( Size i = 1; i <= f->get_length(); i++ ) {
		//~1 at 1, ~0.05 at f->get_length, 0.5 at 0.7*f->get_length()
		Real sigmoid_weight( 1 / ( 1 + exp( (10*( (Real) i ) / f->get_length()) - 7 ) ) );

		totalScore += sigmoid_weight*values[i];
	}

	totalScore /= (Real) f->get_length();

	// Real best80 = ( static_cast< Real > ( f->get_length() ) ) * 0.8;
	//  Size best80i =  static_cast< Size > ( best80 );

	//  Real totalScore = 0.0;
	//  for (Size i = 1; i<= best80i; i++) {
	//   totalScore += values[i];
	//  }
	//  totalScore /= (Real) best80i;

	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}


bool PartialSecondarySimilarity::cached_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	//REMOVED CACHING FOR NOW!
	score(f, empty_map);

	return true;
	// std::string & tmp = f->get_chunk()->chunk_key();
	//  if (tmp.compare(cached_scores_id_) != 0) {
	//   do_caching(f->get_chunk());
	//  }

	//  Real totalScore = cache_[f->get_length()][f->get_first_index_in_query()][f->get_first_index_in_vall()];

	//  totalScore /= (Real) f->get_length();

	//  empty_map->set_score_component(totalScore, id_);
	//  if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
	//   return false;
	//  return true;
}

PartialSecondarySimilarity::PartialSecondarySimilarity(Size priority, Real lowest_acceptable_value, bool use_lowest,
	core::fragment::SecondaryStructureOP query_prediction, std::string prediction_name,
	Size sequence_length, Size longest_vall_chunk) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
	"PartialSecondarySimilarity") , prediction_name_(prediction_name) {
	query_len_ = sequence_length;
	//query_ss_ = query_prediction;

	//Normalize Secondary Structure Scores
	// so: 0.5 H 0.5 E 0.0 L ---> 1.0 H 1.0 E 0.0 L
	norm_query_H_.resize(query_len_);
	norm_query_E_.resize(query_len_);
	norm_query_L_.resize(query_len_);
	for ( Size i = 1; i <= query_len_; ++i ) {
		Real highest_ss_pred(query_prediction->helix_fraction(i));

		if ( highest_ss_pred < query_prediction->strand_fraction(i) ) {
			highest_ss_pred = query_prediction->strand_fraction(i);
		}
		if ( highest_ss_pred < query_prediction->loop_fraction(i) ) {
			highest_ss_pred = query_prediction->loop_fraction(i);
		}

		norm_query_H_[i] = query_prediction->helix_fraction(i) / highest_ss_pred;
		norm_query_E_[i] = query_prediction->strand_fraction(i) / highest_ss_pred;
		norm_query_L_[i] = query_prediction->loop_fraction(i) / highest_ss_pred;
	}


	//scores_ vector = query_size x vall_size for per residue scores
	//prow vector = query_size x 3, used to store 1->0 scores for the 0->1 probabilities
	for ( Size i = 1; i <= query_len_; ++i ) {
		utility::vector1<Real> row(longest_vall_chunk);
		scores_.push_back(row);
		utility::vector1<Real> prow(3);
		prow[1] = 1 - norm_query_H_[i];
		prow[2] = 1 - norm_query_E_[i];
		prow[3] = 1 - norm_query_L_[i];
		raw_probs_.push_back(prow);
	}

}

} // scores
} // frag_picker
} // protocols


