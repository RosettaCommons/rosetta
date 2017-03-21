// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/SecondarySimilarity.cc
/// @brief  scores a fragment by secondary structure similarity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/SecondarySimilarity.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/fragment_scoring_utilities.hh>

#include <core/fragment/SecondaryStructure.hh>


// project headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

static THREAD_LOCAL basic::Tracer trSecondarySimilarity(
	"protocols.frag_picker.scores.SecondarySimilarity");

bool SecondarySimilarity::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	core::Real totalScore = 0.0;
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		core::Real ss_weight(0.0);
		VallChunkOP chunk = f->get_chunk();

		char s(chunk->at(f->get_first_index_in_vall() + i - 1)->ss());
		if ( s == 'H' ) {
			ss_weight = 1 - query_ss_->helix_fraction(f->get_first_index_in_query() + i - 1);
		}
		if ( s == 'E' ) {
			ss_weight = 1 - query_ss_->strand_fraction(f->get_first_index_in_query() + i - 1);
		}
		if ( s == 'L' ) {
			ss_weight = 1 - query_ss_->loop_fraction(f->get_first_index_in_query() + i - 1);
		}

		totalScore += ss_weight;
	}
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}


void SecondarySimilarity::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;

	do_caching_simple(chunk);
	for ( core::Size fl=1; fl<=cache_.size(); fl++ ) {
		if ( cache_[fl].size() != 0 ) {
			trSecondarySimilarity.Trace << "caching secondary score for " << chunk->get_pdb_id()
				<< " of size " << chunk->size() << " for fragment size "<<fl<<std::endl;
			rolling_score(scores_,fl,cache_[fl]);
		}
	}
}

void SecondarySimilarity::do_caching_simple(VallChunkOP chunk) {

	debug_assert(query_ss_);
	utility::vector1<core::Size> chunk_ss_id( chunk->size() );
	for ( core::Size j = 1; j <= chunk->size(); ++j ) {
		char s(chunk->at(j)->ss());
		if ( s == 'H' ) chunk_ss_id[j] = 1;
		if ( s == 'E' ) chunk_ss_id[j] = 2;
		if ( s == 'L' ) chunk_ss_id[j] = 3;
	}

	for ( core::Size i = 1; i <= query_len_; ++i ) {
		for ( core::Size j = 1; j <= chunk->size(); ++j ) {
			scores_[i][j] = raw_probs_[i][chunk_ss_id[j]];
		}
	}
	trSecondarySimilarity.Debug << "precomputed matrix of scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}


bool SecondarySimilarity::cached_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	/*
	std::string & tmp = f->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
	do_caching(f->get_chunk());
	*/
	core::Real totalScore = cache_[f->get_length()][f->get_first_index_in_query()][f->get_first_index_in_vall()];

	totalScore /= (core::Real) f->get_length();

	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

SecondarySimilarity::SecondarySimilarity(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
	core::fragment::SecondaryStructureOP query_prediction, std::string prediction_name,
	core::Size sequence_length, utility::vector1<core::Size> & frag_sizes, core::Size longest_vall_chunk) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
	"SecondarySimilarity") , prediction_name_(prediction_name) {

	query_len_ = sequence_length;
	query_ss_ = query_prediction;

	runtime_assert( query_prediction->total_residue() == query_len_ );

	for ( core::Size i = 1; i <= query_len_; ++i ) {
		utility::vector1<core::Real> row(longest_vall_chunk);
		scores_.push_back(row);
		utility::vector1<core::Real> prow(3);
		prow[1] = 1 - query_prediction->helix_fraction(i);
		prow[2] = 1 - query_prediction->strand_fraction(i);
		prow[3] = 1 - query_prediction->loop_fraction(i);
		raw_probs_.push_back(prow);
	}

	create_cache(frag_sizes,query_len_,longest_vall_chunk,cache_);
}

} // scores
} // frag_picker
} // protocols


