// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.cc
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScoreL1.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>

// utils
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer trProfScoreL1(
	"protocols.frag_picker.scores.ProfileScoreL1");

ProfileScoreL1::~ProfileScoreL1() {}

ProfileScoreL1::ProfileScoreL1(Size priority, Real lowest_acceptable_value, bool use_lowest,
	sequence::SequenceProfileOP query_profile, utility::vector1<Size> & frag_sizes,
	Size longest_vall_chunk) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
	"ProfileScoreL1") {
	if ( query_profile->size() != query_profile->length() ) {
		utility_exit_with_message("ProfileScoreL1 needs a valid sequence profile." );
	}
	query_profile_ = query_profile;

	for ( Size i = 1; i <= query_profile->length(); ++i ) {
		utility::vector1<Real> row(longest_vall_chunk);
		scores_.push_back(row);
	}
	create_cache(frag_sizes,query_profile->length(),longest_vall_chunk,cache_);
	if ( trProfScoreL1.visible() ) {
		trProfScoreL1 << "Created cache for fraglen:";
		for ( Size i=1; i<=cache_.size(); i++ ) {
			if ( cache_[i].size() > 0 ) {
				trProfScoreL1 <<" "<<i;
			}
		}
		trProfScoreL1<<std::endl;
	}
}


void ProfileScoreL1::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;

	do_caching_simple(chunk);
	for ( Size fl=1; fl<=cache_.size(); fl++ ) {
		if ( cache_[fl].size() != 0 ) {
			trProfScoreL1.Trace << "caching profile score for " << chunk->get_pdb_id()
				<< " of size " << chunk->size() << " for fragment size "<<fl<<std::endl;
			rolling_score(scores_,fl,cache_[fl]);
		}
	}
}


void ProfileScoreL1::do_caching_simple(VallChunkOP chunk) {

	Size size_q = query_profile_->length();

	trProfScoreL1.Debug << "caching profile score for " << chunk->get_pdb_id()
		<< " of size " << chunk->size() << std::endl;
	for ( Size i = 1; i <= size_q; ++i ) {
		utility::vector1<Real> query_prof_row = query_profile_->prof_row(i);
		for ( Size j = 1; j <= chunk->size(); ++j ) {
			utility::vector1<Real> tmplt_prof_row = chunk->at(j)->profile();
			Real score = 0.0;
			for ( Size k = 1; k <= 20; k++ ) {
				score += std::abs(tmplt_prof_row[k] - query_prof_row[k]);
			}
			scores_[i][j] = score;
		}
	}

	trProfScoreL1.Debug << "precomputed matrix of L1 scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}

bool ProfileScoreL1::cached_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	std::string tmp = f->get_chunk()->chunk_key();
	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(f->get_chunk());
	}

	Real totalScore = cache_[f->get_length()][f->get_first_index_in_query()][f->get_first_index_in_vall()];

	totalScore /= (Real) f->get_length();

	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

bool ProfileScoreL1::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	Real totalScore = 0.0;
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		utility::vector1<Real> query_prof_row = query_profile_->prof_row(f->get_first_index_in_query() + i - 1);
		VallChunkOP chunk = f->get_chunk();
		utility::vector1<Real> tmplt_prof_row = chunk->at(f->get_first_index_in_vall() + i - 1)->profile();
		for ( Size k = 1; k <= 20; k++ ) {
			totalScore += std::abs(tmplt_prof_row[k] - query_prof_row[k]);
		}
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}


FragmentScoringMethodOP MakeProfileScoreL1::make(Size priority,
	Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	Size len = picker->get_vall()->get_largest_chunk_size();
	trProfScoreL1 << "Profile scoring method is: L1" << std::endl;
	return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ProfileScoreL1(priority,
		lowest_acceptable_value, use_lowest, picker->get_query_seq(),
		picker->frag_sizes_,len) );
}

} //scores
} // frag_picker
} // protocols

/*
bool ProfileScoreL1::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

std::string tmp = f->get_chunk()->chunk_key();
if (tmp.compare(cached_scores_id_) != 0)
do_caching(f->get_chunk());

Real totalScore = 0.0;
for (Size i = 1; i <= f->get_length(); i++) {
assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
assert(f->get_first_index_in_vall()
+ i - 1<= scores_[1].size());
totalScore
+= scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall()
+ i - 1];
}
totalScore /= (Real) f->get_length();
empty_map->set_score_component(totalScore, id_);
if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
return false;
return true;
}
*/
