// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.cc
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScoreBlosum62.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/chemical/AA.hh>

#include <utility/file/FileName.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>

// utils
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;
	using namespace core::chemical;

static thread_local basic::Tracer trProfScoreBlosum62(
		"protocols.frag_picker.scores.ProfileScoreBlosum62");

ProfileScoreBlosum62::~ProfileScoreBlosum62() {}

ProfileScoreBlosum62::ProfileScoreBlosum62(Size priority, Real lowest_acceptable_value, bool use_lowest,
			sequence::SequenceProfileOP query_profile, Size longest_vall_chunk) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
				"ProfileScoreBlosum62")
{
	if( query_profile->size() != query_profile->length() ) {
		utility_exit_with_message("ProfileScoreBlosum62 needs a valid sequence profile.");
	}
	query_profile_ = query_profile;

	for (Size i = 1; i <= query_profile->length(); ++i) {
		utility::vector1<Real> row(longest_vall_chunk);
		scores_.push_back(row);
	}

	sequence::MatrixScoringScheme blosum_reader;

	utility::file::FileName const fn("/work/rvernon/yangshen/adjust_tables/BLOSUM62");

	blosum_reader.read_from_file(fn);

	blosum_matrix_ = blosum_reader.scoring_matrix();

}


void ProfileScoreBlosum62::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;
	Size size_q = query_profile_->length();

	trProfScoreBlosum62.Debug << "caching profile score for " << chunk->get_pdb_id()
			<< " of size " << chunk->size() << std::endl;
	PROF_START( basic::FRAGMENTPICKING_PROFILE_CAHING );
	for (Size i = 1; i <= size_q; ++i) {
		utility::vector1<Real> query_prof_row = query_profile_->prof_row(i);

		for (Size j = 1; j <= chunk->size(); ++j) {

			//AA aa = aa_from_oneletter_code( chunk->at(j)->aa() );

			utility::vector1<Real> tmplt_prof_row = chunk->at(j)->profile();
			Real score = 0.0;
			for (Size k1 = 1; k1 <= 20; k1++) {
				//score -= query_prof_row[k]*blosum_matrix_[k][aa];//*tmplt_prof_row[k];
				for (Size k2 = 1; k2 <= 20; k2++) {

					Real diff = query_prof_row[k1]*sqrt(tmplt_prof_row[k2]);

					Real sign(0.0);

					if ( blosum_matrix_[k1][k2] < 0 ) {
						sign = -1.0;
					} else {
						sign = 1.0;
					}

					Real blosum( pow( std::abs(blosum_matrix_[k1][k2]), 1/4) );

					score -= sign/(1+exp((-10*diff*blosum)+5));

					//score -= query_prof_row[k1]*tmplt_prof_row[k2]*blosum_matrix_[k1][k2];//*tmplt_prof_row[k];;
				}
			}
			scores_[i][j] = score;
		}


		//for (Size j = 1; j <= chunk->size(); ++j) {
		//	utility::vector1<Real> tmplt_prof_row = chunk->at(j)->profile();
		//	Real score = 0.0;
		//	for (Size k = 1; k <= 20; k++){
		//		score += std::abs(tmplt_prof_row[k] - query_prof_row[k]);
		//	}
		//	scores_[i][j] = score;
		//}
	}

	PROF_STOP( basic::FRAGMENTPICKING_PROFILE_CAHING );
	trProfScoreBlosum62.Debug << "precomputed matrix of scores " << scores_.size()
			<< "x" << chunk->size() << std::endl;
}

bool ProfileScoreBlosum62::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	std::string & tmp = f->get_chunk()->chunk_key();

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

bool ProfileScoreBlosum62::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	return cached_score( f, empty_map);

// 	PROF_START( basic::FRAGMENTPICKING_PROFILE_SCORE );

// 	Real totalScore = 0.0;
// 	for (Size i = 1; i <= f->get_length(); i++) {
// 		utility::vector1<Real> query_prof_row = query_profile_->prof_row(f->get_first_index_in_query() + i - 1);
// 		VallChunkOP chunk = f->get_chunk();
// 		utility::vector1<Real> tmplt_prof_row = chunk->at(f->get_first_index_in_vall() + i - 1)->profile();
// 		for (Size k = 1; k <= 20; k++){
// 		    totalScore += std::abs(tmplt_prof_row[k] - query_prof_row[k]);
// 		}
// 	}
// 	totalScore /= (Real) f->get_length();
// 	empty_map->set_score_component(totalScore, id_);
// 	PROF_STOP( basic::FRAGMENTPICKING_PROFILE_SCORE );
// 	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
// 		return false;
// 	return true;
}


FragmentScoringMethodOP MakeProfileScoreBlosum62::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	Size len = picker->get_vall()->get_largest_chunk_size();
	trProfScoreBlosum62 << "Profile scoring method is: Blosum62" << std::endl;
	return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ProfileScoreBlosum62(priority,
			lowest_acceptable_value, use_lowest, picker->get_query_seq(), len) );
}

} //scores
} // frag_picker
} // protocols
