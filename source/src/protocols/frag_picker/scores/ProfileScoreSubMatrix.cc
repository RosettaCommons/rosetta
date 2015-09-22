// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScoreSubMatrix.cc
/// @brief  scores a fragment by substitution matrix (e.g. Blosum62) profile score
/// @author Dan Kulp (dwkulp@gmail.com), based on code from Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScoreSubMatrix.hh>

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

static THREAD_LOCAL basic::Tracer trProfScoreSubMatrix(
	"protocols.frag_picker.scores.ProfileScoreSubMatrix");

ProfileScoreSubMatrix::~ProfileScoreSubMatrix() {}

ProfileScoreSubMatrix::ProfileScoreSubMatrix(Size priority, Real lowest_acceptable_value, bool use_lowest,
	std::string sequence,Size longest_vall_chunk,std::string subMatrixFile) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
	"ProfileScoreSubMatrix")
{
	// Store local copies
	sequence_                = sequence;
	subMatrixFile_           = subMatrixFile;

	// Setup scores to be the proper size, each sequence position has an entry for each chunk position.
	for ( Size i = 1; i <= sequence.size(); ++i ) {
		utility::vector1<Real> row(longest_vall_chunk);
		scores_.push_back(row);
	}

	// Read in the SubMatrix, should be BLOSUM62 format?
	sequence::MatrixScoringScheme sub_matrix_reader;
	sub_matrix_reader.read_from_file(subMatrixFile_);
	sub_matrix_ = sub_matrix_reader.scoring_matrix();


	// Convert score matrix to probability matrix ?

}


void ProfileScoreSubMatrix::do_caching(VallChunkOP chunk) {

	std::string & tmp = chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;
	Size size_q = sequence_.size();

	trProfScoreSubMatrix.Debug << "caching profile score for " << chunk->get_pdb_id()
		<< " of size " << chunk->size() << std::endl;
	PROF_START( basic::FRAGMENTPICKING_PROFILE_CAHING );

	// For each position in sequence
	for ( Size i = 1; i <= size_q; ++i ) {

		AA seqAA = aa_from_oneletter_code( sequence_[i-1] );

		// For each position in chunk
		for ( Size j = 1; j <= chunk->size(); ++j ) {

			AA chunkAA = aa_from_oneletter_code( chunk->at(j)->aa() );

			// For BLOSUM62, a positive score means likely to mutate (F -> Y is 3)
			// Here lower score is better, so invert the substition matrix
			scores_[i][j] = - sub_matrix_[seqAA][chunkAA];
		}
	}

	PROF_STOP( basic::FRAGMENTPICKING_PROFILE_CAHING );
	trProfScoreSubMatrix.Debug << "precomputed matrix of scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}

bool ProfileScoreSubMatrix::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	std::string & tmp = f->get_chunk()->chunk_key();

	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(f->get_chunk());
	}

	Real totalScore = 0.0;
	for ( Size i = 1; i <= f->get_length(); i++ ) {

		// Check sizes of fragment/chunk and scores
		assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		assert(f->get_first_index_in_vall() + i - 1<= scores_[1].size());

		totalScore += scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1];
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);

	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

bool ProfileScoreSubMatrix::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {
	return cached_score( f, empty_map);
}

// MISC = Substitution Matrix File Name ... the reader will spit out an error if the file does not exist.
FragmentScoringMethodOP MakeProfileScoreSubMatrix::make(Size priority, Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string misc) {

	Size len = picker->get_vall()->get_largest_chunk_size();

	trProfScoreSubMatrix << "Profile scoring method is: SubMatrix" << std::endl;

	return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ProfileScoreSubMatrix(
		priority,
		lowest_acceptable_value,
		use_lowest,
		picker->get_query_seq()->sequence(),
		len,
		misc  // is the substituion matrix file location
		) );
}

} //scores
} // frag_picker
} // protocols
