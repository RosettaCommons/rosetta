// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SimpleScoringScheme.cc
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Scoring scheme based on comparing single characters from two sequences
/// with three types of score:
/// - score for a match (characters are equal)
/// - score for a mismatch (characters are not equal)
/// - affine gap penalties of the form penalty = A + Bk
/// @author James Thompson

#include <core/types.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SimpleScoringScheme.hh>

#include <utility/exit.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

void SimpleScoringScheme::read_from_file( utility::file::FileName const & /*fn*/ ) {
	utility_exit_with_message( "SimpleScoringScheme::read_from_file method stubbed out!");
}

Real SimpleScoringScheme::match_score() const {
	return match_score_;
}

Real SimpleScoringScheme::mismatch_score() const {
	return mismatch_score_;
}

Real SimpleScoringScheme::score(
	SequenceOP seq1, SequenceOP seq2,
	core::Size pos1, core::Size pos2
) {
	runtime_assert( pos1 <= seq1->length() );
	runtime_assert( pos2 <= seq2->length() );

	if ( (*seq1)[pos1] == (*seq2)[pos2] ) return match_score_;
	else                                  return mismatch_score_;
}

} // sequence
} // core
