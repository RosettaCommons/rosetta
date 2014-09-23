// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CompositeScoringScheme.cc
/// @brief class definition for a given scoring scheme for an alignment.
/// @detailed Scoring scheme based on comparing single characters from two sequences
/// with three types of score:
/// - score for a match (characters are equal)
/// - score for a mismatch (characters are not equal)
/// - affine gap penalties of the form penalty = A + Bk
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/CompositeSequence.hh>
#include <core/sequence/CompositeScoringScheme.hh>

#include <utility/exit.hh>

// AUTO-REMOVED #include <iostream>
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {


void CompositeScoringScheme::read_from_file( utility::file::FileName const & /*fn*/ ) {
	utility_exit_with_message( "CompositeScoringScheme::read_from_file method stubbed out!");
}

void CompositeScoringScheme::add_scoring_scheme( ScoringSchemeOP scheme ) {
	scoring_schemes_.push_back( scheme->clone() );
}

Real CompositeScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	core::Size pos1,
	core::Size pos2
) {
	runtime_assert( pos1 <= seq1->length() );
	runtime_assert( pos2 <= seq2->length() );

	Real total_score( 0.0 );
	CompositeSequenceOP cs1( utility::pointer::dynamic_pointer_cast< core::sequence::CompositeSequence > ( seq1 ) );
	CompositeSequenceOP cs2( utility::pointer::dynamic_pointer_cast< core::sequence::CompositeSequence > ( seq2 ) );
	runtime_assert( cs1 != 0 );
	runtime_assert( cs2 != 0 );
	runtime_assert( cs1->n_seqs() == cs2->n_seqs() );
	runtime_assert( count()       == cs2->n_seqs() );

	for ( Size idx = 1; idx <= count(); ++idx ) {
		total_score += scoring_schemes_[idx]->score(
			cs1->seq(idx), cs2->seq(idx), pos1, pos2
		);
	}
	return total_score;
}

} // sequence
} // core
