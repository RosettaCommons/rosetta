// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SimpleScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Scoring scheme based on comparing single characters from two sequences
/// with three types of score:
/// - score for a match (characters are equal)
/// - score for a mismatch (characters are not equal)
/// - affine gap penalties of the form penalty = A + Bk
/// @author James Thompson

#ifndef INCLUDED_core_sequence_SimpleScoringScheme_hh
#define INCLUDED_core_sequence_SimpleScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>

#include <utility/file/FileName.fwd.hh>

namespace core {
namespace sequence {

class SimpleScoringScheme : public ScoringScheme {

public:

	/// @brief ctor
	SimpleScoringScheme(
		Real match_score    =  4,
		Real mismatch_score =  1,
		Real gap_open_in    = -4,
		Real gap_extend_in  = -1
	) :
		match_score_( match_score ),
		mismatch_score_( mismatch_score )
	{
		gap_open( gap_open_in );
		gap_extend( gap_extend_in );
		type("Simple");
	}

	ScoringSchemeOP clone() const {
		return ScoringSchemeOP( new SimpleScoringScheme(
			match_score(),
			mismatch_score(),
			gap_open(),
			gap_extend()
		) );
	}

	/// @brief dtor
	virtual ~SimpleScoringScheme() {}

	virtual void read_from_file( utility::file::FileName const & fn );

	Real match_score() const;

	Real mismatch_score() const;

	virtual Real score( SequenceOP seq1, SequenceOP seq2, core::Size pos1, core::Size pos2 );

private:
	Real match_score_;
	Real mismatch_score_;
}; // class SimpleScoringScheme

} // sequence
} // core

#endif
