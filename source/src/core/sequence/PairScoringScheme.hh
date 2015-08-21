// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PairScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Scoring scheme based on comparing sequence columns based on
/// the additive combination of several ScoringScheme objects.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_PairScoringScheme_hh
#define INCLUDED_core_sequence_PairScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/PairScoringScheme.fwd.hh>


#include <utility/vector1_bool.hh>


namespace core {
namespace sequence {

class PairScoringScheme : public ScoringScheme {

public:
	PairScoringScheme() {
		gap_open  ( -4 );
		gap_extend( -1 );
		type("Pair");
	}

	/// @brief ctor
	PairScoringScheme(
		Real gap_open_in,
		Real gap_extend_in,
		utility::vector1< utility::vector1< core::Real > > pairs
	)
	: pairs_( pairs )
	{
		gap_open( gap_open_in );
		gap_extend( gap_extend_in );
		type("Pair");
	}

	ScoringSchemeOP clone() const {
		// maybe clone the scoring_schemes() if object re-use ever causes a weird
		// problem.
		return ScoringSchemeOP( new PairScoringScheme(
			gap_open(),
			gap_extend(),
			pairs()
			) );
	}

	utility::vector1< utility::vector1< core::Real > > pairs() const {
		return pairs_;
	}

	void add_scored_pair(
		core::Size const res1,
		core::Size const res2,
		core::Real const score
	) {
		if ( pairs_.size() <= res1 ) {
			pairs_.resize( res1 );
		}

		if ( pairs_[res1].size() <= res2 ) {
			pairs_[res1].resize( res2 );
		}

		pairs_[res1][res2] = score;
	}

	/// @brief dtor
	virtual ~PairScoringScheme() {}

	virtual void read_from_file( utility::file::FileName const & fn );

	virtual Real score(
		SequenceOP seq1,
		SequenceOP seq2,
		core::Size pos1,
		core::Size pos2
	);

private:
	utility::vector1< utility::vector1< core::Real > > pairs_;
}; // class PairScoringScheme

} // sequence
} // core

#endif
