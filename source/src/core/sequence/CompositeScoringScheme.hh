// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CompositeScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Scoring scheme based on comparing sequence columns based on
/// the additive combination of several ScoringScheme objects.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_CompositeScoringScheme_hh
#define INCLUDED_core_sequence_CompositeScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/CompositeScoringScheme.fwd.hh>

#include <utility/vector1_bool.hh>

namespace core {
namespace sequence {

class CompositeScoringScheme : public ScoringScheme {

public:
	CompositeScoringScheme() {
		gap_open  ( -4 );
		gap_extend( -1 );
		type("Composite");
	}

	/// @brief ctor
	CompositeScoringScheme(
		Real gap_open_in,
		Real gap_extend_in,
		utility::vector1< ScoringSchemeOP > schemes
	) :
		scoring_schemes_( schemes )
	{
		gap_open( gap_open_in );
		gap_extend( gap_extend_in );
		type("Composite");
	}

	ScoringSchemeOP clone() const {
		// maybe clone the scoring_schemes() if object re-use ever causes a weird
		// problem.
		return ScoringSchemeOP( new CompositeScoringScheme(
			gap_open(),
			gap_extend(),
			scoring_schemes()
			) );
	}

	utility::vector1< ScoringSchemeOP > scoring_schemes() const {
		return scoring_schemes_;
	}

	Size count() const {
		return scoring_schemes_.size();
	}

	void add_scoring_scheme( ScoringSchemeOP scheme );

	/// @brief dtor
	virtual ~CompositeScoringScheme() {}

	virtual void read_from_file( utility::file::FileName const & fn );

	virtual Real score( SequenceOP seq1, SequenceOP seq2, core::Size pos1, core::Size pos2 );

private:
	utility::vector1< ScoringSchemeOP > scoring_schemes_;
}; // class CompositeScoringScheme

} // sequence
} // core

#endif
