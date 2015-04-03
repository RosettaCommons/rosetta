// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PairScoringScheme.cc
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/PairScoringScheme.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.fwd.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

void PairScoringScheme::read_from_file( utility::file::FileName const & /*fn*/ ) {
	utility_exit_with_message(
		"PairScoringScheme::read_from_file method stubbed out!"
	);
}

Real PairScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	core::Size pos1,
	core::Size pos2
) {
	runtime_assert( pos1 <= seq1->length() );
	runtime_assert( pos2 <= seq2->length() );

	// looks unsafe, but isn't because of short-circuit operators.
	if ( pos1 <= pairs_.size() && pos2 <= pairs_[pos1].size() ) {
		return pairs_[pos1][pos2];
	} else {
		return 0.0;
	}
}

} // sequence
} // core
