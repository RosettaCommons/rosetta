// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MCAligner.hh
/// @brief class definition for a class that aligns two Sequence objects using a
/// the Needleman-Wunsch alignment algorithm.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_MCAligner_hh
#define INCLUDED_core_sequence_MCAligner_hh

#include <core/types.hh>

#include <core/sequence/Aligner.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

namespace core {
namespace sequence {

class MCAligner : public Aligner {

public:

	/// @brief constructors
	MCAligner()          : kT_( 1 ) {}
	MCAligner( Real kT ) : kT_( kT ) {}

	/// @brief destructor
	~MCAligner() override = default;

	/// @brief Sets the kT used in the align method. A optimal value of kT means
	/// acceptance of less optimal decisions along the dynamic programming matrix.
	void kT( Real new_kT ) {
		kT_ = new_kT;
	}

	/// @brief Returns the kT used in the align method.
	Real kT() const {
		return kT_;
	}

	/// @brief Align these two Sequences using the given ScoringScheme. Rather
	/// than finding an optimal alignment, MCAligner uses a stochastic algorithm
	/// to generate an alignment at a given kT.
	/// @details The Needleman-Wunsch algorithm uses dynamic programming to
	/// generate an optimal alignment between two scoring sequences under a given
	/// scoring scheme.  Rather than making the best decision at each element of
	/// the dynamic programming matrix, MCAligner makes a stochastic decision
	/// between introducing a gap in seq_y, introducing a gap in seq_x, or aligning two
	/// characters between the two sequences. The decision is made by transforming the
	/// scores for each of the three possible decisions into probabilities using
	/// Boltzmann weighting of the scores for each possibility at a given kT. The kT
	/// is stored as a member variable of the MCAligner class, and accessor methods
	/// are provided above.


	SequenceAlignment align(
		SequenceOP seq_y,
		SequenceOP seq_x,
		ScoringSchemeOP ss
	) override;

private:
	Real kT_;
}; // class MCAligner

} // sequence
} // core

#endif
