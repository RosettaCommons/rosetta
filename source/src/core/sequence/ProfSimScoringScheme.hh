// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProfSimScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Simply based on comparing single profiles from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_ProfSimScoringScheme_hh
#define INCLUDED_core_sequence_ProfSimScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace sequence {

class ProfSimScoringScheme : public ScoringScheme {
public:

	/// @brief constructor
	ProfSimScoringScheme(
		Real open   = -4,
		Real extend = -1
	)
	{
		gap_open  ( open );
		gap_extend( extend );
		initialize_parameters();
		type("ProfSim");
	}

	/// @brief destructor
	virtual ~ProfSimScoringScheme() {}

	/// @brief Initialize log-probabilities of occurence for each amino acid.
	void initialize_parameters();

	/// @brief returns owning pointer to a new object with a deep copy of this
	/// object's values.
	virtual ScoringSchemeOP clone() const {
		return ScoringSchemeOP( new ProfSimScoringScheme(
			gap_open(),
			gap_extend()
		) );
	}

	/// @brief ProfSim profile-profile similarity metric based on information theory.
	/// Published by Yona and Levitt in JMB, 2002 in a paper titled "Within the
	/// Twilight Zone: A Sensitive Profile-Profile Comparison Tool Based on
	/// Information Theory."
	/// @details The basic idea for this score is that it incorporates both
	/// divergence of probability distributions at each position and the significance
	/// of that divergence in order to construct a position-specific
	/// profile-profile score. The divergence score is the J-S divergence between
	/// the two probability distributions at probability position, and the
	/// significance score is the J-S divergence between:
	/// 1. the average of the two probability distributions at this position
	/// 2. a prior probability distribution over all allowed characters at this
	/// position.
	/// J-S divergence between two distributions is defined as:
	/// D( p1, p2 ) = 0.5 * sum( p1[i] * log( p1[i] / p2[i] ) ) +
	///               0.5 * sum( p2[i] * log( p2[i] / p1[i] ) )
	virtual Real score(
		SequenceOP seq1,
		SequenceOP seq2,
		Size pos1,
		Size pos2
	);
private:
	utility::vector1< core::Real > prior_probs_;
}; // class ProfSimScoringScheme

} // sequence
} // core

#endif
