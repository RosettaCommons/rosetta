// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ChemicalShiftScoringScheme.hh
/// @author James Thompson

#ifndef INCLUDED_core_sequence_ChemicalShiftScoringScheme_hh
#define INCLUDED_core_sequence_ChemicalShiftScoringScheme_hh

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>

namespace core {
namespace sequence {

class ChemicalShiftScoringScheme : public ScoringScheme {
public:

	/// @brief ctor
	ChemicalShiftScoringScheme(
		Real open   = -4,
		Real extend = -1
	)
	{
		gap_open  ( open );
		gap_extend( extend );
		type("ChemicalShift");
	}

	ScoringSchemeOP clone() const {
		return ScoringSchemeOP( new ChemicalShiftScoringScheme(
			gap_open(),
			gap_extend()
			) );
	}

	/// @brief dtor
	virtual ~ChemicalShiftScoringScheme() {}

	virtual Real score( SequenceOP seq1, SequenceOP seq2, Size pos1, Size pos2 );
}; // class ChemicalShiftScoringScheme

} // sequence
} // core

#endif
