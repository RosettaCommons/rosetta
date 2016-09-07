// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Aligner.hh
/// @brief class definition for a class that aligns two Sequence objects using an
/// dynamic programming and an arbitrary scoring scheme.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_Aligner_hh
#define INCLUDED_core_sequence_Aligner_hh

#include <core/types.hh>

#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/sequence/SequenceAlignment.fwd.hh>

namespace core {
namespace sequence {

class Aligner : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor
	Aligner() {}

	/// @brief dtor
	~Aligner() override = default;

	virtual
	SequenceAlignment align(
		SequenceOP seq_y,
		SequenceOP seq_x,
		ScoringSchemeOP ss
	) = 0;

	void validate_input(
		SequenceOP seq_y,
		SequenceOP seq_x,
		ScoringSchemeOP ss
	);

	SequenceAlignment traceback(
		SequenceOP seq_x,
		SequenceOP seq_y,
		DP_Matrix matrix,
		CellOP start
	);

}; // class Aligner

} // sequence
} // core

#endif
