// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/JumpOutMover.hh
/// @brief JumpOutMover
/// @author

#ifndef INCLUDED_protocols_moves_JumpOutMover_hh
#define INCLUDED_protocols_moves_JumpOutMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/JumpOutMover.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {


class JumpOutMover : public Mover {
public:
	// default constructor
	JumpOutMover();

	JumpOutMover(
		MoverOP first_mover_in,
		MoverOP second_mover_in,
		core::scoring::ScoreFunctionCOP scorefxn_in,
		core::Real tolerance_in
	);

	~JumpOutMover();

	/// @brief Applies a move, and conditionally applies a second move.
	/// @details This mover applies the first move, and checks if the score
	/// difference between the initial score and the score after the first move is
	/// within a tolerance. If the score after the first move is within the
	/// tolerance, the second move is applied.
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	MoverOP first_mover_, second_mover_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Real tolerance_;
};


} // moves
} // protocols


#endif
