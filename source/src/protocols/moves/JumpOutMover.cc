// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file JumpOutMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/JumpOutMover.fwd.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// tracer
#include <basic/Tracer.hh>

#include <utility>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;
// C++ Headers

// ObjexxFCL Headers
//#include <ObjexxFCL/string.functions.hh>

namespace protocols {
namespace moves {

/// JumpOutMover

JumpOutMover::JumpOutMover() : Mover("JumpOutMover"), tolerance_( 100000.0 ) {}


JumpOutMover::JumpOutMover(
	MoverOP first_mover_in,
	MoverOP second_mover_in,
	core::scoring::ScoreFunctionCOP scorefxn_in,
	core::Real tolerance_in
) : Mover("JumpOutMover"),
	first_mover_(std::move(first_mover_in)),
	second_mover_(std::move(second_mover_in)),
	scorefxn_(std::move(scorefxn_in)),
	tolerance_(tolerance_in)
{}

JumpOutMover::~JumpOutMover() = default;

void
JumpOutMover::apply( core::pose::Pose & pose ) {
	using core::scoring::total_score;

	core::Real const initial_score ( pose.energies().total_energy() );
	first_mover_->apply( pose );
	(*scorefxn_)(pose);

	core::Real const move_score( pose.energies().total_energy() );
	if ( (move_score - initial_score) < tolerance_ ) {
		second_mover_->apply( pose );

		// why isn't this done in the constructor?
		type( first_mover_->type() + std::string("+") + second_mover_->type() );
	}
} // apply

std::string
JumpOutMover::get_name() const {
	return "JumpOutMover";
}


} // moves
} // protocols

