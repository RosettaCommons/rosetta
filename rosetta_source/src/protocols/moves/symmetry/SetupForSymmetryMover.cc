// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/symmetric_docking/SymDockingInitialPerturbation.hh>

// Package headers
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

namespace protocols {
namespace moves {
namespace symmetry {

static basic::Tracer TR("protocols.moves.symmetry.SetupForSymmetryMover");

SetupForSymmetryMover::SetupForSymmetryMover()
	: Mover("SetupForSymmetryMover"), slide_(false) {}

SetupForSymmetryMover::~SetupForSymmetryMover(){}

void
SetupForSymmetryMover::apply( core::pose::Pose & pose )
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;

	core::pose::symmetry::make_symmetric_pose( pose );
	assert( core::pose::symmetry::is_symmetric( pose ) );

	// (Optionally) set rigid-body dofs from file
	//    SymDockingInitialPerturbation's behavior is controlled by flags and does nothing by default
	protocols::moves::MoverOP symdock =
		new protocols::symmetric_docking::SymDockingInitialPerturbation(slide_);
	symdock->apply( pose );
}

std::string
SetupForSymmetryMover::get_name() const {
	return "SetupForSymmetryMover";
}


} // symmetry
} // moves
} // protocols
