// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <protocols/fibril/SetupForFibrilMover.hh>

// Package headers
#include <core/pose/symmetry/util.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/fibril/fibril_util.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace fibril {

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.symmetry.SetupForFibrilMover" );

SetupForFibrilMover::SetupForFibrilMover()
: Mover("SetupForFibrilMover") {}

SetupForFibrilMover::~SetupForFibrilMover()= default;

void
SetupForFibrilMover::apply( core::pose::Pose & pose )
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;
	protocols::fibril::make_symmetric_fibril( pose );
	assert( core::pose::symmetry::is_symmetric( pose ) );
}

std::string
SetupForFibrilMover::get_name() const {
	return "SetupForFibrilMover";
}

void
SetupForFibrilMover::align(
	core::pose::Pose & pose,
	core::pose::Pose & monomer_pose,
	protocols::loops::Loops core,
	protocols::loops::Loops ref_core
)
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;
	std::cout<<"align: from "<<core<<" to " <<ref_core<<std::endl;
	protocols::fibril::superimpose_pose_on_subset_bb( pose, monomer_pose, core, ref_core );
}

} // fibril
} // protocols
