// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/TorsionSetMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/simple_moves/TorsionSetMover.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.TorsionSetMover" );

using namespace core;

namespace protocols {
namespace simple_moves {

//Constructor
TorsionSetMover::TorsionSetMover(
	utility::vector1< id::TorsionID >  torsion_ids,
	utility::vector1< Real > torsion_values ):
	torsion_ids_( torsion_ids ),
	torsion_values_( torsion_values )
{
	runtime_assert( torsion_ids_.size() == torsion_values_.size() );
}

//Destructor
TorsionSetMover::~TorsionSetMover()
{}

void
TorsionSetMover::apply( core::pose::Pose & pose ) {
	for ( Size n = 1; n <= torsion_ids_.size(); n++ ) {
		pose.set_torsion( torsion_ids_[ n ], torsion_values_[ n ] );
	}
}

} //simple_moves
} //protocols
