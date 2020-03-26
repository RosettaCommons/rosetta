// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/TorsionSetMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/simple_moves/TorsionSetMover.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.simple_moves.TorsionSetMover" );

using namespace core;

namespace protocols {
namespace simple_moves {

/// @brief Constructor.
TorsionSetMover::TorsionSetMover(
	utility::vector1< id::TorsionID >  torsion_ids,
	utility::vector1< Real > torsion_values
):
	torsion_ids_(std::move( torsion_ids )),
	torsion_values_(std::move( torsion_values ))
{
	runtime_assert( torsion_ids_.size() == torsion_values_.size() );
}

/// @brief Destructor.
TorsionSetMover::~TorsionSetMover() = default;

void
TorsionSetMover::apply( core::pose::Pose & pose ) {
	for ( core::Size n = 1; n <= torsion_ids_.size(); n++ ) {
		pose.set_torsion( torsion_ids_[ n ], torsion_values_[ n ] );
	}
}

/// @brief Clone function: create a copy of this object and return an owning pointer to the copy.
protocols::moves::MoverOP
TorsionSetMover::clone() const {
	return utility::pointer::make_shared< TorsionSetMover >( *this );
}

} //simple_moves
} //protocols
