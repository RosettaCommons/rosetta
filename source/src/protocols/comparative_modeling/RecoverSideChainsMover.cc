// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief simple mover that applies another Mover, and then recovers the input
/// sidechains from the input Pose.
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/SequenceMapping.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>
#include <protocols/comparative_modeling/RecoverSideChainsMover.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

RecoverSideChainsMover::RecoverSideChainsMover(
	Mover & mover
) :
	mover_( mover )
{}

void RecoverSideChainsMover::apply( core::pose::Pose & pose ) {
	using core::Size;
	using namespace core::id;
	core::pose::Pose const original_pose( pose );

	// apply mover
	mover_.apply( pose );

	SequenceMapping map(
		SequenceMapping::identity( original_pose.total_residue() )
	);
	StealSideChainsMover sc_mover( original_pose, map );
	sc_mover.apply( pose );
} // apply

std::string
RecoverSideChainsMover::get_name() const {
	return "RecoverSideChainsMover";
}

} // comparative_modeling
} // protocols
