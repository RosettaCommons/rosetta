// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/sugar/SugarVirtualizeMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/sugar/SugarVirtualizeMover.hh>
#include <core/chemical/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.sugar.SugarVirtualizeMover" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {


//Constructor
SugarVirtualizeMover::SugarVirtualizeMover( Size const moving_res ):
	moving_res_( moving_res )
{}

//Destructor
SugarVirtualizeMover::~SugarVirtualizeMover()
{}

void
SugarVirtualizeMover::apply( pose::Pose & pose ){
	if ( !pose.residue_type( moving_res_ ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) {
		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, moving_res_ ); //May 31, 2010
	}
	if ( !pose.residue_type( moving_res_ ).has_variant_type( core::chemical::VIRTUAL_O2PRIME_HYDROGEN ) ) {
		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, moving_res_ );
	}
}

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols
