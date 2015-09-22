// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/bulge/BulgeApplyMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/bulge/BulgeApplyMover.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.bulge.BulgeApplyMover" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace bulge {

//Constructor
BulgeApplyMover::BulgeApplyMover( Size const moving_res ):
	moving_res_( moving_res )
{}

//Destructor
BulgeApplyMover::~BulgeApplyMover()
{}

////////////////////////////////////////////////////////////////////////////////////////
void
BulgeApplyMover::apply( core::pose::Pose & pose ) {
	Pose const pose_save = pose;
	runtime_assert( !is_virtual_base( pose.residue( moving_res_ ) ) );
	core::pose::rna::apply_virtual_rna_residue_variant_type( pose, moving_res_, true );
	protocols::stepwise::modeler::map_constraints_from_original_pose( pose_save, pose );
}


} //bulge
} //rna
} //modeler
} //stepwise
} //protocols
