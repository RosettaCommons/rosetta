// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/bulge/BulgeUnApplyMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/bulge/BulgeUnApplyMover.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.rna.bulge.BulgeUnApplyMover" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace bulge {

	//Constructor
	BulgeUnApplyMover::BulgeUnApplyMover( Size const moving_res ):
		moving_res_( moving_res )
	{}

	//Destructor
	BulgeUnApplyMover::~BulgeUnApplyMover()
	{}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	BulgeUnApplyMover::apply( core::pose::Pose & pose ) {
		runtime_assert( is_virtual_base( pose.residue( moving_res_ ) ) );
		remove_virtual_rna_residue_variant_type( pose, moving_res_ );
	}

} //bulge
} //rna
} //sampling
} //stepwise
} //protocols
