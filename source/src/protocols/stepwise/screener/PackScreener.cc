// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/PackScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/PackScreener.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/packer/SideChainCopier.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.PackScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	PackScreener::PackScreener( pose::Pose & pose,
																						modeler::packer::StepWisePackerOP stepwise_packer ):
		SampleApplier( pose ),
		stepwise_packer_( stepwise_packer )
	{
	}

	//Destructor
	PackScreener::~PackScreener()
	{}

	bool
	PackScreener::check_screen() {
		stepwise_packer_->apply( pose_ );
		return true;
	}

	/////////////////////////////////////////
	void
	PackScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		update_mover->add_mover(  new modeler::packer::SideChainCopier( pose_,
																																		 stepwise_packer_->previous_working_pack_res(),
																																		 stepwise_packer_->pack_o2prime_hydrogens() )  );
		restore_mover->add_mover( 0 ); // original choice.
		//		restore_mover->add_mover( SideChainMover( *pose_original_ ) );
	}


} //screener
} //stepwise
} //protocols
