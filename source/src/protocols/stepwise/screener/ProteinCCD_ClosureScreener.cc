// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ProteinCCD_ClosureScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/ProteinCCD_ClosureScreener.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/simple_moves/TorsionSetMover.hh>
#include <protocols/moves/CompositionMover.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ProteinCCD_ClosureScreener" );

using protocols::simple_moves::TorsionSetMover;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
  ProteinCCD_ClosureScreener::ProteinCCD_ClosureScreener( sampling::protein::loop_close::StepWiseProteinCCD_CloserOP ccd_closer,
																													pose::Pose & screening_pose ):
		SampleApplier( screening_pose ), // sets up pose_
		ccd_closer_( ccd_closer )
	{
		ccd_closer_->init( pose_ );
	}

	//Destructor
	ProteinCCD_ClosureScreener::~ProteinCCD_ClosureScreener()
	{}

	/////////////////////////////////////////
	bool
	ProteinCCD_ClosureScreener::check_screen() {
		ccd_closer_->get_closure_solution( pose_ );
		return ccd_closer_->closed_loop();
	}

	/////////////////////////////////////////
	void
	ProteinCCD_ClosureScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		update_mover->add_mover(  new TorsionSetMover( ccd_closer_->which_torsions(), ccd_closer_->main_chain_torsion_set() ) );
		restore_mover->add_mover( new TorsionSetMover( ccd_closer_->which_torsions(), ccd_closer_->main_chain_torsion_set_save() ) );
	}

} //screener
} //stepwise
} //protocols
