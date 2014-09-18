// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/screener/ProteinAtrRepScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/screener/ProteinAtrRepScreener.hh>
#include <protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.screener.ProteinAtrRepScreener" );

////////////////////////////////////////////////////////////////////////////////////////////
//
// This is going to be deprecated soon, in favor of PartitionContactScreener.
//
////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	//Constructor
	ProteinAtrRepScreener::ProteinAtrRepScreener( pose::Pose & pose_atr_rep_screen,
																								modeler::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker ):
		stepwise::screener::SampleApplier( pose_atr_rep_screen ),
		atr_rep_checker_( atr_rep_checker )
	{}

	//Destructor
	ProteinAtrRepScreener::~ProteinAtrRepScreener()
	{}

	bool
	ProteinAtrRepScreener::check_screen(){
		bool const pass_screen = ( atr_rep_checker_->check_screen( pose_ ) );
		return pass_screen;
	}

} //screener
} //legacy
} //stepwise
} //protocols
