// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/legacy/ProteinAtrRepScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/legacy/ProteinAtrRepScreener.hh>
#include <protocols/stepwise/sampling/protein/checker/ProteinAtrRepChecker.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ProteinAtrRepScreener" );

////////////////////////////////////////////////////////////////////////////////////////////
//
// This is going to be deprecated soon, in favor of PartitionContactScreener.
//
////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	ProteinAtrRepScreener::ProteinAtrRepScreener( pose::Pose & pose_atr_rep_screen,
																								sampling::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker ):
		SampleApplier( pose_atr_rep_screen ),
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
} //stepwise
} //protocols
