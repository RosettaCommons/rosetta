// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/legacy/screener/RNA_AtrRepScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/screener/RNA_AtrRepScreener.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility>

static basic::Tracer TR( "protocols.stepwise.legacy.screener.RNA_AtrRepScreener" );

using namespace protocols::stepwise::modeler::rna::checker;
using namespace core;

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
RNA_AtrRepScreener::RNA_AtrRepScreener( RNA_AtrRepCheckerOP atr_rep_checker,
	pose::Pose & screening_pose ):
	atr_rep_checker_(std::move( atr_rep_checker )),
	screening_pose_( screening_pose ),
	exit_on_fail_( false )
{}

//Destructor
RNA_AtrRepScreener::~RNA_AtrRepScreener() = default;

////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_AtrRepScreener::check_screen(){
	bool const pass_screen = ( atr_rep_checker_->check_screen( screening_pose_ ) );
	if ( !pass_screen && exit_on_fail_ ) exit( 0 ); // this was for debugging.
	return pass_screen;
}

} //screener
} //legacy
} //stepwise
} //protocols
