// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/AtrRepScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/AtrRepScreener.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.stepwise.screener.AtrRepScreener" );

using namespace protocols::stepwise::sampling::rna::checker;
using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	AtrRepScreener::AtrRepScreener( AtrRepCheckerOP atr_rep_checker,
																	pose::Pose & screening_pose ):
		atr_rep_checker_( atr_rep_checker ),
		screening_pose_( screening_pose ),
		exit_on_fail_( false )
	{}

	//Destructor
	AtrRepScreener::~AtrRepScreener()
	{}

	////////////////////////////////////////////////////////////////////////////////////////
	bool
	AtrRepScreener::check_screen(){
		bool const pass_screen = ( atr_rep_checker_->check_screen( screening_pose_ ) );
		//		TR << pass_screen
		//			 << " " << atr_rep_checker_->delta_atr_score() << " " << atr_rep_checker_->delta_rep_score()
		//			 << " " << atr_rep_checker_->base_atr_score() << " " << atr_rep_checker_->base_rep_score()  << std::endl;
		//		screening_pose_.dump_pdb( "BLAH"+ObjexxFCL::string_of( exit_on_fail_ )+".pdb" );
		if ( !pass_screen && exit_on_fail_ ) exit( 0 );
		return pass_screen;
	}

} //screener
} //stepwise
} //protocols
