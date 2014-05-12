// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosableGeometryScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/screener/RNA_ChainClosableGeometryScreener.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.RNA_ChainClosableGeometryScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	RNA_ChainClosableGeometryScreener::RNA_ChainClosableGeometryScreener( sampling::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
																							  pose::PoseOP screening_pose,
																								bool const finer_sampling_at_chain_closure ):
		StepWiseResiduePairScreener( chain_closable_geometry_checker->five_prime_chain_break_res(),
																 chain_closable_geometry_checker->three_prime_chain_break_res() ),
		chain_closable_geometry_checker_( chain_closable_geometry_checker ),
		screening_pose_( screening_pose ),
		finer_sampling_at_chain_closure_( finer_sampling_at_chain_closure )
	{}


	//Destructor
	RNA_ChainClosableGeometryScreener::~RNA_ChainClosableGeometryScreener()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////
	///////////////RNA_Chain_break_screening -- distance cut                     /////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChainClosableGeometryScreener::check_screen() {
		bool const ok = chain_closable_geometry_checker_->check_screen( *screening_pose_, finer_sampling_at_chain_closure_ );
		//		TR << "RES " << chain_closable_geometry_checker_->dist_squared() << " " << ok << std::endl;
		return ok;
	}


} //screener
} //stepwise
} //protocols
