// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/AlignRMSD_Screener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/AlignRMSD_Screener.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.AlignRMSD_Screener" );
using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	AlignRMSD_Screener::AlignRMSD_Screener( pose::Pose const & align_pose,
																						pose::Pose const & screening_pose,
																						utility::vector1< core::Size > const & moving_res_list,
																						core::Real const rmsd_cutoff,
																						bool const do_screen /* = true */ ):
		align_pose_( align_pose ),
		screening_pose_( screening_pose ),
		moving_res_list_( moving_res_list ),
		rmsd_cutoff_( rmsd_cutoff ),
		do_screen_( do_screen ),
		pass_count_( 0 ) // only incremented if do_screen true, and screen passes.
	{
		runtime_assert( moving_res_list_.size() > 0 );
		pose_aligner_ = modeler::align::StepWisePoseAlignerOP( new modeler::align::StepWisePoseAligner( align_pose_ ) );
		pose_aligner_->set_user_defined_calc_rms_res( moving_res_list_ );
		// this is extra -- I was originally doings checks of alignment based on RMSD screens.
		//		pose_aligner_->set_root_partition_res( modeler::figure_out_root_partition_res( screening_pose_, moving_res_list ) );
		pose_aligner_->initialize( screening_pose );
	}

	//Destructor
	AlignRMSD_Screener::~AlignRMSD_Screener()
	{}

	bool
	AlignRMSD_Screener::check_screen(){
		if ( !do_screen_ ) return true;

		core::Real const rmsd = pose_aligner_->get_rmsd_no_superimpose( screening_pose_,
																																		false /* check alignment */ );
		if ( rmsd > rmsd_cutoff_ ) return false;

		pass_count_++;
		return true;
	}

} //screener
} //stepwise
} //protocols
