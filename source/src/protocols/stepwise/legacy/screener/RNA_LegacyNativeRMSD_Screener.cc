// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_LegacyNativeRMSD_Screener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/RNA_LegacyNativeRMSD_Screener.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.legacy.screener.RNA_LegacyNativeRMSD_Screener" );

using namespace protocols::stepwise::sampling::rna;

//////////////////////////////////////////////////////////////////////////
// currently RNA-centric
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	//Constructor
	RNA_LegacyNativeRMSD_Screener::RNA_LegacyNativeRMSD_Screener( pose::Pose const & native_pose,
																						pose::Pose & screening_pose,
																						StepWiseWorkingParametersCOP working_parameters,
																						Real const native_screen_rmsd_cutoff,
																						bool const do_screen /* = true */ ):
		native_pose_( native_pose ),
		screening_pose_( screening_pose ),
		working_parameters_( working_parameters ),
		native_screen_rmsd_cutoff_( native_screen_rmsd_cutoff ),
		do_screen_( do_screen ),
		pass_count_( 0 ) // only incremented if do_screen true, and screen passes.
	{}

	//Destructor
	RNA_LegacyNativeRMSD_Screener::~RNA_LegacyNativeRMSD_Screener()
	{}

	bool
	RNA_LegacyNativeRMSD_Screener::check_screen(){
		if ( !do_screen_ ) return true;

		if ( !screening_pose_.residue_type( working_parameters_->actually_moving_res() ).is_RNA() ) return true;
		Real const suite_rmsd_ = suite_rmsd( native_pose_, screening_pose_,
																				 working_parameters_->actually_moving_res(), working_parameters_->is_prepend(), true );
		if ( suite_rmsd_ > native_screen_rmsd_cutoff_ ) return false;

		if ( rmsd_over_residue_list( native_pose_, screening_pose_, working_parameters_, true ) > native_screen_rmsd_cutoff_ ) return false; //Oct 14, 2010

		pass_count_++;
		return true;
	}


} //screener
} //legacy
} //stepwise
} //protocols
