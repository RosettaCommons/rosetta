// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_LegacyPoseSelectionScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/RNA_LegacyPoseSelectionScreener.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.legacy.screener.RNA_LegacyPoseSelectionScreener" );

using namespace protocols::stepwise::sampling::rna;
using namespace core::pose;
using namespace core;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	//Constructor
	RNA_LegacyPoseSelectionScreener::RNA_LegacyPoseSelectionScreener( StepWiseRNA_PoseSelectionOP pose_selection,
																								pose::Pose const & pose,
																								TagDefinitionOP tag_definition,
																								bool const output_to_silent_file,
																								std::string const silent_file,
																								PoseCOP native_pose,
																								StepWiseWorkingParametersCOP working_parameters ):
		pose_selection_( pose_selection ),
		pose_( pose ),
		tag_definition_( tag_definition ),
		output_to_silent_file_( output_to_silent_file ),
		silent_file_( silent_file ),
		native_pose_( native_pose ),
		working_parameters_( working_parameters)
	{
	}

	//Destructor
	RNA_LegacyPoseSelectionScreener::~RNA_LegacyPoseSelectionScreener()
	{}

	bool
	RNA_LegacyPoseSelectionScreener::check_screen(){
		// the original reason for this copy was that we might apply a variant, which can produce thread conflicts with graphics.
		// also adding tags -- may want to leave original pose unspoiled.
		Pose selected_pose = pose_;

		std::string const & tag = tag_definition_->tag();
		pose_selection_->pose_selection_by_full_score( selected_pose, tag );

		if ( output_to_silent_file_ ) sampling::rna::output_data( silent_file_, tag, true,
																															 selected_pose, native_pose_, working_parameters_ );

		return true;
	}

} //screener
} //legacy
} //stepwise
} //protocols
