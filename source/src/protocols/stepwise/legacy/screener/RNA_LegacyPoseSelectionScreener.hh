// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_LegacyPoseSelectionScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_LegacyPoseSelectionScreener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_LegacyPoseSelectionScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/TagDefinition.fwd.hh>
#include <protocols/stepwise/screener/RNA_LegacyPoseSelectionScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>

using namespace protocols::stepwise::sampling::working_parameters;
using namespace protocols::stepwise::sampling::rna;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	class RNA_LegacyPoseSelectionScreener: public stepwise::screener::StepWiseScreener {

	public:

		//constructor
		RNA_LegacyPoseSelectionScreener(  StepWiseRNA_PoseSelectionOP pose_selection,
														Pose const & pose,
														TagDefinitionOP tag_definition,
														bool const output_to_silent_file,
														std::string const silent_file,
														PoseCOP native_pose,
														StepWiseWorkingParametersCOP working_parameters );

		//destructor
		~RNA_LegacyPoseSelectionScreener();

	public:

		std::string
		name() const { return "RNA_LegacyPoseSelectionScreener"; }

		stepwise::screener::StepWiseScreenerType
		type() const { return POSE_SELECTION; }

		bool
		check_screen();

	private:

		StepWiseRNA_PoseSelectionOP pose_selection_;
		Pose const & pose_;
		TagDefinitionOP tag_definition_;
		bool const output_to_silent_file_;
		std::string const silent_file_;
		PoseCOP native_pose_;
		StepWiseWorkingParametersCOP working_parameters_;

	};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
