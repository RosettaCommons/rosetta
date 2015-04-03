// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetupFromCommandLine.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_PoseSetupFromCommandLine_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_PoseSetupFromCommandLine_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh> // might remove after we move apply_chi_cst out.

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

	void apply_chi_cst( core::pose::Pose & pose, core::pose::Pose const & ref_pose );

	std::string get_working_directory();

	utility::vector1< core::Size >	get_fixed_res( core::Size const nres );

	bool is_nonempty_input_silent_file( std::string const input_silent_file, std::string const exit_key_string );

	utility::vector1< core::Size > get_input_res( core::Size const nres, std::string const pose_num );

	utility::vector1< std::string >	get_silent_file_tags();

	core::scoring::ScoreFunctionOP create_scorefxn();

	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP setup_simple_full_length_rna_working_parameters();

	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP setup_rna_working_parameters( bool check_for_previously_closed_cutpoint_with_input_pose = false );

	void setup_copy_DOF_input( StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup );

	StepWiseRNA_PoseSetupOP setup_pose_setup_class( stepwise::modeler::working_parameters::StepWiseWorkingParametersOP & working_parameters, bool const copy_DOF = true );

	// Undefined, commenting out to fix PyRosetta build  void check_if_silent_file_exists();

	bool
	get_tag_and_silent_file_for_struct( std::string & swa_silent_file,
																			std::string & out_tag,
																			core::Size const & n,
																			bool const & multiple_shots,
																			std::string const & silent_file );
	void
	ensure_directory_for_out_silent_file_exists();


} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
