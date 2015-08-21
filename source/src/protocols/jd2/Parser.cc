// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/Parser.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Interface base class Parser
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/Parser.hh>

///Package headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <sstream>

namespace protocols {
namespace jd2 {

protocols::jd2::Parser::~Parser(){}

/// @details the impetus for this function is that Parser is a friend of the
/// InnerJob class and can modify the pose - actual mover generation/pose
/// updating is handled by derived classes. Now, if the generate_mover_from_pose
/// function modifies the input pose as part of the APPLY_TO_POSE block, then
/// the modified Pose will be stored in the input Job, but only so long as
/// allow_job_update is true.
void
protocols::jd2::Parser::generate_mover_from_job(
	JobOP job,
	core::pose::Pose & pose,
	protocols::moves::MoverOP & mover,
	bool new_input,
	bool allow_job_update,
	bool guarantee_new_mover
){


	std::stringstream err_msg;
	err_msg
		<< "Attempting to initiate job distribution for "
		<< "input job '" << job->input_tag() << "', "
		<< "but the generated pose has no residues."
		<< "make sure you specified a valid input PDB, silent file "
		<< "or database.";

	runtime_assert_string_msg( pose.total_residue() > 0, err_msg.str() );

	// generate_mover_from_pose returns true if there was a pose change
	// (i.e. NOT if only the mover changed)
	if ( generate_mover_from_pose( job, pose, mover, new_input,
			"" /*empty xml_fname, this means go to options system*/,
			guarantee_new_mover ) ) {

		if ( allow_job_update && new_input ) {
			// Store the modified pose into the job.  This is done only if new_input is true even if
			// generate_mover_from_pose modified the input pose and if allow_job_update is true.

			job->inner_job_nonconst()->set_pose( pose.clone() );
		}

	}

	return;
}

}//jd2
}//protocols
