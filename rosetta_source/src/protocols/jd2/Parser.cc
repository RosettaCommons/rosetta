// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

//Auto Headers
#include <utility/exit.hh>


///Utility headers

///C++ headers

namespace protocols {
namespace jd2 {

protocols::jd2::Parser::~Parser(){}

///@details the impetus for this function is that Parser is a friend of the InnerJob class and can modify the pose - actual mover generation/pose updating is handled by derived classes.
void
protocols::jd2::Parser::generate_mover_from_job( JobOP job, protocols::moves::MoverOP & mover, bool new_input ){

	//unpackage job
	core::pose::Pose pose( *(job->get_pose()) );
	runtime_assert(pose.total_residue()); //if this is an empty pose we're in trouble

	//returns true if there was a pose change (NOT if the mover changed)
	if ( generate_mover_from_pose( job, pose, mover, new_input, ""/*xml_fname, this means go to options system*/ ) ){
		//repackage pose into job
		job->inner_job_nonconst()->set_pose( core::pose::PoseCOP( new core::pose::Pose(pose) ) );
	}

	return;
}

}//jd2
}//protocols
