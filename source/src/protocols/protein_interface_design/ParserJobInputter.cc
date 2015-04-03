// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PDBJobInputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class PDBJobInputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/protein_interface_design/ParserJobInputter.hh>
#include <protocols/protein_interface_design/ParserJobInputterCreator.hh>
#include <protocols/protein_interface_design/read_patchdock.hh>
#include <protocols/jd2/Job.hh>

///Project headers

#ifdef WIN32
// required for VS2005 build
#include <core/conformation/Residue.hh>
#endif

///Utility headers
#include <basic/Tracer.hh>

///C++ headers
#include <string>

// option key includes

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.protein_interface_design.ParserJobInputter" );

namespace protocols {
namespace protein_interface_design {

using namespace protocols::jd2;

ParserJobInputter::ParserJobInputter(){
	TR << "Instantiate ParserJobInputter" << std::endl;
}

ParserJobInputter::~ParserJobInputter(){
}

/// @details This function will first see if the pose already exists in the Job.  If not, it will read it into the pose reference, and hand a COP cloned from that pose to the Job. If the pose pre-exists it just copies the COP's pose into it.
void
ParserJobInputter::pose_from_job( core::pose::Pose & pose, JobOP job){
	static protocols::protein_interface_design::PatchdockReader pd_reader;

	TR << "ParserJobInputter::pose_from_job" << std::endl;

/// normally the job would contain the saved pose, but here, the saved pose
/// would be transformed according to patchdock reading. So, we're letting
/// the static patchdock reader take care of saving the pose in memory

	std::string input_tag( job->input_tag() );
	pd_reader.read_poses( pose, input_tag );
	load_pose_into_job(pose, job);
}

//CREATOR SECTION
std::string
ParserJobInputterCreator::keyname() const
{
        return "ParserJobInputter";
}

protocols::jd2::JobInputterOP
ParserJobInputterCreator::create_JobInputter() const {
        return protocols::jd2::JobInputterOP( new ParserJobInputter );
}

}//protein_interface_design
}//protocols
