// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobInputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Interface base class JobInputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/JobInputter.hh>

///Project headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>


///Utility headers

///C++ headers

namespace protocols {
namespace jd2 {

protocols::jd2::JobInputter::~JobInputter(){}

/// @brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
void JobInputter::load_pose_into_job( core::pose::Pose const & pose, JobOP job ){
	job->inner_job_nonconst()->set_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) );
}

/// @brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
void JobInputter::load_pose_into_job( core::pose::PoseCOP pose, JobOP job ){
	job->inner_job_nonconst()->set_pose( pose );
}

/// @brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
core::Size JobInputter::get_nstruct( ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ run::shuffle ]() ) {
		return option[ out::shuffle_nstruct ]();
	} else {
		return option[ out::nstruct ]();
	}
}

/// @brief This function is only called by certain JobInputters to update the jobs list after it has already been created.
/// @details An example case would be the LargeNstructJobInputter, which uses this function to load additional jobs after the first N have started to come back.
void JobInputter::update_jobs_list( JobsContainerOP /*jobs*/ ) {
	//Do nothing by default, unless implemented in a derived class.
	return;
}

std::string
JobInputter::job_inputter_input_source_to_string(
	JobInputterInputSource::Enum source
) {
	switch(source) {
	case JobInputterInputSource::NONE : return "None";
	case JobInputterInputSource::UNKNOWN : return "Uknown";
	case JobInputterInputSource::POSE : return "Pose";
	case JobInputterInputSource::SILENT_FILE : return "SilentFile";
	case JobInputterInputSource::PDB_FILE : return "PdbFile";
	case JobInputterInputSource::ATOM_TREE_FILE : return "AtomTreeFile";
	case JobInputterInputSource::DATABASE : return "Database";
	case JobInputterInputSource::MAKE_ROT_LIB : return "MakeRotLib";
	case JobInputterInputSource::RESOURCE_MANAGED_JOB : return "ResourceManagerJob";
	case JobInputterInputSource::SCREENING_FILE : return "ScreeningJobInputter";
	default :
		utility_exit_with_message("Unrecognized JobInputterInputSource");
	}
}

} // jd2
} // protocols
