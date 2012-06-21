// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

//STL headers

#include <map>

namespace protocols {
namespace jd2 {

protocols::jd2::JobInputter::~JobInputter(){}

///@brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
void JobInputter::load_pose_into_job( core::pose::Pose const & pose, JobOP job ){
	job->inner_job_nonconst()->set_pose( core::pose::PoseCOP( new core::pose::Pose(pose) ) );
}

///@brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
void JobInputter::load_pose_into_job( core::pose::PoseCOP pose, JobOP job ){
	job->inner_job_nonconst()->set_pose( pose );
}

///@brief this code is here to restrict the use of inner_job_nonconst (this class is a friend class and can do it)
core::Size JobInputter::get_nstruct( ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ run::shuffle ]() ) {
		return option[ out::shuffle_nstruct ]();
	} else {
		return option[ out::nstruct ]();
	}
}

core::pose::Pose JobInputter::get_input_structure_from_cache(std::string const & input_tag) const
{
	std::map<std::string, core::pose::Pose>::const_iterator it(input_structure_cache_.find(input_tag));

	//If you hit this assert it's because you didn't check the value of is_input_structure_in_cache first
	assert(it != input_structure_cache_.end());
	return it->second;

}


bool JobInputter::is_input_structure_in_cache(std::string const & input_tag) const
{
	std::map<std::string, core::pose::Pose>::const_iterator it(input_structure_cache_.find(input_tag));
	if(it != input_structure_cache_.end())
	{
		return true;
	}else
	{
		return false;
	}
}


void JobInputter::insert_input_structure_into_cache(std::string const & input_tag, core::pose::Pose const & pose)
{
	std::map<std::string, core::pose::Pose>::const_iterator it(input_structure_cache_.find(input_tag));

	//If you hit this assert it's because you didn't check the value of is_input_structure_in_cache first
	assert(it == input_structure_cache_.end());

	input_structure_cache_.insert(std::make_pair(input_tag,pose));
}


} // jd2
} // protocols
