// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/Job.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Job classes
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/InnerJob.hh>

///Project headers

// AUTO-REMOVED #include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/JobInputter.hh>

#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <utility/exit.hh>

///C++ headers
#include <string>

#include <utility/vector1.hh>

static basic::Tracer TR("protocols.jd2.InnerJob");

namespace protocols {
namespace jd2 {

/////////////////////////////InnerJob/////////////////////////////
InnerJob::InnerJob( std::string const & input_tag, core::Size nstruct_max )
	: input_tag_(input_tag), nstruct_max_(nstruct_max), pose_(NULL), bad_( false )
{
	//commented these out... somehow they don't respond to -mute
	//TR.Trace << "Using InnerJob (base class) for JobDistributor" << std::endl;
}

InnerJob::InnerJob( core::pose::PoseCOP pose, std::string const & input_tag, core::Size nstruct_max )
	: input_tag_(input_tag), nstruct_max_(nstruct_max), pose_( pose ), bad_( false )
{
	//TR.Trace << "Using InnerJob (base class) for JobDistributor" << std::endl;
}

InnerJob::~InnerJob(){}

std::string const & InnerJob::input_tag() const { return input_tag_; }

core::Size InnerJob::nstruct_max() const { return nstruct_max_; }

core::pose::PoseCOP InnerJob::get_pose() const { return pose_; }

void InnerJob::set_pose( core::pose::PoseCOP pose ) { pose_ = pose; }

} // jd2
} // protocols
