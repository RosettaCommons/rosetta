// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Job classes
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/InnerMultiThreadingJob.hh>

///Project headers
#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>

///C++ headers
#include <string>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.InnerMultiThreadingJob" );

namespace protocols {
namespace jd2 {

InnerMultiThreadingJob::InnerMultiThreadingJob(
	core::pose::PoseCOP pose,
	std::string const & input_tag,
	core::Size nstruct_max,
	utility::vector1< core::sequence::SequenceAlignment > const & alns,
	utility::vector1< core::pose::Pose > const & templates
) :
	InnerJob(pose,input_tag,nstruct_max),
	alns_(alns),
	template_poses_(templates)
{}

utility::vector1< core::sequence::SequenceAlignment > const &
InnerMultiThreadingJob::alignments() const {
	return alns_;
}

utility::vector1< core::pose::Pose > const &
InnerMultiThreadingJob::template_poses() const {
	return template_poses_;
}


InnerMultiThreadingJob::~InnerMultiThreadingJob() = default;

} // jd2
} // protocols
