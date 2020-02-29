// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/EmptyPoseJobInputter.cc
/// @brief  Base class EmptyPoseJobInputter
/// @author Danny Farrell

// Unit headers
#include <protocols/jd2/EmptyPoseJobInputter.hh>
#include <protocols/jd2/EmptyPoseJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <core/pose/symmetry/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/sequence/util.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


static basic::Tracer tr( "protocols.jd2.EmptyPoseJobInputter" );

namespace protocols {
namespace jd2 {

EmptyPoseJobInputter::EmptyPoseJobInputter() {
	tr.Debug << "Instantiate EmptyPoseJobInputter" << std::endl;
}

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void EmptyPoseJobInputter::pose_from_job( core::pose::Pose& pose, protocols::jd2::JobOP job) {
	using std::string;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "EmptyPoseJobInputter::pose_from_job" << std::endl;
	if ( !job->inner_job()->get_pose() ) {
		pose.clear();
	} else {
		pose = *(job->inner_job()->get_pose());
		tr.Debug << "filling pose from saved copy " << job->inner_job()->input_tag() << std::endl;
	}
}

/// @details this function determines what jobs exist from -s/-l
void EmptyPoseJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs ){
	tr.Debug << "EmptyPoseJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	core::Size const nstruct( get_nstruct () );

	//note that we are not really using the second and third fields in this implementation
	using basic::options::OptionKeys::jd2::generic_job_name; //This option defaults to 'S' for original behavior
	protocols::jd2::InnerJobOP ijob( new protocols::jd2::InnerJob( basic::options::option[ generic_job_name ].value() , nstruct ) );

	for ( core::Size index = 1; index <= nstruct; ++index ) {
		jobs.push_back( utility::pointer::make_shared< protocols::jd2::Job >( ijob, index ) );
		tr.Trace << "create job index " << index << std::endl;
	}
}

/// @brief Return the type of input source that the EmptyPoseJobInputter is currently using
/// @return Always <em>POSE</em>.
protocols::jd2::JobInputterInputSource::Enum EmptyPoseJobInputter::input_source() const {
	return protocols::jd2::JobInputterInputSource::POSE;
}

//CREATOR SECTION
std::string
EmptyPoseJobInputterCreator::keyname() const
{
	return "EmptyPoseJobInputter";
}

protocols::jd2::JobInputterOP
EmptyPoseJobInputterCreator::create_JobInputter() const {
	return utility::pointer::make_shared< EmptyPoseJobInputter >();
}

}// jd2
}// protocols
