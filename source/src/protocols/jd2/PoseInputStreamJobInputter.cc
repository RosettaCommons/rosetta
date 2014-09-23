// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PoseInputStreamJobInputter.cc
/// @brief
/// @author James Thompson

#include <protocols/jd2/PoseInputStreamJobInputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/PoseInputStreamJobInputterCreator.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <ObjexxFCL/string.functions.hh>

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


static thread_local basic::Tracer tr( "protocols.jd2.PoseInputStreamJobInputter" );

namespace protocols {
namespace jd2 {

PoseInputStreamJobInputter::PoseInputStreamJobInputter() {
	tr.Debug << "Instantiate PoseInputStreamJobInputter" << std::endl;
}

PoseInputStreamJobInputter::~PoseInputStreamJobInputter() {}

void PoseInputStreamJobInputter::pose_from_job(
	core::pose::Pose & pose, JobOP /*job*/
) {
	tr.Debug << "PoseInputStreamJobInputter::pose_from_job" << std::endl;
	pose.clear();
	if ( input_.has_another_pose() ) {
		input_.fill_pose( pose, *rsd_set_ );
	} else {
		utility_exit_with_message( "Error: no more poses!" );
	}
}

/// @details this function determines what jobs exist from -in::file::silent and
/// -in::file::tags
void PoseInputStreamJobInputter::fill_jobs( Jobs & jobs ) {
	tr.Debug << "PoseInputStreamJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	using core::Size;
	using std::string;
	using utility::vector1;
	using core::pose::Pose;
	using core::chemical::rsd_set_from_cmd_line;
	using namespace core::import_pose::pose_stream;

	input_ = streams_from_cmd_line();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Size const nstruct( option[ out::nstruct ]() );
	Size const shuffle_nstruct( option[ out::shuffle_nstruct ]() );
	vector1< JobOP > inner_jobs;
	rsd_set_ = core::chemical::ResidueTypeSetCOP( rsd_set_from_cmd_line() );

	for ( core::Size index = 1; index <= shuffle_nstruct; ++index ) {
		InnerJobOP ijob( new InnerJob( "job_" + ObjexxFCL::string_of(index), nstruct ) );
		jobs.push_back( JobOP( new Job( ijob, index ) ) );
		tr.Trace << "pushing " << ijob->input_tag() << " nstruct index " << index << std::endl;
	}

} // fill_jobs

JobInputterInputSource::Enum PoseInputStreamJobInputter::input_source() const {
	return JobInputterInputSource::UNKNOWN;
}

//CREATOR SECTION
std::string
PoseInputStreamJobInputterCreator::keyname() const
{
        return "PoseInputStreamJobInputter";
}

protocols::jd2::JobInputterOP
PoseInputStreamJobInputterCreator::create_JobInputter() const {
        return protocols::jd2::JobInputterOP( new PoseInputStreamJobInputter );
}

} // jd2
} // protocols
