// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/SerializedPoseJobInputter.cc
/// @brief  parses inputs as serialized poses
/// @author Jack Maguire, jackmaguire1444@gmail.com

///Unit headers
#include <protocols/jd2/SerializedPoseJobInputter.hh>
#include <protocols/jd2/SerializedPoseJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers

#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/version.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

///C++ headers
#include <string>

#include <core/import_pose/import_pose.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>
#endif

static basic::Tracer TR( "protocols.jd2.SerializedPoseJobInputter" );

namespace protocols {
namespace jd2 {

SerializedPoseJobInputter::SerializedPoseJobInputter(){
#ifndef SERIALIZATION
	utility_exit_with_message( "please build with extras=serialization if you want to use the SerializedPoseJobInputter" );
#endif

	TR << "Instantiate SerializedPoseJobInputter" << std::endl;
}

SerializedPoseJobInputter::~SerializedPoseJobInputter()= default;

/// @details This function will first see if the pose already exists in the Job.  If not,
/// it will read it into the pose reference, and hand a COP cloned from that pose to the
/// Job. If the pose pre-exists it just copies the COP's pose into it.
void SerializedPoseJobInputter::pose_from_job( core::pose::Pose & pose, JobOP job){

	TR << "SerializedPoseJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
		pose_from_file( job->input_tag(), pose );
		load_pose_into_job(pose, job);
	} else {
		pose = * ( job->inner_job()->get_pose() );
	}
}

#ifdef SERIALIZATION
void SerializedPoseJobInputter::pose_from_file( std::string const & filename, core::pose::Pose & pose ) const {
	std::string const contents = utility::file_contents( filename );
	std::istringstream iss( contents );

	//First, check for correct commit
	std::string const current_commit = utility::Version::commit();
	std::string saved_commit;
	std::getline( iss, saved_commit, '\n' );
	if ( saved_commit != current_commit ) {
		if ( basic::options::option[ basic::options::OptionKeys::in::file::srlz_override ].value() ) {
			TR << "Pose " << filename << " was saved by a different version of Rosetta but in:file:srlz_override"
				" was set to true so we will attempt to deserialize it anyways." << std::endl;
			TR << "Original Rosetta Sha1: " << saved_commit << std::endl;
		} else {
			utility_exit_with_message( "Error! Can Not Load File: " + filename + "!\n"
				"File was created with a previous version of Rosetta.\n"
				"Commit '" + saved_commit + "' is required to load this file." );
		}
	}

	cereal::BinaryInputArchive arc( iss );
	arc( pose );
}
#else
void SerializedPoseJobInputter::pose_from_file( std::string const &, core::pose::Pose & ) const {
}
#endif

/// @details this function determines what jobs exist from -s/-l
void SerializedPoseJobInputter::fill_jobs( JobsContainer & jobs ){
	TR << "SerializedPoseJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	utility::vector1< std::string > const inputs( basic::options::start_files() );
	core::Size const nstruct( get_nstruct() );

	for ( core::Size i(1); i <= inputs.size(); ++i ) {

		//protocols::jobdist::BasicJob = InnerJob
		//note that we are not really using the second and third fields in this implementation
		InnerJobOP ijob( new InnerJob( inputs[i], nstruct ) );

		for ( core::Size index(1); index <= nstruct; ++index ) {
			jobs.push_back( JobOP( new Job( ijob, index ) ) );
			TR.Debug << "pushing " << inputs[i] << " nstruct index " << index << std::endl;
		}//loop over nstruct
		TR << "pushed " << inputs[i] << " nstruct "
			<< ( ( int(nstruct)==1 ) ? "index " : "indices 1 - ")
			<< nstruct << std::endl;
	} //loop over inputs
} //fill_jobs


//CREATOR SECTION
std::string
SerializedPoseJobInputterCreator::keyname() const
{
	return "SerializedPoseJobInputter";
}

JobInputterOP
SerializedPoseJobInputterCreator::create_JobInputter() const {
	return utility::pointer::make_shared< SerializedPoseJobInputter >();
}

}//jd2
}//protocols
