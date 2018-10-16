// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/SerializedPoseJobOutputter.cc
/// @brief  outputs serialized Poses that should only be read in by applications from the same git checkout on the same computer
/// @author Jack Maguire, jackmaguire1444@gmail.com

///Unit headers
#include <protocols/jd2/SerializedPoseJobOutputter.hh>
#include <protocols/jd2/SerializedPoseJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/io/util.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/memory.hh>
#include <core/types.hh>
#include <basic/options/option.hh>

///C++ headers
#include <string>
#include <map>
#include <sstream>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/version.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>
#endif

static basic::Tracer TR( "protocols.jd2.SerializedPoseJobOutputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::SerializedPoseJobOutputter::SerializedPoseJobOutputter() :
	protocols::jd2::wwPDBJobOutputter()
{
#ifndef SERIALIZATION
	utility_exit_with_message( "please build with extras=serialization if you want to use the SerializedPoseJobOutputter" );
#endif

	using namespace basic::options::OptionKeys;
	using basic::options::option;

	TR.Debug << "Using SerializedPoseJobOutputter for JobDistributor" << std::endl;

	set_extension( ".srlz" );
}

protocols::jd2::SerializedPoseJobOutputter::~SerializedPoseJobOutputter() = default;

/// @details private function (just prevents code duplication) to fill ozstream
#ifdef SERIALIZATION
void protocols::jd2::SerializedPoseJobOutputter::dump_pose(
	JobCOP,
	core::pose::Pose const & pose,
	utility::io::ozstream & out,
	std::string const & /* filename is an optional label in the score data table */
) {
	//serialize pose, stolen from jd3::MPIWorkPoolJobDistributor::serialize_larval_job
	std::ostringstream oss;
	oss << utility::Version::commit() << '\n';
	cereal::BinaryOutputArchive arc( oss );
	arc( pose );
	out << oss.str();
}
#else
void protocols::jd2::SerializedPoseJobOutputter::dump_pose(
	JobCOP,
	core::pose::Pose const &,
	utility::io::ozstream &,
	std::string const & /* filename is an optional label in the score data table */
) {}
#endif

//CREATOR SECTION
std::string
SerializedPoseJobOutputterCreator::keyname() const
{
	return "SerializedPoseJobOutputter";
}

protocols::jd2::JobOutputterOP
SerializedPoseJobOutputterCreator::create_JobOutputter() const {
	return utility::pointer::make_shared< SerializedPoseJobOutputter >();
}

}//jd2
}//protocols
