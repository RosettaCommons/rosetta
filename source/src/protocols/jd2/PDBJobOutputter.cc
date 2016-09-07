// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/PDBJobOutputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - PDBJobOutputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/PDBJobOutputter.hh>
#include <protocols/jd2/PDBJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/util.hh>
#include <core/io/pdb/pdb_writer.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
//#include <utility/file/FileName.hh>
#include <core/types.hh>
#include <basic/options/option.hh>

///Basic headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

///C++ headers
#include <string>
#include <map>
#include <sstream>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/io/pdb/build_pose_as_is.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.PDBJobOutputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::PDBJobOutputter::PDBJobOutputter()
: parent()
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	TR.Debug << "Using PDBJobOutputter for JobDistributor" << std::endl;

	set_extension(".pdb");
	if ( option[ out::pdb_gz ] ) {
		set_extension(".pdb.gz");
	}

	if ( option[ out::path::pdb ].user() ) {
		set_path( option[ out::path::pdb ]().path());
	}
}

protocols::jd2::PDBJobOutputter::~PDBJobOutputter()= default;

/// @details private function (just prevents code duplication) to fill ozstream
void protocols::jd2::PDBJobOutputter::dump_pose(
	JobCOP job,
	core::pose::Pose const & pose,
	utility::io::ozstream & out,
	std::string const & filename /* filename is an optional label in the score data table */
)
{
	core::io::pdb::dump_pdb(
		pose,
		extract_data_from_Job(job),
		out,
		filename /* filename is an optional label in the score data table */
	);

}

/////////////////////////////////state of output functions/////////////////////////////////
/*
bool protocols::jd2::PDBJobOutputter::job_has_completed( JobCOP job ){
bool complete = utility::file::file_exists( path_ + extended_name(job) );
if ( TR.Debug.visible() && complete ) {
TR.Debug << "Skipping job " << output_name(job) << " because the output file already exists on disk." << std::endl;
}
return complete;
}

std::string
protocols::jd2::PDBJobOutputter::extended_name( JobCOP job, std::string const & suffix )
{
return output_name(job) + std::string(suffix.empty() ? "" : "_") + suffix + extension_;
}

////////////////////////////////////////score-related functions///////////////////////////////////


/// @brief this function extracts the pose's scores and outputs them as a string to be packaged in an output structure.
/// @details Refactored in the 2016 Chemical XRW (eXtreme Rosetta Workshop) by Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in] job Const-access owning pointer to the job from which the data will be extracted.
/// @return A string in which the data will be stored, that can later be passed to whatever container wants it.
std::string
protocols::jd2::PDBJobOutputter::extract_data_from_Job(
JobCOP job
) {
//TR << "protocols::jd2::PDBJobOutputter::extract_data_from_Job" << std::endl;

std::stringstream out;

for ( Job::Strings::const_iterator it(job->output_strings_begin()), end(job->output_strings_end());
it != end;
++it ) {
out << *it << std::endl;
//TR << *it << std::endl;
}

for ( Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin()), end(job->output_string_string_pairs_end());
it != end;
++it ) {
out << it->first << " " << it->second << std::endl;
//TR << it->first << " " << it->second << std::endl;
}

for ( Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end());
it != end;
++it ) {
out << it->first << " " << it->second << std::endl;
//TR << it->first << " " << it->second << std::endl;
}

return out.str();
}
*/

//CREATOR SECTION
std::string
PDBJobOutputterCreator::keyname() const
{
	return "PDBJobOutputter";
}

protocols::jd2::JobOutputterOP
PDBJobOutputterCreator::create_JobOutputter() const {
	return protocols::jd2::JobOutputterOP( new PDBJobOutputter );
}

}//jd2
}//protocols
