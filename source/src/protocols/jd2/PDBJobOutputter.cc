// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PDBJobOutputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - PDBJobOutputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/PDBJobOutputter.hh>
#include <protocols/jd2/PDBJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>

///Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>

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

#include <core/io/pdb/file_data.hh>
#include <utility/vector1.hh>




static basic::Tracer TR("protocols.jd2.PDBJobOutputter");

namespace protocols {
namespace jd2 {

protocols::jd2::PDBJobOutputter::PDBJobOutputter()
	: parent(), extension_(".pdb")
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	TR.Debug << "Using PDBJobOutputter (base class) for JobDistributor" << std::endl;

	if ( option[ out::pdb_gz ] ) {
		extension_ = ".pdb.gz";
	}

	if ( option[ out::path::pdb ].user() ) {
		path_ = option[ out::path::pdb ]().path();
	} else if (option[ out::path::all ].user() ) {
		path_ = option[ out::path::all ]().path();
	} else {
		path_ = "";
	}
}

protocols::jd2::PDBJobOutputter::~PDBJobOutputter(){}

//////////////////////////////creating output functions/////////////////////////////////////////

void protocols::jd2::PDBJobOutputter::final_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag
)
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	TR.Debug << "PDBJobOutputter::final_pose" << std::endl;

	call_output_observers( pose, job );
	std::string const file(path_ + extended_name(job, tag));
	utility::io::ozstream out( file );
	if ( !out.good() ) utility_exit_with_message( "Unable to open file: " + file + "\n" );
	dump_pose(job, pose, out);
	out.close();
	scorefile(job, pose, "", (tag.empty() ? "" : std::string("_") + tag));
}

void protocols::jd2::PDBJobOutputter::other_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int /*copy_count*/, /*default -1 */
	bool /* score_only*/ /*default false*/
){
	TR.Debug << "PDBJobOutputter::other_pose" << std::endl;
	runtime_assert( !tag.empty() ); //else you'll overwrite your pdb when the job finishes

	call_output_observers( pose, job );
	std::string const file(path_ + tag + "_" + extended_name(job));
	utility::io::ozstream out( file );
	if ( !out.good() ) utility_exit_with_message( "Unable to open file: " + file + "\n" );
	dump_pose(job, pose, out);
	out.close();

	//these are separate options because leaving the default on other_pose_scorefile is totally valid, but you can't both specify it on the command line and leave it blank
	//THIS FUNCTIONALITY IS GOING TO BE DEPRECATED SOON
	if( basic::options::option[ basic::options::OptionKeys::run::other_pose_to_scorefile ].value() ){
		scorefile(job, pose, tag, "", basic::options::option[ basic::options::OptionKeys::run::other_pose_scorefile ].value());
	}
}

///@details private function (just prevents code duplication) to fill ozstream
void protocols::jd2::PDBJobOutputter::dump_pose(
	JobCOP job,
	core::pose::Pose const & pose,
	utility::io::ozstream & out
)
{
	core::io::pdb::FileData::dump_pdb( pose, out );
	extract_scores(pose, out);
	extract_extra_scores( pose, out );
	extract_data_from_Job( job, out );
}

/////////////////////////////////state of output functions/////////////////////////////////
bool protocols::jd2::PDBJobOutputter::job_has_completed( JobCOP job ){
	return utility::file::file_exists( path_ + extended_name(job) );
}

std::string protocols::jd2::PDBJobOutputter::output_name( JobCOP job ){
	return affixed_numbered_name( job );
}

std::string protocols::jd2::PDBJobOutputter::extended_name( JobCOP job, std::string const suffix ){
	return output_name(job) + std::string(suffix.empty() ? "" : "_") + suffix + extension_;
}

////////////////////////////////////////score-related functions///////////////////////////////////

///@brief this function extracts the pose's scores for printing
/// @details In order for this to work as expected, the Pose's cached energies
/// must match up with the (current) conformation.
/// A good time to do this is at the end of your protocol's apply() method:
///   scorefxn( pose );
void protocols::jd2::PDBJobOutputter::extract_scores(
	core::pose::Pose const & pose,
	utility::io::ozstream & out
)
{
	core::io::pdb::extract_scores( pose, out );
}

//THIS FUNCTION WILL MOVE HIGHER IN THE HIERARCHY AT SOME POINT
///@brief this function extracts the pose's scores and outputs them into the PDB
void protocols::jd2::PDBJobOutputter::extract_data_from_Job( JobCOP job, utility::io::ozstream & out ){
	//TR << "protocols::jd2::PDBJobOutputter::extract_data_from_Job" << std::endl;

	for( Job::Strings::const_iterator it(job->output_strings_begin()), end(job->output_strings_end());
			 it != end;
			 ++it) {
		out << *it << std::endl;
		//TR << *it << std::endl;
	}

	for( Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin()), end(job->output_string_string_pairs_end());
			 it != end;
			 ++it) {
		out << it->first << " " << it->second << std::endl;
		//TR << it->first << " " << it->second << std::endl;
	}

	for( Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end());
			 it != end;
			 ++it) {
		out << it->first << " " << it->second << std::endl;
		//TR << it->first << " " << it->second << std::endl;
	}
}


///@brief This function extracts the pose's extra data/scores and outputs them into the PDB
///@details YOU are responsible for putting things into the pose's DataCache using core::pose::util::setPoseExtraScore().
/// Both string_real and string_real pairs can be stored using setPoseExtraScore().
void protocols::jd2::PDBJobOutputter::extract_extra_scores(
	core::pose::Pose const & pose,
	utility::io::ozstream & out
)
{
	// ARBITRARY_STRING_DATA
	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		basic::datacache::CacheableStringMapCOP data
			= dynamic_cast< basic::datacache::CacheableStringMap const * >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		assert( data.get() != NULL );

		for(std::map< std::string, std::string >::const_iterator it( data->map().begin() ), end( data->map().end() );
		  it != end;
		  ++it) {
			//TR << it->first << " " << it->second << std::endl;
			out << it->first << " " << it->second << std::endl;
		}
	}

	// ARBITRARY_FLOAT_DATA
	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		basic::datacache::CacheableStringFloatMapCOP data
			= dynamic_cast< basic::datacache::CacheableStringFloatMap const * >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		assert( data.get() != NULL );

		for(std::map< std::string, float >::const_iterator it( data->map().begin() ), end( data->map().end() );
		  it != end;
		  ++it) {
			//TR << it->first << " " << it->second << std::endl;
			out << it->first << " " << it->second << std::endl;
		}
	}
}

//CREATOR SECTION
std::string
PDBJobOutputterCreator::keyname() const
{
        return "PDBJobOutputter";
}

protocols::jd2::JobOutputterOP
PDBJobOutputterCreator::create_JobOutputter() const {
        return new PDBJobOutputter;
}

}//jd2
}//protocols
