// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/wwPDBJobOutputter.cc
/// @brief  wwPDBJobOutputter class, interstitial between FileJobOutputter and PDBJobOutputter or mmCIFJobOutputter to share functionality.
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/wwPDBJobOutputter.hh>
//#include <protocols/jd2/wwPDBJobOutputterCreator.hh> //abstract class b/c dump_pose=0, no need for creator
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


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.wwPDBJobOutputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::wwPDBJobOutputter::wwPDBJobOutputter()
: parent(), path_(""), extension_("")
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	//Old logic was out:path:PDB > out:path:all.  Since base class ctors run before child class, checking all here A) works for PDB or mmCIF shared, and B) maintains old logic.
	if ( option[ out::path::all ].user() ) {
		path_ = option[ out::path::all ]().path();
	}
}

protocols::jd2::wwPDBJobOutputter::~wwPDBJobOutputter(){}

//////////////////////////////creating output functions/////////////////////////////////////////

/// @details originally written for PDBJobOutputter, but works as-is for mmCIF as well.  Moved to interstitial wwPDBJO.
void protocols::jd2::wwPDBJobOutputter::final_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag
)
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	TR.Debug << "wwPDBJobOutputter::final_pose" << std::endl;

	call_output_observers( pose, job );
	std::string const base_file ( extended_name(job, tag) );
	std::string const file(path_ + base_file );
	utility::io::ozstream out( file );
	if ( !out.good() ) utility_exit_with_message( "Unable to open file: " + file + "\n" );
	dump_pose(job, pose, out, base_file );
	out.close();
	scorefile(job, pose, "", (tag.empty() ? "" : std::string("_") + tag));
}

/// @details originally written for PDBJobOutputter, but works as-is for mmCIF as well.  Moved to interstitial wwPDBJO.
void protocols::jd2::wwPDBJobOutputter::other_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int /*copy_count*/, /*default -1 */
	bool /* score_only*/ /*default false*/
){
	TR.Debug << "wwPDBJobOutputter::other_pose" << std::endl;
	runtime_assert( !tag.empty() ); //else you'll overwrite your pdb when the job finishes

	call_output_observers( pose, job );
	std::string const base_file( tag + "_" + extended_name(job) );
	std::string const file(path_ + base_file );
	utility::io::ozstream out( file );
	if ( !out.good() ) utility_exit_with_message( "Unable to open file: " + file + "\n" );
	dump_pose(job, pose, out, base_file );
	out.close();

	//these are separate options because leaving the default on other_pose_scorefile is totally valid, but you can't both specify it on the command line and leave it blank
	//THIS FUNCTIONALITY IS GOING TO BE DEPRECATED SOON. WHY???
	if ( basic::options::option[ basic::options::OptionKeys::run::other_pose_to_scorefile ].value() ) {
		scorefile(job, pose, tag, "", basic::options::option[ basic::options::OptionKeys::run::other_pose_scorefile ].value());
	}
}

/////////////////////////////////state of output functions/////////////////////////////////
/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb/cif with the job's name already in existence.
bool protocols::jd2::wwPDBJobOutputter::job_has_completed( JobCOP job ){
	return utility::file::file_exists( path_ + extended_name(job) );
}

/// @brief determines the unique output identifier for a job
std::string protocols::jd2::wwPDBJobOutputter::output_name( JobCOP job ){
	return affixed_numbered_name( job );
}

std::string
protocols::jd2::wwPDBJobOutputter::extended_name( JobCOP job, std::string const & suffix )
{
	return output_name(job) + std::string(suffix.empty() ? "" : "_") + suffix + extension_;
}

////////////////////////////////////////score-related functions///////////////////////////////////
/// @brief this function extracts the pose's scores and outputs them as a string to be packaged in an output structure.
/// @details Refactored in the 2016 Chemical XRW (eXtreme Rosetta Workshop) by Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in] job Const-access owning pointer to the job from which the data will be extracted.
/// @return A string in which the data will be stored, that can later be passed to whatever container wants it.
std::string
protocols::jd2::wwPDBJobOutputter::extract_data_from_Job(
	JobCOP job
) {
	//TR << "protocols::jd2::wwPDBJobOutputter::extract_data_from_Job" << std::endl;

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


///////////////////protected: child classes PDBJO and mmCIFJO need to write these data////////////////////////////////

//Note this function is protected (as opposed to private or public)
/// @brief setter for output file paths, in case child class needs to override -out:path:all with -out:path:[PDB/mmCIF]
void wwPDBJobOutputter::set_path(std::string const & path) {path_ = path;}

//Note this function is protected (as opposed to private or public)
/// @brief getter for output file path
std::string const & wwPDBJobOutputter::get_path() {return path_;}

//Note this function is protected (as opposed to private or public)
/// @brief setter for output file extensions
void wwPDBJobOutputter::set_extension(std::string const & extension) {extension_ = extension;}

//Note this function is protected (as opposed to private or public)
/// @brief getter for output file extension
std::string const & wwPDBJobOutputter::get_extension() {return extension_;}


}//jd2
}//protocols
