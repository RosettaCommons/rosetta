// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/DatabaseJobInputter.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/jd2/DatabaseJobInputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/mysql.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>
#include <cppdb/errors.h>


// C++ headers
#include <string>

static basic::Tracer tr("protocols.jd2.DatabaseJobInputter");

namespace protocols {
namespace jd2 {

using std::string;
using std::endl;
using core::Size;
using core::pose::Pose;
using core::pose::symmetry::is_symmetric;
using core::pose::symmetry::make_asymmetric_pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using protocols::features::ProteinSilentReport;
using utility::file::FileName;
using utility::vector1;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using cppdb::result;



DatabaseJobInputter::DatabaseJobInputter() :
	scfxn_(new ScoreFunction()),
	protein_silent_report_(new ProteinSilentReport())

{
	tr.Debug << "Instantiate DatabaseJobInputter" << endl;
	load_options_from_option_system();
}

DatabaseJobInputter::~DatabaseJobInputter() {}


void
DatabaseJobInputter::load_options_from_option_system(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (option.has(inout::database_filename) &&
		option[inout::database_filename].user()){
		set_database_fname(option[inout::database_filename]);
	}

	// The in::file::tags option was created for the silent file
	// system--but using it makes sense here because, it serves the same
	// purpose: specify which structures to use from the data source.

	if (option.has(in::file::tags) && option[in::file::tags].user()){
		set_tags(option[in::file::tags]);
	}
	if(option[in::file::tags].user() && option[in::select_structures_from_database].user()) {
		utility_exit_with_message("you cannot use -in:file:tags and -in:select_structures_from_database simultaniously");
	}

	if (option[in::select_structures_from_database].user()) {
		set_tags_from_sql(option[in::select_structures_from_database]);
	}
}

void
DatabaseJobInputter::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( inout::database_filename );
	option.add_relevant( in::file::tags );
}

void
DatabaseJobInputter::set_database_fname(
	string const & database_fname
) {
	database_fname_ = database_fname;
}

std::string
DatabaseJobInputter::get_database_fname() const {
	if(database_fname_ == ""){
		utility_exit_with_message(
			"To use the DatabaseJobInputter, please specify the database "
			"where thinput is data is stored, eg. via the -inout:database_filename "
			"<database_fname> option system flag.");
	}
	return database_fname_;
}

/// @brief Get score function
ScoreFunctionOP
DatabaseJobInputter::get_scorefunction(){
	return scfxn_;
}

/// @brief Set score function
void
DatabaseJobInputter::set_scorefunction(ScoreFunctionOP scorefunction ){
	scfxn_ = scorefunction;
}


/// @details The specified tags indicate which structures should be
/// used.  If no tags are specified, then all tags will be used.  Tags
/// are short strings that uniquely identify a structure.  WARNING: As
/// of Jan '11, gives a pose there are several places one my look for
/// such a string identifier:
///    * If the pose has a PDBInfo object: pose.pdb_info()->name()
///    * tag_from_pose(pose) <=> pose.data().get(JOBDIST_OUTPUT_TAG)
///    * JobDistributor::get_instance()->current_job()->input_tag()
///
/// The last two will work, and the first will work when pdb info
/// stored with the structures in the database.

void
DatabaseJobInputter::set_tags(
	vector1< string > const & tags
) {
	tags_ = tags;
}

void DatabaseJobInputter::set_tags_from_sql(utility::vector1<std::string> const & sql)
{
	//first do some basic validation, make sure this is a SELECT command that is selecting the tag or structures.tag
	if(sql[1] != "SELECT" && !(sql[2] == "tag" || sql[2] == "structures.tag"))
	{
		utility_exit_with_message("you must provide an SQL SELECT command that selects the tag or structures.tag column");
	}

	std::string sql_command(utility::join(sql, " "));

	sessionOP db_session(basic::database::get_db_session(database_fname_));

	while(true)
	{
		try
		{
			result res = (*db_session) << sql_command;
			while(res.next()){
				string tag;
				res >> tag;
				tags_.push_back(tag);
			break;
			}

		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
}

void
DatabaseJobInputter::get_tags(
	vector1< string > & tags
) {
	tags = tags_;
}



/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void
DatabaseJobInputter::pose_from_job(
	Pose & pose,
	JobOP job
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	tr.Debug << "DatabaseJobInputter::pose_from_job" << std::endl;
	string tag(job->input_tag());
	pose.clear();

	if ( !job->inner_job()->get_pose() ) {
		tr.Debug << "filling pose from Database (tag = " << tag	<< ")" << endl;
		sessionOP db_session(basic::database::get_db_session(database_fname_));
		while(true)
		{
			try
			{
				protein_silent_report_->load_pose(db_session, tag, pose);
				break;
			}catch(cppdb::cppdb_error & )
			{
				usleep(10);
				continue;
			}
		}

	} else {
		tr.Debug << "filling pose from saved copy (tag = " << tag << ")" << endl;
		pose = *(job->inner_job()->get_pose());
	}

	// TODO: Move to pose.clear()
	if (is_symmetric(pose)) make_asymmetric_pose( pose );
}

/// @details this function determines what jobs exist
void protocols::jd2::DatabaseJobInputter::fill_jobs( Jobs & jobs ){
	tr.Debug << "DatabaseJobInputter::fill_jobs" << std::endl;
	jobs.clear(); //should already be empty anyway

	Size const nstruct(get_nstruct());

	if(!tags_.size()){
		sessionOP db_session(basic::database::get_db_session(database_fname_));

		while(true)
		{
			try
			{
				result res = (*db_session) << "SELECT tag FROM structures;";
				while(res.next()){
					string tag;
					res >> tag;
					tags_.push_back(tag);
				}
				break;
			}catch(cppdb::cppdb_error &)
			{
				usleep(10);
				continue;
			}
		}

	}

	vector1< InnerJobOP > inner_jobs;
	//save list of all inner_jobs first... this allows better sampling
	//of jobs in case of unfinished runs:
	// input1_0001
	// input2_0001
	// ...
	// inputn_0001
	// input1_0002
	// input2_0002
	// ....
	tr.Debug << "reserve memory for InnerJob List " << tags_.size() << endl;
	inner_jobs.reserve( tags_.size() );
	tr.Debug
		<< "fill list with " << tags_.size()
		<< " InnerJob Objects" << endl;

	foreach(string tag, tags_) inner_jobs.push_back(new InnerJob(tag, nstruct));

	tr.Debug
		<< "reserve list for " << inner_jobs.size() * nstruct
		<< " Job Objects" << endl;

	jobs.reserve(inner_jobs.size() * nstruct);

	tr.Debug << "fill job list with... " << endl;
	for ( Size index = 1; index <= nstruct; ++index ) {
		foreach(InnerJobOP ijob, inner_jobs){
			jobs.push_back(new Job(ijob, index));
			tr.Trace
				<< "pushing " << ijob->input_tag() << " nstruct index " << index	<< std::endl;
		}
	}
}

/// @brief Return the type of input source that the
///  DatabaseJobInputter is currently using.
/// @return Always <em>DATABASE</em>.
JobInputterInputSource::Enum DatabaseJobInputter::input_source() const {
	return JobInputterInputSource::DATABASE;
}

} // namespace
} // namespace
