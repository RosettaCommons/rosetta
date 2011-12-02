// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/DatabaseJobOutputter.cc
/// @brief  Job Outputter to a database
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <protocols/features/DatabaseJobOutputter.hh>
#include <protocols/features/DatabaseJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>
// AUTO-REMOVED #include <protocols/jd2/util.hh>

// Project Headers
#include <protocols/features/ProteinSilentReport.hh>

// Platform Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>

///Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/database/sql_utils.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>
// AUTO-REMOVED #include <cppdb/errors.h>

// C++ Headers
#include <string>

//Auto Headers
#include <core/scoring/ScoreFunction.fwd.hh>
static basic::Tracer tr("protocols.features.DatabaseJobOutputter");

namespace protocols {
namespace features {

using std::string;
using core::Size;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using protocols::features::ProteinSilentReport;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using cppdb::result;


DatabaseJobOutputter::DatabaseJobOutputter() :
	protein_silent_report_(new ProteinSilentReport())
{
	load_options_from_option_system();

	//sessionOP db_session(basic::database::get_db_session(database_fname_));

	//protein_silent_report_->write_schema_to_db(db_session);
}

DatabaseJobOutputter::~DatabaseJobOutputter() {
	//DO NOT PUT THINGS HERE - it is not guaranteed to get called - use flush below instead.
}

void
DatabaseJobOutputter::load_options_from_option_system(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (option.has(inout::database_filename) &&
		option[inout::database_filename].user()){
		set_database_fname(option[inout::database_filename]);
	}

}

void
DatabaseJobOutputter::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( inout::database_filename );
}

void
DatabaseJobOutputter::set_database_fname(
	string const & database_fname
) {
	database_fname_ = database_fname;
}

std::string
DatabaseJobOutputter::get_database_fname() const {
	if(database_fname_ == ""){
		utility_exit_with_message(
			"To use the DatabaseJobInputter, please specify the database "
			"where the input is data is stored, eg. via the -inout:database_filename "
			"<database_fname> option system flag.");
	}
	return database_fname_;
}

void DatabaseJobOutputter::flush() {
}

void DatabaseJobOutputter::final_pose(
	protocols::jd2::JobCOP job, core::pose::Pose const & pose
) {
	sessionOP db_session(basic::database::get_db_session(database_fname_));
	protein_silent_report_->apply(pose, db_session, output_name(job));
}

/// @brief this function is intended for saving mid-protocol poses; for example
/// the final centroid structure in a combined centroid/fullatom protocol.
void DatabaseJobOutputter::other_pose(
	protocols::jd2::JobCOP,
	core::pose::Pose const & pose,
	std::string const & tag,
	int copy_count, /*default -1 */
	bool score_only /*default false*/
) {

	sessionOP db_session(basic::database::get_db_session(database_fname_));
	protein_silent_report_->apply(pose, db_session, tag);

}

/////////////////////////////////state of output functions/////////////////////////////////
bool DatabaseJobOutputter::job_has_completed( protocols::jd2::JobCOP job ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// did we complete the job later ?
	if ( job->completed() ) {
		return true;
	}
	sessionOP db_session(basic::database::get_db_session(database_fname_));

	//It is possible for the mpi distributor to call this function
	//before the database has even been initialized
	//if this is the case, return false
	if(!protein_silent_report_->is_initialized())
	{
		return false;
	}


	result res;
	while(true)
	{
		try
		{
			res = (*db_session) <<
				"SELECT count(*) FROM structures WHERE tag=?;" << output_name(job);
			break;
		}catch(cppdb::cppdb_error &)
		{
			#ifndef WIN32
				usleep(10);
			#endif
			continue;
		}
	}

	res.next();
	Size already_written;
	res >> already_written;

	return already_written;


}

/// @details
/// Database tags should preserve the FULL NAME such that we don't end up with
/// duplicate tags. This will cause problems on BOINC if changed.
std::string DatabaseJobOutputter::output_name( protocols::jd2::JobCOP job ) {
	return affixed_numbered_name( job );
}

//CREATOR SECTION
std::string
DatabaseJobOutputterCreator::keyname() const
{
	return "DatabaseJobOutputter";
}

protocols::jd2::JobOutputterOP
DatabaseJobOutputterCreator::create_JobOutputter() const {
	return new DatabaseJobOutputter;
}

} // namespace features
} // namespace protocols
