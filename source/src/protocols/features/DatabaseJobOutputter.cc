// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

// Project Headers
#include <protocols/features/ProteinSilentReport.hh>

// Platform Headers
#include <core/pose/util.hh>

///Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>

//Auto Headers
#include <core/scoring/ScoreFunction.fwd.hh>
static thread_local basic::Tracer tr( "protocols.features.DatabaseJobOutputter" );

namespace protocols {
namespace features {

using std::string;
using core::Size;
using core::pose::Pose;
using protocols::jd2::JobCOP;
using protocols::jd2::JobOP;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using basic::database::get_db_session;
using cppdb::result;


DatabaseJobOutputter::DatabaseJobOutputter() :
	protocols::jd2::FileJobOutputter(),
	protein_silent_report_(protocols::features::ProteinSilentReportOP( new ProteinSilentReport() )),
	database_name_(),
	database_pq_schema_()
{
	load_options_from_option_system();
	sessionOP db_session(
		get_db_session(path_ + database_name_, database_pq_schema_));
	protein_silent_report_->initialize(db_session);

}

DatabaseJobOutputter::~DatabaseJobOutputter() {
	//DO NOT PUT THINGS HERE - it is not guaranteed to get called - use flush below instead.
}

void
DatabaseJobOutputter::load_options_from_option_system(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option.has(inout::dbms::database_name) &&
			option[inout::dbms::database_name].user() ) {
		set_database_name(option[inout::dbms::database_name]);
	}

	if ( option.has(inout::dbms::pq_schema) &&
			option[inout::dbms::pq_schema].user() ) {
		set_database_pq_schema(option[inout::dbms::pq_schema]);
	}

	if ( option[ out::path::db ].user() ) {
		path_ = option[ out::path::db ]().path();
	} else if ( option[ inout::dbms::path ].user() ) {
		path_ = option[ inout::dbms::path ]().path();
	} else if ( option[ out::path::all ].user() ) {
		path_ = option[ out::path::all ]().path();
	} else {
		path_ = "";
	}

}

void
DatabaseJobOutputter::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( inout::dbms::database_name );
	option.add_relevant( inout::dbms::pq_schema );
	option.add_relevant( inout::dbms::host );
	option.add_relevant( inout::dbms::user );
	option.add_relevant( inout::dbms::password );
	option.add_relevant( inout::dbms::port );
	option.add_relevant( inout::dbms::separate_db_per_mpi_process );
	option.add_relevant( out::resume_batch);
	option.add_relevant( out::path::db);
	option.add_relevant( out::path::all);
	option.add_relevant( inout::dbms::path);

}

void
DatabaseJobOutputter::set_database_pq_schema(
	string const & database_pq_schema
) {
	database_name_ = database_pq_schema;
}

string
DatabaseJobOutputter::get_database_pq_schema() const {
	return database_pq_schema_;
}

void
DatabaseJobOutputter::set_database_name(
	string const & database_name
) {
	database_name_ = database_name;
}

std::string
DatabaseJobOutputter::get_database_name() const {
	if ( database_name_ == "" ) {
		utility_exit_with_message(
			"To use the DatabaseJobInputter, please specify the database "
			"where the input is data is stored, eg. via the -inout:dbms:database_name "
			"<database_name> option system flag.");
	}
	return database_name_;
}


void DatabaseJobOutputter::flush() {
}

void DatabaseJobOutputter::final_pose(
	JobOP job,
	Pose const & pose,
	std::string const & /*tag*/
) {

	call_output_observers( pose, job );

	// If this is bottle neck, consider hanging on to the db_session
	// rather than recreating it each time.

	sessionOP db_session(get_db_session(path_ + database_name_, database_pq_schema_));
	protein_silent_report_->apply(pose, db_session, output_name(job));
}

/// @brief this function is intended for saving mid-protocol poses; for example
/// the final centroid structure in a combined centroid/fullatom protocol.
void DatabaseJobOutputter::other_pose(
	JobOP job,
	Pose const & pose,
	string const & tag,
	int , /*default -1 */
	bool /*default false*/
) {

	call_output_observers( pose, job );

	sessionOP db_session(get_db_session(path_ + database_name_, database_pq_schema_));
	protein_silent_report_->apply(pose, db_session, tag);

}

/////////////////////////////////state of output functions/////////////////////////////////
bool DatabaseJobOutputter::job_has_completed(
	JobCOP job
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// did we complete the job later ?
	if ( job->completed() ) {
		return true;
	}
	sessionOP db_session(get_db_session(path_ + database_name_, database_pq_schema_));

	//It is possible for the mpi distributor to call this function
	//before the database has even been initialized

	protein_silent_report_->initialize(db_session);

	result res;
	if ( option[out::resume_batch].user() ) {
		utility::vector1<core::Size> batch_ids(option[out::resume_batch].value());
		core::Size placeholder_count = batch_ids.size();
		std::string placeholder_block= "(?";
		for ( Size j = 1; j < placeholder_count; ++j ) {
			placeholder_block += ",?";
		}
		placeholder_block += ")";

		std::string job_completion_string = "SELECT count(*) FROM sampled_structures WHERE tag=? AND batch_id IN " +placeholder_block+";";
		cppdb::statement job_completion_statement(basic::database::safely_prepare_statement(job_completion_string,db_session));
		job_completion_statement.bind(1,output_name(job));
		for ( Size i = 1; i <= batch_ids.size(); ++i ) {
			core::Size column_index =i+1;
			job_completion_statement.bind(column_index,batch_ids[i]);
		}
		res = basic::database::safely_read_from_database(job_completion_statement);

	} else {
		std::string job_completion_string = "SELECT count(*) FROM sampled_structures WHERE tag=? and batch_id = ?;";
		cppdb::statement job_completion_statement(basic::database::safely_prepare_statement(job_completion_string,db_session));
		job_completion_statement.bind(1,output_name(job));
		job_completion_statement.bind(2,protein_silent_report_->get_batch_id());
		res = basic::database::safely_read_from_database(job_completion_statement);
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
	return protocols::jd2::JobOutputterOP( new DatabaseJobOutputter );
}

} // namespace features
} // namespace protocols
