// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/DatabaseJobInputter.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/DatabaseJobInputter.hh>
#include <protocols/features/DatabaseJobInputterCreator.hh>
#include <protocols/features/ProteinSilentReport.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>


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
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>


// C++ headers
#include <string>

static basic::Tracer tr("protocols.features.DatabaseJobInputter");

namespace protocols {
namespace features {

using std::string;
using std::endl;
using core::Size;
using core::pose::initialize_disulfide_bonds;
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

	//TODO allow list of struct ids (it's not immediately clear how to do this with SQLite due to the fact that there is no simple "unhex" functionality
//	if (option.has(in::file::tags) && option[in::file::tags].user()){
//		set_tags(option[in::file::tags]);
//	}
	if(option[in::file::tags].user() && option[in::select_structures_from_database].user()) {
		utility_exit_with_message("you cannot use -in:file:tags and -in:select_structures_from_database simultaniously");
	}

	if (option[in::select_structures_from_database].user()) {
		set_struct_ids_from_sql(option[in::select_structures_from_database]);
	}

	//TODO do we want this still?
//	input_protocol_id_ = option[in::database_protocol];

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


/// @details The specified struct_ids indicate which structures should be
/// used.  If no ids are specified, then all will be used.  Unless a tag column
/// is specified in the SQL statement, the job name (and
/// consequently, the file output name) will be an ASCII hexadecimal representation
/// of the struct_id (a boost UUID). If a tag column is given, then the file name will
/// be the tag associated with the given row.
void DatabaseJobInputter::set_struct_ids_from_sql(utility::vector1<std::string> const & sql)
{
	std::string sql_command(utility::join(sql, " "));
	basic::database::check_statement_sanity(sql_command);

	sessionOP db_session(basic::database::get_db_session(database_fname_));

	result res;
	while(true)
	{
		try
		{
			res = (*db_session) << sql_command;
			break;
		}catch(cppdb::cppdb_error &)
		{
			#ifndef WIN32
				usleep(10);
			#endif
			continue;
		}
	}
	
	bool res_nums_specified = false;
	if(res.find_column("resnum") != -1){res_nums_specified=true;}
	
	bool tags_specified = false;
	if(res.find_column("tag") != -1){tags_specified=true;}
	
	if(res.find_column("struct_id") != -1){	
		while(res.next()){
			boost::uuids::uuid struct_id;
			res.fetch("struct_id", struct_id);
			
			std::string tag;
			if(tags_specified){
				res.fetch("tag", tag);
				if(tag_structures_.count(tag) > 0 && tag_structures_[tag] != struct_id){
					utility_exit_with_message("You have specified non-unque input tags which can cause ambigous output. Please make input tags unique");
				}				
			}
			else{
				tag = to_string(struct_id);
			}
			tag_structures_[tag] = struct_id;
			
			if(res_nums_specified){
				core::Size resnum;
				res.fetch("resnum", resnum);
				tag_residues_[tag].insert(resnum);
			}
		}
		if(!tag_structures_.size()){
			utility_exit_with_message("The provided SQL query did not produce any struct_ids");
		}
	}
	else{
		utility_exit_with_message("Must provide an SQL SELECT command that selects the struct_id column from the structures table");
	}
}

//void
//DatabaseJobInputter::get_tags(
//	vector1< string > & tags
//) {
//	tags = tags_;
//}



/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void
DatabaseJobInputter::pose_from_job(
	Pose & pose,
	protocols::jd2::JobOP job
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	tr.Debug << "DatabaseJobInputter::pose_from_job" << std::endl;
	string tag(job->input_tag());
	pose.clear();

	if ( !job->inner_job()->get_pose() ) {
		tr.Debug << "filling pose from Database (input tag = " << tag << ")" << endl;
		sessionOP db_session(basic::database::get_db_session(database_fname_));

		boost::uuids::uuid struct_id = struct_id=tag_structures_[tag];

		if(!tag_residues_.size()){
			protein_silent_report_->load_pose(db_session, struct_id, pose);
		}
		else{
			tr << "Residues list size " << tag_residues_[tag].size() << std::endl;
			protein_silent_report_->load_pose(db_session, struct_id, tag_residues_[tag], pose);
		}

	} else {
		tr.Debug << "filling pose from saved copy (input tag = " << tag << ")" << endl;
		pose = *(job->inner_job()->get_pose());
	}

	// TODO: Move to pose.clear()
	if (is_symmetric(pose)) make_asymmetric_pose( pose );


	initialize_disulfide_bonds(pose);

}

/// @details this function determines what jobs exist
void protocols::features::DatabaseJobInputter::fill_jobs( protocols::jd2::Jobs & jobs ){
	tr.Debug << "DatabaseJobInputter::fill_jobs" << std::endl;
	jobs.clear(); //should already be empty anyway

	Size const nstruct(get_nstruct());

	if(!tag_structures_.size()){
		sessionOP db_session(basic::database::get_db_session(database_fname_));

		result res;
		while(true)
		{
			try
			{
				res = (*db_session) << "SELECT struct_id FROM structures;";
				break;
			}catch(cppdb::cppdb_error &)
			{
				#ifndef WIN32
					usleep(10);
				#endif
				continue;
			}
		}
		while(res.next()){
			boost::uuids::uuid struct_id;
			res >> struct_id;
			tag_structures_[to_string(struct_id)]=struct_id;
		}

	}

	vector1< protocols::jd2::InnerJobOP > inner_jobs;
	//save list of all inner_jobs first... this allows better sampling
	//of jobs in case of unfinished runs:
	// input1_0001
	// input2_0001
	// ...
	// inputn_0001
	// input1_0002
	// input2_0002
	// ....
	tr.Debug << "reserve memory for InnerJob List " << tag_structures_.size() << endl;
	inner_jobs.reserve( tag_structures_.size() );
	tr.Debug
		<< "fill list with " << tag_structures_.size()
		<< " InnerJob Objects" << endl;

	for(std::map<std::string, boost::uuids::uuid>::const_iterator iter=tag_structures_.begin(); iter!=tag_structures_.end(); ++iter){
		inner_jobs.push_back(new protocols::jd2::InnerJob(iter->first, nstruct));
	} 

	tr.Debug
		<< "reserve list for " << inner_jobs.size() * nstruct
		<< " Job Objects" << endl;

	jobs.reserve(inner_jobs.size() * nstruct);

	tr.Debug << "fill job list with... " << endl;
	for ( Size index = 1; index <= nstruct; ++index ) {
		foreach(protocols::jd2::InnerJobOP ijob, inner_jobs){
			jobs.push_back(new protocols::jd2::Job(ijob, index));
			tr.Trace
				<< "pushing " << ijob->input_tag() << " nstruct index " << index	<< std::endl;
		}
	}
}

/// @brief Return the type of input source that the
///  DatabaseJobInputter is currently using.
/// @return Always <em>DATABASE</em>.
protocols::jd2::JobInputterInputSource::Enum DatabaseJobInputter::input_source() const {
	return protocols::jd2::JobInputterInputSource::DATABASE;
}

//CREATOR SECTION
std::string
DatabaseJobInputterCreator::keyname() const
{
	return "DatabaseJobInputter";
}

protocols::jd2::JobInputterOP
DatabaseJobInputterCreator::create_JobInputter() const {
	return new DatabaseJobInputter;
}

} // namespace features
} // namespace protocols
