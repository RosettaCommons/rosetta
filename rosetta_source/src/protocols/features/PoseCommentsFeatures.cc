// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/PoseCommentsFeatures.cc
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/PoseCommentsFeatures.hh>

// Project Headers
#include <core/pose/util.hh>

//External
#include <boost/uuid/uuid.hpp>

// Platform Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>


// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <string>
#include <map>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::map;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::pose::get_all_comments;
using core::pose::add_comment;
using core::conformation::Residue;
using core::chemical::num_canonical_aas;
using basic::database::table_exists;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;

string
PoseCommentsFeatures::type_name() const { return "PoseCommentsFeatures"; }

string
PoseCommentsFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS pose_comments (\n"
			"	struct_id BLOB,\n"
			"	comment_key TEXT,\n"
			"	value TEXT,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, comment_key));";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS pose_comments (\n"
			"	struct_id BINARY(36),\n"
			"	comment_key VARCHAR(255),\n"
			"	value TEXT,\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id),\n"
			"	PRIMARY KEY(struct_id, comment_key));";
	}else
	{
		return "";
	}

}

utility::vector1<std::string>
PoseCommentsFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}


Size
PoseCommentsFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & /*relevant_residues*/,
	boost::uuids::uuid struct_id,
	sessionOP db_session
){

	typedef map< string, string >::value_type kv_pair;
	//cppdb::transaction transact_guard(*db_session);
	std::string statement_string = "INSERT INTO pose_comments VALUES (?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	foreach(kv_pair const & kv, get_all_comments(pose)){
		stmt.bind(1,struct_id);
		stmt.bind(2,kv.first);
		stmt.bind(3,kv.second);
		basic::database::safely_write_to_database(stmt);
	}
	//transact_guard.commit();
	return 0;
}

void PoseCommentsFeatures::delete_record(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session
) {

	std::string statement_string = "DELETE FROM pose_comments where struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(stmt);

}

void
PoseCommentsFeatures::load_into_pose(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose){

	if(!table_exists(db_session, "pose_comments")) return;


	std::string statement_string =
		"SELECT\n"
		"	comment_key,\n"
		"	value\n"
		"FROM\n"
		"	pose_comments\n"
		"WHERE\n"
		"	pose_comments.struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(stmt));

	while(res.next()){
		string key, value;
		res >> key >> value;
		add_comment(pose, key, value);
	}
}

} //namesapce
} //namespace
