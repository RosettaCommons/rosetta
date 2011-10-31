// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProtocolFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/StructureFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>

// C++
#include <string>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::pose::Pose;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StructureFeatures::StructureFeatures(){}

StructureFeatures::StructureFeatures( StructureFeatures const & ) :
	FeaturesReporter()
{}

StructureFeatures::~StructureFeatures(){}

string
StructureFeatures::type_name() const { return "StructureFeatures"; }

string
StructureFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);
	if(db_mode == "sqlite3")
		{
		return
			"CREATE TABLE IF NOT EXISTS structures (\n"
			"	struct_id INTEGER PRIMARY KEY AUTOINCREMENT,\n"
			"	protocol_id INTEGER,\n"
			"	tag TEXT,\n"
			"	input_tag TEXT,\n"
			"	UNIQUE (protocol_id, tag)"
			"	FOREIGN KEY (protocol_id)\n"
			"		REFERENCES protocols (protocol_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED);";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS structures (\n"
			"	struct_id INTEGER PRIMARY KEY AUTO_INCREMENT,\n"
			"	protocol_id INTEGER REFERENCES protocols(protocol_id),\n"
			"	tag VARCHAR(255),\n"
			"	input_tag VARCHAR(255),\n"
			"	UNIQUE (protocol_id, tag));";
			//"	FOREIGN KEY (protocol_id) REFERENCES protocols (protocol_id));";
	}else
	{
		return "";
	}
}

Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size protocol_id,
	sessionOP db_session
){
	std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
	return report_features(
		pose, relevant_residues, protocol_id, db_session, find_tag(pose),input_tag);
}

Size
StructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	Size protocol_id,
	sessionOP db_session
){
	std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
	return report_features(
		pose, relevant_residues, struct_id, protocol_id, db_session, find_tag(pose),input_tag);
}


Size
StructureFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	Size protocol_id,
	sessionOP db_session,
	string const & tag,
	string const & input_tag
){
	std::string statement_string = "INSERT INTO structures VALUES (?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	stmt.bind_null(1);
	stmt.bind(2,protocol_id);
	stmt.bind(3,tag);
	stmt.bind(4,input_tag);
	basic::database::safely_write_to_database(stmt);

	return stmt.last_insert_id();
}
//%TODO you stopped here
Size
StructureFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	Size struct_id,
	Size protocol_id,
	sessionOP db_session,
	string const & tag,
	string const & input_tag
){

	std::string statement_string = "INSERT INTO structures VALUES (?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	stmt.bind(1,struct_id);
	stmt.bind(2,protocol_id);
	stmt.bind(3,tag);
	stmt.bind(4,input_tag);
	basic::database::safely_write_to_database(stmt);

	return stmt.last_insert_id();
}
void StructureFeatures::delete_record(
	core::Size struct_id,
	utility::sql_database::sessionOP db_session
){

	std::string statement_string = "DELETE FROM structures WHERE struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(stmt);

}

void
StructureFeatures::load_into_pose(
	sessionOP db_session,
	Size struct_id,
	Pose & pose
){
	load_tag(db_session, struct_id, pose);
}

void
StructureFeatures::load_tag(
	sessionOP db_session,
	Size struct_id,
	Pose & pose) {

	std::string statement_string =
		"SELECT\n"
		"	tag\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.struct_id=?";

	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session))	;
	stmt.bind(1,struct_id);


	result res(basic::database::safely_read_from_database(stmt));
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with struct_id '"
			<< struct_id << "'." << endl;
		utility_exit_with_message(error_message.str());
	}
	string tag;
	res >> tag;
	tag_into_pose(pose,tag);
}

Size
StructureFeatures::get_struct_id(
	sessionOP db_session,
	string const & tag
){

	std::string statement_string =
		"SELECT\n"
		"	struct_id\n"
		"FROM\n"
		"	structures\n"
		"WHERE\n"
		"	structures.tag=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,tag);

	result res(basic::database::safely_read_from_database(stmt));
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with tag '"<<tag<<"'."<<endl;
		utility_exit_with_message(error_message.str());
	}
	Size struct_id;
	res >> struct_id;
	return struct_id;
}

} // namesapce
} // namespace
