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
#include <basic/database/sql_utils.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

// C++
#include <string>
#include <sstream>
#include <iostream>

//Basic Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.features.StructureFeatures");

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using basic::database::safely_prepare_statement;
using boost::hash_value;
using core::Size;
using core::io::silent::BinaryProteinSilentStruct;
using core::pose::Pose;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StructureFeatures::StructureFeatures(){}

StructureFeatures::StructureFeatures( StructureFeatures const & src) :
	FeaturesReporter()
{}

StructureFeatures::~StructureFeatures(){}

string
StructureFeatures::type_name() const { return "StructureFeatures"; }

string
StructureFeatures::schema() const {
	using namespace basic::database::schema_generator;
	
    //Don't autoincrement the struct_id because it is a UUID generated here
    Column struct_id("struct_id",DbUUID(), false /*not null*/, false /*don't autoincrement*/);
    Column batch_id("batch_id",DbInteger());
    Column tag("tag", DbText(255));
    Column input_tag("input_tag", DbText());
    
    /***structures***/
    Schema structures("structures", PrimaryKey(struct_id));
    
    structures.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true /*defer*/));
    structures.add_column( tag );
    structures.add_column( input_tag );
    
    /***sampled_structures***/
    Schema sampled_structures("sampled_structures");
    sampled_structures.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true /*defer*/));
    
    sampled_structures.add_column( tag );
    sampled_structures.add_column( input_tag );
    
    utility::vector1<Column> unique_cols;
    unique_cols.push_back(tag);
    unique_cols.push_back(batch_id);
    sampled_structures.add_constraint(new UniqueConstraint(unique_cols));
    
    return structures.print() + sampled_structures.print();
}

utility::vector1<std::string>
StructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}


//@details missing struct_id and input/output tags
boost::uuids::uuid
StructureFeatures::report_features(
	vector1< bool > const & relevant_residues,
	Size batch_id,
	sessionOP db_session
){
	string const output_tag(protocols::jd2::JobDistributor::get_instance()->current_output_name());
	string const input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
    boost::uuids::uuid struct_id(
		report_features(relevant_residues, batch_id,
			db_session, output_tag, input_tag));
	return struct_id;
}

//@details missing struct_id and input/output tags
boost::uuids::uuid
StructureFeatures::report_features(
	vector1< bool > const & relevant_residues,
	Size batch_id,
	sessionOP db_session,
	string const & tag,
	string const & input_tag
){
	string statement_string = "INSERT INTO structures (struct_id, batch_id, tag, input_tag) VALUES (?,?,?,?);";
	statement structure_stmt(safely_prepare_statement(statement_string,db_session));

    boost::uuids::uuid struct_id = boost::uuids::random_generator()();
//    boost::uuids::uuid const struct_id(to_string(struct_uuid));
    
//    std::stringstream bytes_stream;
//    std::copy(struct_uuid.begin(), struct_uuid.end(),std::ostream_iterator<char>(bytes_stream, ""));    
//    
//    structure_stmt.bind(1, bytes_stream);
    
    structure_stmt.bind(1, struct_id);
    structure_stmt.bind(2, batch_id);
    structure_stmt.bind(3, tag);
    structure_stmt.bind(4, input_tag);
    
    basic::database::safely_write_to_database(structure_stmt);
    
    return struct_id;
    
//	BinaryProteinSilentStruct silent_struct(pose, "");
//	stringstream pose_string;
//	silent_struct.print_conformation(pose_string);
//	boost::uuids::uuid const struct_id = hash_value(pose_string.str());
//        
//    //Check to see if we've reported this structure before    
//    std::string select_string =
//    "SELECT *\n"
//    "FROM\n"
//    "	structures\n"
//    "WHERE\n"
//    "   struct_id = ?;";
//    cppdb::statement select_stmt(basic::database::safely_prepare_statement(select_string,db_session));
//    select_stmt.bind(1, struct_id);
//    
//    cppdb::result res(basic::database::safely_read_from_database(select_stmt));
//    if(!res.next()) {
//        TR << "No existing structure found, adding the new one" << endl;
//        
//        structure_stmt.bind(1, struct_id);
//        structure_stmt.bind(2, tag);
//        structure_stmt.bind(3, input_tag);
//        
//        TR << "struct id: " << struct_id << "\nbatch_id: " << batch_id << "\ntag: " << tag << "\ninputtag: " << input_tag << endl;
//        basic::database::safely_write_to_database(structure_stmt);
//    }
//    
//    std::string batch_structures_string = "INSERT INTO batch_structures (struct_id, batch_id) VALUES (?,?);";
//    statement batch_structures_stmt(safely_prepare_statement(batch_structures_string, db_session));
//    batch_structures_stmt.bind(1, struct_id);
//    batch_structures_stmt.bind(2, batch_id);
//    basic::database::safely_write_to_database(batch_structures_stmt);
}

void StructureFeatures::mark_structure_as_sampled(
	core::Size batch_id,
	std::string const & tag,
	std::string const & input_tag,
	utility::sql_database::sessionOP db_session
){
	std::string insert_deleted_structure_string = "INSERT INTO sampled_structures (batch_id, tag, input_tag) VALUES (?,?,?);";
	statement insert_statement(safely_prepare_statement(insert_deleted_structure_string,db_session));
	insert_statement.bind(1,batch_id);
	insert_statement.bind(2,tag);
	insert_statement.bind(3,input_tag);
	basic::database::safely_write_to_database(insert_statement);
}

void StructureFeatures::delete_record(
	boost::uuids::uuid struct_id,
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
	boost::uuids::uuid struct_id,
	Pose & pose
){
	load_tag(db_session, struct_id, pose);
}

void
StructureFeatures::load_tag(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
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

boost::uuids::uuid
StructureFeatures::get_struct_id(
	sessionOP db_session,
	string const & tag,
	core::Size const & protocol_id
){

	std::string statement_string =
		"SELECT\n"
		"	structures.struct_id\n"
		"FROM\n"
		"	protocols\n"
        "JOIN batches ON\n"
        "   protocols.protocol_id = batches.protocol_id\n"
        "JOIN structures ON\n"
        "   batches.batch_id = structures.batch_id\n"
		"WHERE\n"
		"	structures.tag=? AND protocols.protocol_id=?;";

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,tag);
	stmt.bind(2,protocol_id);

	result res(basic::database::safely_read_from_database(stmt));
	if(!res.next()){
		stringstream error_message;
		error_message << "Unable to locate structure with tag '"<<tag<<"'."<<endl;
		utility_exit_with_message(error_message.str());
	}
    boost::uuids::uuid struct_id;
	res >> struct_id;
	return struct_id;
}

} // namesapce
} // namespace
