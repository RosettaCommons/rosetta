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
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <numeric/random/random.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#ifndef __native_client__
#include <boost/lexical_cast.hpp>
#endif

// C++
#include <string>
#include <sstream>
#include <iostream>
#include <limits>

//Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "protocols.features.StructureFeatures" );

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using basic::database::safely_prepare_statement;
using boost::hash_value;
using core::Size;
using core::io::silent::BinarySilentStruct;
using core::pose::Pose;
using core::pose::tag_from_pose;
using core::pose::tag_into_pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowData;
using basic::database::insert_statement_generator::RowDataBaseOP;

StructureFeatures::StructureFeatures(){}

StructureFeatures::StructureFeatures( StructureFeatures const & ) :
	FeaturesReporter()
{}

StructureFeatures::~StructureFeatures(){}

string
StructureFeatures::type_name() const { return "StructureFeatures"; }

void
StructureFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	TR.Debug << "Writing StructureFeatures schema." << std::endl;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false /*not null*/, true);

	if (db_session->is_db_partitioned())
	{
		// If database is partitioned struct_id prefix (32 high bits in structure id) should be
		// set to the database partition identifier.
		//
		// Set autoincrement base value to bit shifted partition id.
		runtime_assert(db_session->get_db_partition() >= 0 && db_session->get_db_partition() < (platform::SSize) std::numeric_limits<uint32_t>::max());

		StructureID structure_prefix = db_session->get_db_partition();
		structure_prefix = structure_prefix << 32;

		TR.Debug << "Setting struct_id autoincrement prefix for partitioned DB. Partition: " << db_session->get_db_partition() << " Prefix: " << structure_prefix << std::endl;

		struct_id = Column("struct_id", DbDataTypeOP( new DbBigInt() ), false /*not null*/, true, structure_prefix + 1);
	}

	Column batch_id("batch_id", DbDataTypeOP( new DbInteger() ));
	Column tag("tag", DbDataTypeOP( new DbText(255) ));
	Column input_tag("input_tag", DbDataTypeOP( new DbText() ));

	/***structures***/
	Schema structures("structures", PrimaryKey(struct_id));

	structures.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true /*defer*/));
	structures.add_column( tag );
	structures.add_column( input_tag );

	structures.write(db_session);

	/***sampled_structures***/
	Schema sampled_structures("sampled_structures");
	sampled_structures.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true /*defer*/));

	sampled_structures.add_column( tag );
	sampled_structures.add_column( input_tag );

	utility::vector1<Column> unique_cols;
	unique_cols.push_back(tag);
	unique_cols.push_back(batch_id);
	sampled_structures.add_constraint(ConstraintOP( new UniqueConstraint(unique_cols) ));

	sampled_structures.write(db_session);
}

utility::vector1<std::string>
StructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}

//@details missing struct_id
StructureID
StructureFeatures::report_features(
	Size batch_id,
	sessionOP db_session,
	string const & tag,
	string const & input_tag
){
	InsertGenerator structures_insert("structures");
	structures_insert.add_column("batch_id");
	structures_insert.add_column("tag");
	structures_insert.add_column("input_tag");

	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id",batch_id) );
	RowDataBaseOP tag_data( new RowData<string>("tag",tag) );
	RowDataBaseOP input_tag_data( new RowData<string>("input_tag",input_tag) );

	structures_insert.add_row(utility::tools::make_vector(batch_id_data,tag_data,input_tag_data));

	long long int structure_sequence_id;
	structures_insert.write_to_database(db_session, structure_sequence_id, "structures_struct_id_seq");

	StructureID inserted_struct_id(structure_sequence_id);

	return inserted_struct_id;
}

void StructureFeatures::mark_structure_as_sampled(
	core::Size batch_id,
	std::string const & tag,
	std::string const & input_tag,
	utility::sql_database::sessionOP db_session
){

	InsertGenerator sampled_insert("sampled_structures");
	sampled_insert.add_column("batch_id");
	sampled_insert.add_column("tag");
	sampled_insert.add_column("input_tag");

	RowDataBaseOP batch_id_data( new RowData<Size>("batch_id",batch_id) );
	RowDataBaseOP tag_data( new RowData<string>("tag",tag) );
	RowDataBaseOP input_tag_data( new RowData<string>("input_tag",input_tag) );

	sampled_insert.add_row(utility::tools::make_vector(batch_id_data,tag_data,input_tag_data));
	sampled_insert.write_to_database(db_session);
}

void StructureFeatures::delete_record(
	StructureID struct_id,
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
	StructureID struct_id,
	Pose & pose
){
	load_tag(db_session, struct_id, pose);
}

void
StructureFeatures::load_tag(
	sessionOP db_session,
	StructureID struct_id,
	Pose & pose) {

#ifndef __native_client__
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
#endif

}

StructureID
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
		"	protocols.protocol_id = batches.protocol_id\n"
		"JOIN structures ON\n"
		"	batches.batch_id = structures.batch_id\n"
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
	StructureID struct_id;
	res >> struct_id;
	return struct_id;
}

} // namesapce
} // namespace
