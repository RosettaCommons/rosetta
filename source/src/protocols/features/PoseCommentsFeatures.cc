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

// Platform Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>


#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <string>
#include <map>
#include <sstream>

namespace protocols {
namespace features {

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
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

string
PoseCommentsFeatures::type_name() const { return "PoseCommentsFeatures"; }

void
PoseCommentsFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	//******pose_comments******//
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ), false);
	Column comment_key("comment_key", DbDataTypeOP( new DbTextKey() ), false);
	Column value("value", DbDataTypeOP( new DbText() ), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(comment_key);

	Schema pose_comments("pose_comments", PrimaryKey(pkey_cols));
	pose_comments.add_column(struct_id);
	pose_comments.add_column(comment_key);
	pose_comments.add_column(value);

	pose_comments.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));
	pose_comments.write(db_session);

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
	StructureID struct_id,
	sessionOP db_session
){

	InsertGenerator pose_comments_insert("pose_comments");
	pose_comments_insert.add_column("struct_id");
	pose_comments_insert.add_column("comment_key");
	pose_comments_insert.add_column("value");


	typedef map< string, string >::value_type kv_pair;

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );

	BOOST_FOREACH ( kv_pair const & kv, get_all_comments(pose) ) {

		RowDataBaseOP comment_key_data( new RowData<string>("comment_key",kv.first) );
		RowDataBaseOP value_data( new RowData<string>("value",kv.second) );

		pose_comments_insert.add_row(utility::tools::make_vector(struct_id_data,comment_key_data,value_data));
	}

	pose_comments_insert.write_to_database(db_session);

	return 0;
}

void PoseCommentsFeatures::delete_record(
	StructureID struct_id,
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
	StructureID struct_id,
	Pose & pose){

	if ( !table_exists(db_session, "pose_comments") ) return;


	std::string statement_string =
		"SELECT\n"
		"\tcomment_key,\n"
		"\tvalue\n"
		"FROM\n"
		"\tpose_comments\n"
		"WHERE\n"
		"\tpose_comments.struct_id = ?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(stmt));

	while ( res.next() ) {
		string key, value;
		res >> key >> value;
		add_comment(pose, key, value);
	}
}

} //namesapce
} //namespace
